import numpy as np
import itertools

from Cylinder import CylinderFit

class Nanotube:

    ELEMENT_MASSES = {'H': 1.0079, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453}
    BOND_LENGTHS = {'H': 109, 'C': 154, 'N': 178.5, 'O': 179, 'F': 134, 'P': 187, 'S': 218, 'Cl': 176}
    UNIT_CONVERSIONS = {'pm': 1., 'angstrom': 100., 'bohr': 52.9177}

    def __init__(self, atoms=None, xyz=None, unit='angstrom'):
        self.atoms = atoms
        self.geom = xyz
        self.unit = unit
        self.func_atoms = []
        self.func_geom = np.empty(0)

    def fromXYZ(self, xyz_path, unit='angstrom'):
        self.unit = unit
        self.atoms, self.geom = self._importXYZ(xyz_path, self.UNIT_CONVERSIONS[unit])

    def toXYZ(self, xyz_path, out_unit=None):
        if out_unit is None:
            out_unit = self.unit
        with open(xyz_path, 'w') as out_file:
            out_file.write('{0:d}\n\n'.format(self.geom.shape[0] + self.func_geom.shape[0]))
            for i in range(self.geom.shape[0]):
                out_file.write('{0:<2s}  {1:>14.5f} {2:>14.5f} {3:14.5f}\n'.format(self.atoms[i], self.geom[i,0]/self.UNIT_CONVERSIONS[out_unit], self.geom[i,1]/self.UNIT_CONVERSIONS[out_unit], self.geom[i,2]/self.UNIT_CONVERSIONS[out_unit]))
            for i in range(self.func_geom.shape[0]):
                out_file.write('{0:<2s}  {1:>14.5f} {2:>14.5f} {3:14.5f}\n'.format(self.func_atoms[i], self.func_geom[i,0]/self.UNIT_CONVERSIONS[out_unit], self.func_geom[i,1]/self.UNIT_CONVERSIONS[out_unit], self.func_geom[i,2]/self.UNIT_CONVERSIONS[out_unit]))

    def functionalize(self, xyz_path, wt_frac, method='random', bond_length=0., normal=np.array([0, 0, 1]), origin=np.array([0, 0, 0])):
        func_atoms, func_geom = self._importXYZ(xyz_path, self.UNIT_CONVERSIONS[self.unit])
        func_mmass = self._molarMass(func_atoms)
        cnt_mmass = self._molarMass(self.atoms)
        n = int(np.round(wt_frac*cnt_mmass/((1 - wt_frac)*func_mmass)))
        if np.isscalar(origin):
            if bond_length is None:
                bond_length = self.BOND_LENGTHS[func_atoms[origin]]
            origin = func_geom[origin]
        if isinstance(bond_length, str):
            bond_length = self.BOND_LENGTHS[bond_length]
        # Translate so that the origin is at [0, 0, 0].
        func_geom = func_geom - origin
        # Rotate so that normal is the +z axis.
        carbon_geom = np.array([x for x in self.geom if 'C' in self.atoms])
        self._cylinder = CylinderFit(X=carbon_geom.T)
        self._cylinder.fit()
        print(self._cylinder.C)
        print(self._cylinder.W)
        z_axis = np.array([0, 0, 1])
        rotation_to_z = self._rotateMatrix(normal, z_axis)
        func_geom = np.transpose(rotation_to_z@func_geom.T)
        # Translate so that the new origin is shifted along the normal (z-axis) direction by bond_length.
        func_geom = func_geom + np.array([0, 0, bond_length])
        # Choose the carbon locations that will bond to functional groups.
        carbon_location = self._chooseCarbons(carbon_geom, n, min_dist=2*self.BOND_LENGTHS['C'], method=method)
        # Position the functional groups
        all_func_coord = list(itertools.repeat(func_geom, n))
        all_func_coord = [self._placeFunctionalGroup(coord, carbon_geom[carbon_location[i],:]) for i, coord in enumerate(all_func_coord)]
        self.func_geom = np.concatenate(all_func_coord, 0)
        self.func_atoms = func_atoms*n

    def _chooseCarbons(self, carbon_geom, n, method='random', min_dist=0):
        if method=='random':
            choose_fn = lambda x, n: np.random.choice(x, n, replace=False)
            num_iter = 0
            max_iter = 100
            num_passed = 0
            locs = np.empty(0, dtype='int32')
            while (num_passed!=n):
                locs = np.union1d(locs, choose_fn(np.setdiff1d(np.arange(carbon_geom.shape[0]), locs, assume_unique=True), n-num_passed))
                dist = np.ma.array(self._distanceBetweenAllPoints(carbon_geom[locs,:]), mask=np.triu(np.ones((n, n))))
                is_close = np.any(dist<min_dist, axis=1).data
                locs = locs[~is_close]
                num_passed = np.sum(~is_close)
                num_iter = num_iter + 1
                if num_iter==max_iter:
                    raise Exception('Could not find substituent locations that satisfy min_dist={} pm with method {}.'.format(min_dist, method))
        elif method=='stride':
            locs = np.round(np.linspace(0, carbon_geom.shape[0], num=n))
        elif method=='kmeans':
            raise NotImplementedError
        else:
            raise NotImplementedError
        return locs

    @staticmethod
    def _distanceBetweenAllPoints(points):
        n = points.shape[0]
        return np.reshape(np.sqrt(np.sum((np.repeat(points, n, axis=0) - np.tile(points, (n,1)))**2, axis=1)), (n, n))

    def _placeFunctionalGroup(self, group_geom, pos):
        normal = (pos - self._cylinder.C) - np.dot((pos - self._cylinder.C), (self._cylinder.W))*self._cylinder.W
        z_axis = np.array([0, 0, 1])
        rotation_to_normal = self._rotateMatrix(z_axis, normal)
        group_geom = np.transpose(rotation_to_normal@group_geom.T) + pos
        return group_geom

    def _rotateMatrix(self, a, b):
        a = a/np.linalg.norm(a)
        b = b/np.linalg.norm(b)
        if np.allclose(a, b):
            rot = np.eye(3)
        elif np.allclose(a, -1*b):
            rot = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
        else:
            v = np.cross(a, b)
            V = self._cylinder.skew_symmetric_matrix(v)
            rot = np.eye(3) + V + np.dot(V,V)*((1 - np.dot(a, b))/np.linalg.norm(v)**2)
        return rot

    @classmethod
    def _molarMass(cls, atoms):
        mmass = 0.
        for at in atoms:
            mmass = mmass + cls.ELEMENT_MASSES[at]
        return mmass

    @staticmethod
    def _importXYZ(xyz_path, conversion):
        file_data = np.genfromtxt(xyz_path, dtype=['U2', 'f', 'f', 'f'], skip_header=2, names=['a', 'x', 'y', 'z'])
        atoms = file_data['a'].tolist()
        geom = (np.hstack((file_data['x'][:,np.newaxis], file_data['y'][:,np.newaxis], file_data['z'][:,np.newaxis])))*conversion
        return atoms, geom