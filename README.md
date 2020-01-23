# nanotube-functionalizer
Adds functional groups to carbon nanotube geometries.

## Idea
There are many tools for generating carbon nanotube (CNT) structures for modelling, but adding functional groups may be tedious. Nanotube Functionalizer automates this process easier by doing the following:
1. Importing CNT and functional group molecular geometries from XYZ files.
2. Fitting a cylinder to the CNT [1].
3. Choosing carbon molecules that will be bonded to functional groups.
4. Calculating the number of substituents for a given weight fraction.
5. Positioning the substituents normal to the CNT surface at the correct bond distance away.
6. Writing an XYZ geometry file of the final, substituted molecule.

## Caveats
* The unfunctionalized CNT geometry should be very clean, i.e. as un-optimized as possible, as this improves the cylinder fit. The XYZ file should only contain the CNT without any other molecules present.

## Prerequisites
* Numpy
* Scipy

## Usage
1. Obtain a nanotube geometry using your favourite software and export in XYZ format. E.g. Avogadro's nanotube builder.
2. Generate a geometry for the functional group of interest and also export in XYZ format. Note where the axes align relative to the functional group.
3. In your script, generate a `Nanotube` object from the XYZ file:

        import Nanotube
        n = Nanotube()
        n.fromXYZ(xyz_path='<PATH_TO_CNT>.xyz', \
                  unit='angstrom')
                  
   Specify the distance unit to be the same as that of the geometry files.
        
4. Add the functional group to the nanotube with the `functionalize()` method:
    
        n.functionalize(xyz_path='<PATH_TO_FUNC>.xyz', \
                        wt_frac=0.1, \
                        method='random', \
                        bond_length=0, \
                        normal=np.array([0, 0, 1]), \
                        origin=np.array([0, 0, 0]))
                
   Where the arguments to `functionalize()` are:
    1. `xyz_path`: path to the functional group geometry file.
    2. `wt_frac`: weight fraction of substituent in the functionalized molecule. Nanotube Functionalizer currently supports groups with H, C, N, O, F, P, S, and Cl atoms. 
    3. `method`: algorithm to choose the carbons bonded to the substituents. Currently, only `random` (randomly chosen) and `stride` (every nth atom) are supported.
    4. `bond_length`: distance from the target carbon to the origin of the functional group. In this example, the origin was already placed 1 C-N bond length away, so `bond_length` can be set to zero.
    5. `normal`: the normal vector for the functional group. In this example, the normal vector is along the +z axis.
    6. `origin`: an offset for the functional group origin.
5. Write the result to an XYZ file:
      
        n.toXYZ('<PATH_TO_OUTPUT>.xyz')

## References
[1] D. Eberly, “Least Squares Fitting of Data by Linear or Quadratic Structures.” Geometric Tools, ch. 7, pp. 38–51, 2019.
