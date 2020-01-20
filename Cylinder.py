import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize

class Cylinder:
    def __init__(self, C=np.array([0, 0, 0]), W=np.array([0, 0, 1]), r=1, h=np.inf):
        self.C = C
        self.W = W
        self.r = r
        self.h = h

    @property
    def W(self):
        return self._W

    @property
    def _P(self):
        return np.eye(3) - np.outer(self._W, self._W.T)

    @W.setter
    def W(self, W):
        self._W = W/np.linalg.norm(W)
        # Choose any right-hand orthonormal set so that: UxV=W, VxW=U, and WxU=V.
        ns = null_space(np.row_stack((self._W, np.zeros(3), np.zeros(3))))
        self._U = ns[:,0]
        self._V = ns[:,1]
        
class CylinderFit(Cylinder):
    def __init__(self, X=None, *args, **kwargs):
        self.X = X
        super().__init__(*args, **kwargs)

    @property
    def X(self):
        if self._X is not None:
            return self._X + self._X_mean[:,np.newaxis]
        else:
            return None

    @property
    def n(self):
        return np.shape(self._X)[1] # X consists of column vectors of coordinate points (3xn).

    @X.setter
    def X(self, X):
        if X is not None:
            self._X_mean = np.mean(X, axis=1)
            self._X = X - self._X_mean[:,np.newaxis]
        else:
            self._X = None
            self._X_mean = None

    @property
    def _A(self):
        return (1/self.n)*self._P@self._X@self._X.T@self._P

    @property
    def _Ahat(self):
        S = self.skew_symmetric_matrix(self.W)
        return S@self._A@S.T # Eq. 128 in [1]
    
    def _G(self):
        XPX = np.einsum('ij,ji->i', self._X.T@self._P, self._X)
        G = (1/self.n)*np.sum(np.square(XPX - (1/self.n)*np.sum(XPX) - 2*self._X.T@(self._Ahat/np.trace(self._Ahat@self._A))@((1/self.n)*np.sum(XPX*self._X, axis=1)))) # Eq. 132 in [1]
        return G

    def fit(self, angles=[(0, 0), (np.pi/2, 0), (np.pi/2, np.pi/2)]):
        best_fit = None
        best_score = np.inf
        for ang in angles:
            fitted = minimize(lambda x: self.minimization_handler(x), ang, method='Powell', tol=1e-6)
            if fitted.fun < best_score:
                best_score = fitted.fun
                best_fit = fitted
        XPX = np.einsum('ij,ji->i', self._X.T@self._P, self._X)
        self.W = self.ang_to_direction((best_fit.x[0], best_fit.x[1]))
        self.C = (self._Ahat/np.trace(self._Ahat@self._A))@((1/self.n)*np.sum(XPX*self._X, axis=1)) + self._X_mean # Eq. 130 in [1]
        self.r = np.sqrt((1/self.n)*np.einsum('ij,ji->', (self._P@self.C[:,np.newaxis] - self._P@self._X).T, (self._P@self.C[:,np.newaxis] - self._P@self._X))) # Eq. 112 in [1]
        XdotW = np.dot((self.X - self.C[:,np.newaxis]).T, self.W)
        self.h = np.max(XdotW) - np.min(XdotW)
        print(self.h)

    def minimization_handler(self, ang):
        self.W = self.ang_to_direction(ang)
        return self._G()

    @staticmethod
    def skew_symmetric_matrix(vec):
        S = np.array([[0, -vec[2], vec[1]], [vec[2], 0, -vec[0]], [-vec[1], vec[0], 0]])
        return S

    @staticmethod
    def ang_to_direction(ang):
        return np.array([np.cos(ang[1])*np.sin(ang[0]), np.sin(ang[1])*np.sin(ang[0]), np.cos(ang[0])])