import numpy as np
from scipy import integrate

class Model(object):
    def __init__(
            self,
            f=1.2e-4,
            b_s=0.025,
            B_int=3e3,
            A=7e13,
            nz=100,
            sol_init=None,
            kappa=6e-5,
            dkappa_dz=None,
            psi_so=None,
            z=None,
            psi=None,
            b=None,
            H=None,
    ):
 
        self.f = f
        self.A = A
        self.z = z
        self.zi=np.asarray(np.linspace(-1, 0, nz))
        
        if callable(kappa): 
            self.kappa= lambda z,H: kappa(z*H)/ (H**2 * self.f) # non-dimensionalize (incl. norm. of vertical coordinate)
            if callable(dkappa_dz):  
                self.dkappa_dz=lambda z,H: dkappa_dz(z*H) / (H * self.f)
            else:
                if not self.check_numpy_version():
                    raise ImportError('You need NumPy version 1.13.0 or later if you want to automatically compute dkappa_dz. Please upgrade your NumPy libary or provide functional form of dkappa_dz.')
                self.dkappa_dz=lambda z,H: np.gradient(kappa(z*H), z*H) / (H * self.f)
        else:
            self.kappa= lambda z,H: kappa/ (H**2 * self.f)
            self.dkappa_dz=lambda z,H: 0   
 
        if psi_so:
            self.psi_so=lambda z,H: psi_so(z*H)/ (self.f * H**3) # non-dimensionalize (incl. norm. of vertical coordinate)
        else:
            self.psi_so=lambda z,H: 0
            
        self.b = -b_s / f**2
        self.B_int = B_int
        if sol_init:
            self.sol_init = sol_init
        else:
            sol_init = np.zeros((4, nz))
            sol_init[0,:] = np.ones((nz))
            sol_init[2,:] = -0.1 * np.ones((nz))
            sol_init[3,:] = -self.bz(2000.) * np.ones((nz))
            self.sol_init = sol_init

    def check_numpy_version(self):
        v = [int(i) for i in np.version.version.split('.')]
        if v[0] <= 1 and v[1] < 13:
            return False
        return True
    
    def alpha(self, z, H):
        return H**2 / (self.A * self.kappa(0, H))

    def bz(self, H): 
        return self.B_int / (self.f**3 * H**2 * self.A * self.kappa(0, H))
    
    def bc(self, ya, yb, p):
        return np.array([ya[0], yb[0], ya[1], ya[3] + self.bz(p[0]), yb[2] - self.b/p[0]])

    def ode(self, z, y, p):
        H = p[0]
        return np.vstack((y[1], y[2], y[3], self.alpha(z, H) * y[3] * (y[0] - self.psi_so(z, H) - self.A * self.dkappa_dz(z, H)/(H**2))))

    def solve(self):
        res = integrate.solve_bvp(self.ode, self.bc, self.zi, self.sol_init, p=[2000.])
        self.H = res.p[0]
        # if self.z does not yet exist use mesh from solver
        if self.z is None:
            self.z   = res.x * self.H
            self.psi = res.y[0, :] * self.f * self.H**3 / 1e6
            self.b   =-res.y[2, :] * self.f**2 * self.H
        # if self.z does exist, compute solution at those points (setting all points below z=-H to NaN)
        else:
            self.psi = res.sol(self.z/self.H)[0, :] * self.f * self.H**3 / 1e6
            self.psi[self.z<-self.H]= np.NaN
            self.b   =-res.sol(self.z/self.H)[2, :] * self.f**2 * self.H
            self.b[self.z<-self.H]= np.NaN
 
         