import numpy as np
from scipy import integrate

class Model(object):
    def __init__(
            self,
            f=1.2e-4,
            b_s=0.025,
            kappa_const=6e-5,
            B_int=3e3,
            nz=100,
            A=7e13,
            sol_init=None,
            kappa=None,
            dkappa_dz=None,
            psi_so=None,
    ):
   
        self.f = f
        self.A = A
        self.zi=np.asarray(np.linspace(-1, 0, nz))
        if kappa and dkappa_dz:  
            self.kappa= lambda z,H: kappa(z*H)/ (H**2 * self.f) # non-dimensionalize (incl. norm. of vertical coordinate)
            self.dkappa_dz=lambda z,H: dkappa_dz(z*H) / (H * self.f)
        elif kappa or dkappa_dz:
           raise ValueError('If you want kappa profile you need to provide both kappa and dkappa_dz')
        else:
            self.kappa= lambda z,H: kappa_const/ (H**2 * self.f)
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
        H = res.p[0]
        z = res.x * H
        y = res.y
        y[0, :] = y[0, :] * self.f * H**3 / 1e6
        y[2, :] = -y[2, :] * self.f**2 * H
        return dict({'z': z, 'psi': y, 'H': H})
