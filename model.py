import math
import numpy as np
from scipy import integrate

class Model(object):
    def __init__(
            self,
            f=1.2e-4,
            b_s=0.025,
            kappa=6e-5,
            kappa_back=1e-5,
            kappa_s=3e-5,
            kappa_4k=3e-4,
            H_max_so=2000,
            B_int=3e3,
            psi_so_max=4.,
            nz=100,
            A=7e13,
            sol_init=None,
            diff_type='constant',
    ):
        self.f = f
        self.A = A
        self.zi=np.asarray(np.linspace(-1, 0, nz))

        if diff_type == 'constant':
            self.kappa = lambda z, H: kappa / (f * H**2)
            self.dkappa_dz = lambda z, H: 0
        else:
            self.kappa = lambda z, H: (kappa_back + kappa_s*np.exp(z*H/100)+kappa_4k*np.exp(-z*H/1000 - 4)) / (H**2 * f) # nondim kappa
            self.dkappa_dz = lambda z, H: (kappa_s/100*np.exp(z*H/100)-kappa_4k/1000*np.exp(-z*H/1000-4)) / (H * f) #nondimensional d(kappa)/dz  
        self.b = -b_s / f**2
        self.B_int = B_int
        self.psi_so_max = psi_so_max
        self.H_max_so = H_max_so
        if sol_init:
            self.sol_init = sol_init
        else:
            sol_init = np.zeros((4, nz))
            sol_init[0,:] = np.ones((nz))
            sol_init[2,:] = -0.1 * np.ones((nz))
            sol_init[3,:] = -self.bz(H_max_so + 200) * np.ones((nz))
            self.sol_init = sol_init

    def alpha(self, z, H):
        return H**2 / (self.A * self.kappa(0, H))

    def bz(self, H): 
        return self.B_int / (self.f**3 * H**2 * self.A * self.kappa(0, H))
    
    def bc(self, ya, yb, p):
        return np.array([ya[0], yb[0], ya[1], ya[3] + self.bz(p[0]), yb[2] - self.b/p[0]])

    def psi_so(self, z, H):
        pre = self.psi_so_max * 1e6 / (self.f * H**3)
        return pre * np.sin([-np.pi * max(H * x, -self.H_max_so) / self.H_max_so for x in z])**2

    def ode(self, z, y, p):
        H = p[0]
        return np.vstack((y[1], y[2], y[3], self.alpha(z, H) * y[3] * (y[0] - self.psi_so(z, H) - self.A * self.dkappa_dz(z, H)/(H**2))))

    def solve(self):
        res = integrate.solve_bvp(self.ode, self.bc, self.zi, self.sol_init, p=[self.H_max_so - 200])
        H = res.p[0]
        z = res.x * H
        y = res.y
        y[0, :] = y[0, :] * self.f * H**3 / 1e6
        y[2, :] = -y[2, :] * self.f**2 * H
        return dict({'z': z, 'psi': y, 'H': H})
