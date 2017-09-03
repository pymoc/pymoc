'''
This script defines a model class that can be used to solve a 1D column model
for the AMOC.
 
The model is written in terms of a boundary value problem
 
In dimensional units, the ODE is
d_{zzzz}(\Psi_{N}) = (\kappa A)^{-1}( \Psi_{N} - \Psi_{SO} - A d_z(\kappa)d_{zzz}(\Psi_N)

This equation is solved subject to the boundary conditions:
(1) \Psi_N(0) = 0
(2) \Psi_N(-H) = 0 
(3) d_z\psi_N(-H) = 0
(4) b(0)=-f \partial_{zz} \Psi_N (0) = b_s
(5) d_z b(-H) = -f d_{zzz} \Psi_N (-H) =  (A \kappa(-H))^{-1} B_{int}
Where H is the total depth of the upper cell, which is also solved for

The solution is found by non-dimensionalizing the equations using H and f as length and time scales
The model is then solved between z^*=z/H=0..1
Notice that H then appears as a parameter (to be determined) in the equations. 5 boundary conditions are needed because we are dealing with a 4th order ODE and are also solving for the parameter H
'''

import numpy as np
from scipy import integrate

class Model_Equi(object):
    def __init__(
            self,
            f=1.2e-4,         # Coriolis parameter (input)
            b_s=0.025,        # surface buoyancy (input)
            B_int=3e3,        # integrated downward buoyancy flux at the bottom of the upper cell (input)
            A=7e13,           # Area of the basin (input) 
            nz=100,           # minimum number of vert. layers for numerical solver (input)
            sol_init=None,    # Initial conditions for ODE solver (input)
            kappa=6e-5,       # Diapycnal diffusivity (input; can be const., function, or array on grid given by z)
            dkappa_dz=None,   # Vertical derivative of diffusivity profile (input; function or nothing)
            psi_so=None,      # SO streamfunction (input; function or array on grid given by z)
            z=None,           # vertical grid for I/O (input / output) 
            psi=None,         # streamfunction (output)
            b=None,           # streamfunction (output)
            H=None,           # depth of cell (output)
    ):
 
        self.f = f
        self.A = A
        self.z = z
        self.zi=np.asarray(np.linspace(-1, 0, nz)) # grid for initial conditions for solver        
        
        # Initialize vertical diffusivity profile:
        if callable(kappa): 
            self.kappa= lambda z,H: kappa(z*H)/ (H**2 * self.f) # non-dimensionalize (incl. norm. of vertical coordinate)
            if callable(dkappa_dz):  
                self.dkappa_dz=lambda z,H: dkappa_dz(z*H) / (H * self.f)
            else:
                if not self.check_numpy_version():
                    raise ImportError('You need NumPy version 1.13.0 or later if you want to automatically compute dkappa_dz. Please upgrade your NumPy libary or provide functional form of dkappa_dz.')
                self.dkappa_dz=lambda z,H: np.gradient(kappa(z*H), z*H) / (H * self.f)
        elif isinstance(kappa,np.ndarray):
            self.kappa= lambda z,H: np.interp(z*H,self.z,kappa)/(H**2 * self.f)
            if not self.check_numpy_version():
                raise ImportError('You need NumPy version 1.13.0 or later if you want to automatically compute dkappa_dz. Please upgrade your NumPy libary.')
            dkappa_dz= np.gradient(kappa, z)
            self.dkappa_dz= lambda z,H: np.interp(z*H,self.z,dkappa_dz)/(H * self.f)
        else:
            self.kappa= lambda z,H: kappa /(H**2 * self.f)
            self.dkappa_dz=lambda z,H: 0   
 
        # Initialize Southern Ocean Streamfunction   
        if callable(psi_so):
            self.psi_so=lambda z,H: psi_so(z*H) /(self.f * H**3) # non-dimensionalize (incl. norm. of vertical coordinate)
        elif isinstance(psi_so,np.ndarray):
            self.psi_so=lambda z,H: np.interp(z*H,self.z,psi_so) /(self.f * H**3) 
        else:
            self.psi_so=lambda z,H: 0
            
        # initialize non-dimensional surface buoyancy and abyssal buoyancy flux boundary condition    
        self.b = -b_s / f**2
        self.B_int = B_int
        
        # Set initial conditions for ODE solver
        if sol_init:
            self.sol_init = sol_init
        else:
            sol_init = np.zeros((4, nz))
            sol_init[0,:] = np.ones((nz))
            sol_init[2,:] = 0.0 * np.ones((nz))
            sol_init[3,:] = -self.bz(1500.) * np.ones((nz))
            self.sol_init = sol_init    
    # end of init
        
    
    def check_numpy_version(self):
        # check numpy version (version >= 1.13 needed to automatically compute db/dz)    
        v = [int(i) for i in np.version.version.split('.')]
        if v[0] <= 1 and v[1] < 13:
            return False
        return True
    
    def alpha(self, z, H):
        #return factor on the RHS of ODE
        return H**2 / (self.A * self.kappa(z, H))
    
    def bz(self, H): 
        # return the properly non-dimensionalized stratification at the bottom of the cell
        return self.B_int / (self.f**3 * H**2 * self.A * self.kappa(-1, H))
    
    def bc(self, ya, yb, p):
         #return the boundary conditions for the ODE
        return np.array([ya[0], yb[0], ya[1], ya[3] + self.bz(p[0]), yb[2] - self.b/p[0]])

    def ode(self, z, y, p):
        #return the ODE to be solved 
        H = p[0]
        return np.vstack((y[1], y[2], y[3], self.alpha(z, H) * y[3] * (y[0] - self.psi_so(z, H) - self.A * self.dkappa_dz(z, H)/(H**2))))

    def solve(self):
        #Solve the boundary value problem
        res = integrate.solve_bvp(self.ode, self.bc, self.zi, self.sol_init, p=[1500.])
        self.H = res.p[0]
        # if self.z does not yet exist use mesh from solver:
        if self.z is None:
            self.z   = res.x * self.H
            self.psi = res.y[0, :] * self.f * self.H**3 / 1e6
            self.b   =-res.y[2, :] * self.f**2 * self.H
        # if self.z does exist, compute solution at those points (setting all points below z=-H to NaN):
        else:
            self.psi = res.sol(self.z/self.H)[0, :] * self.f * self.H**3 / 1e6
            self.psi[self.z<-self.H]= np.NaN
            self.b   =-res.sol(self.z/self.H)[2, :] * self.f**2 * self.H
            self.b[self.z<-self.H]= np.NaN
    
