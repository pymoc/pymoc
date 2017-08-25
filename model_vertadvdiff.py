'''
Instances of this class represent 1D representations of buoyancy in a water
column governed by vertical advection and diffusion the vertical velocity
profile is required as an input. The script can either compute the equilibrium
buoyancy profile for a given vertical velocity profile and boundary conditions
or compute the tendency and perform a time-step of given length
(Or so the plan... As of 08/25, the time-stepping is actually not yet implemented!)
For now, BCs have to be buoyancy at the top and bottom
'''

import numpy as np
from scipy import integrate

class Model_VertAdvDiff(object):
    def __init__(
            self,
            z=None,         # grid (input)
            kappa=None,     # diffusivity profile (input)
            dkappa_dz=None, # vertical derivative diffusivity profile (optional input)
            bs=0.025,       # surface buoyancy bound. cond (input)
            bbot=0.0,       # bottom buoyancy boundary condition (input)  
            b=None,         # Buoyancy profile (input, output)
    ):
 
        # initialize grid:
        if isinstance(z,np.ndarray):
            self.z = z
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
                      
        self.kappa=self.make_func(kappa,'kappa')
        if callable(dkappa_dz):  
            self.dkappa_dz=dkappa_dz
        else:
            if not self.check_numpy_version():
                raise ImportError('You need NumPy version 1.13.0 or later if you want to automatically compute dkappa_dz. Please upgrade your NumPy libary or provide functional form of dkappa_dz.')
            self.dkappa_dz=lambda z: np.gradient(self.kappa(z), z)
    
        self.bs=bs
        self.bbot=bbot
         
        if isinstance(b,np.ndarray):
           self.b=b
        elif callable(b):
           self.b=b(self.z)
        elif isinstance(b,float):
           self.b=b+0*self.z 
        else:
           raise TypeError(name,'needs to be either function, numpy array, or float') 
            
        
        if self.check_numpy_version():
           self.bz=np.gradient(self.b, z)
        else:
           self.bz=0.*z   
           
    def check_numpy_version(self):
        # check numpy version (version >= 1.13 needed to automatically compute db/dz)    
        v = [int(i) for i in np.version.version.split('.')]
        if v[0] <= 1 and v[1] < 13:
            return False
        return True
       
    
    def make_func(self,myst,name):
    # turn array or float into callable function (if needed)    
        if callable(myst):
            return myst
        elif isinstance(myst,np.ndarray):
            def funfun(z): return np.interp(z,self.z,myst)
            return funfun
        elif isinstance(myst,float):
            def funfun(z): return myst +0*z
            return funfun
        else:
            raise TypeError(name,'needs to be either function, numpy array, or float') 
     
    def bc(self, ya, yb):
        #return the boundary conditions
        return np.array([ya[0]-self.bbot, yb[0]-self.bs])

    def ode(self, z, y):
        #return the equation to be solved 
        return np.vstack((y[1], (self.w(z)-self.dkappa_dz(z))/self.kappa(z)*y[1]))
                         
    def solve_equi(self,w):
        #Solve for equilibrium solution given vert. vel profile and BCs
        self.w=self.make_func(w,'w')
        sol_init = np.zeros((2, np.size(self.z)));
        sol_init[0,:] = self.b;sol_init[1,:] = self.bz;
        res = integrate.solve_bvp(self.ode, self.bc, self.z, sol_init)
        # interpolate solution for b and db/dz onto original grid
        self.b = res.sol(self.z)[0, :]  
        self.bz = res.sol(self.z)[1, :]  
        
        
        