'''
This script defines a model class that can be used to compute the
overturning circulation in the north of the basin, given buoyancy profiles
in the basin and in the northern deepwater formation region.
 
The model assumes a thermal-wind based equation for the overturning circulation
as in Nikurashin and Vallis (2012):
d_{zz}(\Psi_{N}) = f^{-1} (b_{N} - b_{basin})

This equation is solved subject to the boundary conditions that
\Psi(0) = \Psi(-H) = 0 
'''

import numpy as np
from scipy import integrate

class Model_PsiNA(object):
    def __init__(
            self,
            f=1.2e-4,         # Coriolis parameter (input)
            z=None,           # grid (input)
            sol_init=None,    # Initial conditions for ODE solver (input)
            b_basin=None,     # Buoyancy in the basin (input, output)
            b_N=0.,           # Buoyancy in the deep water formation region (input, output)
            Psi= None,      # Streamfunction (output) 
    ):
 
        self.f = f
        # initialize grid:
        if isinstance(z,np.ndarray):
            self.z = z
            nz = np.size(z) 
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
                        
        self.b_basin=self.make_func(b_basin,'b_basin',self.z)
        self.b_N=self.make_func(b_N,'b_N',self.z)
                 
        # Set initial conditions for BVP solver
        if sol_init:
            self.sol_init = sol_init
        else:
            self.sol_init = np.zeros((2, nz))    
    # end of init
     
    
    def make_func(self,myst,name,zin):
    # turn array or float into callable function (if needed)    
        if callable(myst):
            return myst
        elif isinstance(myst,np.ndarray):
            def funfun(z): return np.interp(z,zin,myst)
            return funfun
        elif isinstance(myst,float):
            def funfun(z): return myst +0*z
            return funfun
        else:
            raise TypeError(name,'needs to be either function, numpy array, or float') 
       
    def bc(self, ya, yb):
        #return the boundary conditions
        return np.array([ya[0], yb[0]])

    def ode(self, z, y):
        #return the equation to be solved 
        return np.vstack((y[1], 1./self.f*(self.b_N(z)-self.b_basin(z))))
                         
    def solve(self):
        #Solve the boundary value problem
        res = integrate.solve_bvp(self.ode, self.bc, self.z, self.sol_init)
        # interpolate solution for overturning circulation onto original grid (and change units to SV)
        self.Psi = res.sol(self.z)[0, :] / 1e6  
       
    
    def update(self,b_basin=None,b_N=None):
        # update buoyancy profiles     
        if b_basin is not None:
            self.b_basin=self.make_func(b_basin,'b_basin',self.z)
        if b_N is not None:
            self.b_N=self.make_func(b_N,'b_N',self.z)
        