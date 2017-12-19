'''
This script defines a model class that can be used to compute the
overturning circulation between two columns, given buoyancy profiles
in the columns.
 
The model assumes a thermal-wind based equation for the overturning circulation
as in Nikurashin and Vallis (2012):
d_{zz}(\Psi) = f^{-1} (b_{2} - b_{1})

This equation is solved subject to the boundary conditions that
\Psi(0) = \Psi(-H) = 0 
'''

import numpy as np
from scipy import integrate

class Model_Thermwind(object):
    def __init__(
            self,
            f=1.2e-4,         # Coriolis parameter (input)
            z=None,           # grid (input)
            sol_init=None,    # Initial conditions for ODE solver (input)
            b1=None,     # Buoyancy in the basin (input, output)
            b2=0.,           # Buoyancy in the deep water formation region (input, output)
            Psi= None,      # Streamfunction (output) 
    ):
 
        self.f = f
        # initialize grid:
        if isinstance(z,np.ndarray):
            self.z = z
            nz = np.size(z) 
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
                        
        self.b1=self.make_func(b1,'b1',self.z)
        self.b2=self.make_func(b2,'b2',self.z)
                 
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
       
    def make_array(self,myst,name):
    # turn mysterious object into array(if needed)    
        if isinstance(myst,np.ndarray):
            return myst
        elif callable(myst):
            return myst(self.z)
        elif isinstance(myst,float):
            return myst+0*self.z
        else:
            raise TypeError(name,'needs to be either function, numpy array, or float') 

    
    def bc(self, ya, yb):
        #return the boundary conditions
        return np.array([ya[0], yb[0]])

    def ode(self, z, y):
        #return the equation to be solved 
        return np.vstack((y[1], 1./self.f*(self.b2(z)-self.b1(z))))
                         
    def solve(self):
        #Solve the boundary value problem
        res = integrate.solve_bvp(self.ode, self.bc, self.z, self.sol_init)
        # interpolate solution for overturning circulation onto original grid (and change units to SV)
        self.Psi = res.sol(self.z)[0, :] / 1e6  
       
    def Psib(self,nb=500):
        b1=self.make_array(self.b1,'b1')
        b2=self.make_array(self.b2,'b2')
        bmin=min(np.min(b1),np.min(b2))
        bmax=max(np.max(b1),np.max(b2))
        self.bgrid=np.linspace(bmin,bmax,nb)
        udydz=self.Psi[1:]-self.Psi[:-1]
        b1mid=0.5*(b1[1:]+b1[:-1])
        b2mid=0.5*(b2[1:]+b2[:-1])
        bgridmid=0.5*(self.bgrid[1:]+self.bgrid[:-1])
        mask=np.ones((len(bgridmid),len(b1mid)))
        for i in range(0,len(b1mid)):
           mask[bgridmid<min(b1mid[i],b2mid[i]),i]=0 
           mask[bgridmid>max(b1mid[i],b2mid[i]),i]=0 
           sumi=sum(mask[:,i])
           if sumi>0:
              mask[:,i]=mask[:,i]/sumi
           else:
             #  if no buoyancy levels falls between b1 and b2, find the closest b-level and put transport there:
             idx = (np.abs(bgridmid-0.5*b1mid[i]-0.5*b2mid[i])).argmin()  
             mask[idx,i]=1.
        transb=np.zeros(np.shape(bgridmid));
        for i in range(0,len(bgridmid)):
            transb[i]=np.sum(udydz*mask[i,:])
        return np.append([0],np.cumsum(transb))
    
    def Psibz(self):
        psib=self.Psib()
        # This would do a linear interploation in b: return np.interp(self.b1,self.bgrid,psib)
        # instead we first estimate the depth levels for the bgrid and then do linear interpolation in z
        z1_of_bgrid=np.interp(self.bgrid,self.b1(self.z),self.z)
        #b2eff=0.01*self.b1+0.99*self.b2;
        z2_of_bgrid=np.interp(self.bgrid,self.b2(self.z),self.z)
        return [np.interp(self.z,z1_of_bgrid,psib),np.interp(self.z,z2_of_bgrid,psib)]
   
    
    def update(self,b1=None,b2=None):
        # update buoyancy profiles     
        if b1 is not None:
            self.b1=self.make_func(b1,'b1',self.z)
        if b2 is not None:
            self.b2=self.make_func(b2,'b2',self.z)
        