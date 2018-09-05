'''
A class to interpolate the buoyancy field in the North along lines
with constant slope - only to make fancy plots
'''

import numpy as np
from scipy.optimize import brenth

class Interpolate_twocol(object):
    def __init__(
            self,            
            y=None,          # y-grid
            z=None,          # z-grid
            bs=None,         # buoyancy profile in the south
            bn=None,         # buoyancy profile in the north
    ):
 
        # initialize grid:
        if isinstance(y,np.ndarray):
            self.y = y
        else:
            raise TypeError('y needs to be numpy array providing grid levels') 
        if isinstance(z,np.ndarray):
            self.z = z
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
        
        self.bs=self.make_func(bs,'bs',self.z)
        self.bn=self.make_func(bn,'bn',self.z)
             
        
    def make_func(self,myst,name,xin):
    # turn array or float into callable function (if needed)    
        if callable(myst):
            return myst
        elif isinstance(myst,np.ndarray):
            def funfun(x): return np.interp(x,xin,myst)
            return funfun
        elif isinstance(myst,float):
            def funfun(x): return myst +0*x
            return funfun
        else:
            raise TypeError(name,'needs to be either function, numpy array, or float') 
           
          
    def __call__(self, y, z):
        l=self.y[-1]
        bsurf=self.make_func(self.y/l*self.bn(0) + (1-self.y/l)*self.bs(0),'bsurf',self.y)
        if z==0 and y==0:
            # slope ill defined at (0,0); evaluate infinitesimally below the surface:
            z=-0.01  
        if z==self.z[0] and y==l:
            # slope also potentially ill defined at (l,-H); evaluate infinitesimally above the bottom:
            z=self.z[0]+0.01  
        def fint(x):
               # function to help determine slope at bottom of vent. region
               return self.bn(0)-self.bs(-x*l)  
        def fup(x):
               # function to help determine slope in vent. region
               return self.bs(z-x*y)-bsurf(y-z/x)    
        def fdeep(x):
               # function to help determine slope below vent. region
               return self.bs(z-x*y)-self.bn(z+x*(l-y))            
        # first determine slope at bottom of vent. region here
        sbot=brenth(fint, 0.,1.) 
        # than set slope for stuff above and below...           
        if z>-sbot*(l-y):
          s=brenth(fup, 1e-10,1.0)
        else:
          s=brenth(fdeep, -1.0,1.0)
        return self.bs(z-s*y)
    
    def gridit(self):
        ny=len(self.y);nz=len(self.z)
        barray=np.zeros((ny,nz))
        for iy in range(0,ny):
            for iz in range(0,nz):
               barray[iy,iz]=self(self.y[iy],self.z[iz])
        return barray
      
    