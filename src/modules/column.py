'''
Instances of this class represent 1D representations of buoyancy in a water
column governed by vertical advection and diffusion. The velocity
profile is required as an input. The script can either compute the equilibrium
buoyancy profile for a given vertical velocity profile and boundary conditions
or compute the tendency and perform a time-step of given length.
BCs have to be fixed buoyancy at the top and either fixed b or db/dz at the bottom
The time-stepping version can also handle horizontal advection
into the column. This is, however, not (yet) implemented for the equilibrium solver
'''

import sys
import numpy as np
from scipy import integrate
sys.path.append('/pymoc/src/utils')
from make_func import make_func 

class Column(object):
    # This module creates an advective-diffusive column
    # Notice that the column here represents a horizontal integral, rather than
    # an average, thus allowing for the area of to be a function of depth
    def __init__(
            self,
            z=None,         # grid (input)
            kappa=None,     # diffusivity profile (input)
            bs=0.025,       # surface buoyancy bound. cond (input)
            bbot=0.0,       # bottom buoyancy boundary condition (input)  
            bzbot=None,     # bottom strat. as alternative boundary condition (input) 
            b=0.0,          # Buoyancy profile (input, output)
            Area=None,      # Horizontal area (can be function of depth)
            N2min=1e-7      # Minimum strat. for conv adjustment
    ):
 
        # initialize grid:
        if isinstance(z,np.ndarray) and len(z) > 0:
            self.z = z
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
                      
        self.kappa=self.make_func(kappa,'kappa')
        self.Area=self.make_func(Area,'Area')
        
        self.bs=bs
        self.bbot=bbot 
        self.bzbot=bzbot  
       
        self.N2min=N2min  
       
        self.b=self.make_array(b,'b')     
        
        if self.check_numpy_version():
           self.bz=np.gradient(self.b, z)
        else:
           self.bz=0.*z   # notice that this is just for initialization of ode solver
           
    def check_numpy_version(self):
        # check numpy version (version >= 1.13 needed to automatically compute db/dz)    
        v = [int(i) for i in np.version.version.split('.')]
        if v[0] <= 1 and v[1] < 13:
            return False
        return True
       
    
    def make_func(self,myst,name): return make_func(myst, self.z, name)
    
    def make_array(self, myst, name):
      # turn mysterious object into array(if needed)    
      if isinstance(myst, np.ndarray):
        return myst
      elif callable(myst):
        return myst(self.z)
      elif isinstance(myst, float):
        return myst + 0 * self.z
      else:
        raise TypeError(name, 'needs to be either function, numpy array, or float') 
    
    def Akappa(self, z):
      return self.Area(z)*self.kappa(z)
    
    def dAkappa_dz(self, z):
      if not self.check_numpy_version():
        raise ImportError('You need NumPy version 1.13.0 or later. Please upgrade your NumPy libary.')
      return np.gradient(self.Akappa(z), z)
    
    
    def bc(self, ya, yb):
      #return the boundary conditions
      if self.bzbot is None:
        return np.array([ya[0] - self.bbot, yb[0] - self.bs])
      else:
        return np.array([ya[1] - self.bzbot, yb[0] - self.bs])
        
    
    def ode(self, z, y):
        #return the equation to be solved 
        return np.vstack((y[1], (self.wA(z) - self.dAkappa_dz(z)) / self.Akappa(z) * y[1]))
                         
    def solve_equi(self,wA):
        #Solve for equilibrium solution given vert. vel profile and BCs
        self.wA = self.make_func(wA,'w')
        sol_init = np.zeros((2, np.size(self.z)))
        sol_init[0,:] = self.b
        sol_init[1,:] = self.bz
        res = integrate.solve_bvp(self.ode, self.bc, self.z, sol_init)
        # interpolate solution for b and db/dz onto original grid
        self.b = res.sol(self.z)[0, :]  
        self.bz = res.sol(self.z)[1, :]  
        
    
    def vertadvdiff(self,wA,dt):
      #upwind vert. adv. and diffusion
      wA = self.make_array(wA,'wA')
      dz = self.z[1:] - self.z[:-1]

      # apply boundary conditions:        
      self.b[-1] = self.bs
      self.b[0] = self.bbot if self.bzbot is None else self.b[1] - self.bzbot*dz[0]

      bz = (self.b[1:] - self.b[:-1]) / dz
      bz_up = bz[1:]
      bz_down = bz[:-1]
      bzz = (bz_up - bz_down) / (0.5*(dz[1:] + dz[:-1]))

      #upwind advection: 
      weff = wA - self.dAkappa_dz(self.z)
      bz = bz_down
      bz[weff[1:-1] < 0] = bz_up[weff[1:-1] < 0]      

      db_dt = (
        -weff[1:-1] * bz / self.Area(self.z[1:-1])
        + self.kappa(self.z[1:-1]) * bzz
      )
      self.b[1:-1] = self.b[1:-1]+dt*db_dt
        
    
    def convect(self):
        # do convective adjustment to minimum strat N2min
        # notice that this parameterization currently only handles convection
        # from the top down which is the only case we really encounter here...
        # it also assumes a fixed surface buoyancy BC, and hence any
        # surface heat flux that's required to adjust the buoyancy of the
        # convecting column:
        ind = self.b > self.b[-1] + self.N2min*self.z
        self.b[ind] = self.b[-1] + self.N2min*self.z[ind]
        # Below is an energy conserving version that could be used 
        # for model formulations without fixed surface b. But for fixed surface
        # b, the simpler version above is preferable as it deals better with long time steps
        # (for infinitesimal time-step and fixed surface b, the two are equivalent)
        # this version does currently also not include adjustment to finite strat.
        # dz=self.z[1:]-self.z[:-1];
        # dz=np.append(dz,dz[-1]);dz=np.insert(dz,0,dz[0])
        # dz=0.5*(dz[1:]+dz[0:-1]);
        # self.b[ind]=(np.mean(self.b[ind]*dz[ind]*self.Area(self.z[ind]))
        #            /np.mean(dz[ind]*self.Area(self.z[ind])) )   
            
    def horadv(self,vdx_in,b_in,dt):
        # upwind horizontal advection:
        # notice that vdx_in is the total transport per unit height
        # into the column (units m^2/s, sign positive for velocity into the column)
        vdx_in = self.make_array(vdx_in,'vdx_in')  
        b_in = self.make_array(b_in,'b_in')  

        adv_idx = vdx_in > 0.0
        db = b_in - self.b

        self.b[adv_idx] = self.b[adv_idx] + dt*vdx_in[adv_idx]*db[adv_idx] / self.Area(self.z[adv_idx])
        
    
    def timestep(self,wA=0.,dt=1.,do_conv=False,vdx_in=None,b_in=None):
        #Integrate buoyancy profile evolution for one time-step
        # do vertical advection and diffusion:
        self.vertadvdiff(wA=wA,dt=dt)
        if vdx_in is not None:  
           # do horizontal advection: (optional)
           if b_in is not None: 
              self.horadv(vdx_in=vdx_in,b_in=b_in,dt=dt)
           else:
              raise TypeError('b_in is needed if vdx_in is provided')
        if do_conv:
           # do convection: (optional)
           self.convect()
        
