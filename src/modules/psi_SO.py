'''
A 1D model of the overturninng transpot in the S.O.,
given the density profile in the basin [b(z)],
the surface buoyancy in the SO [bs(y)],
and the surface wind stress [tau(y)]
'''

import numpy as np
from scipy import integrate, optimize

class Psi_SO(object):
    def __init__(
            self,            
            z=None,         # vertical grid (array, in)
            y=None,         # horizontal grid (array, in)  
            b=None,         # buoyancy profile at northern end of ACC (function, array or float, in)
            bs=None,        # surface buoyancy (function, array or float, in)
            tau=None,       # surface wind stress (function, array or float, in)  
            f=1.2e-4,       # Coriolis parameter (in)
            rho=1030,       # Density of sea water (in)
            L =1e7,         # Zonal length of the ACC (in)
            KGM=1e3,        # GM coefficient (in)
            c = None,       # phase speed for F2010 BVP smoother of GM streamfunction
            bvp_with_Ek=False, # if true, apply boundary condition that Psi_GM=-Psi_EK in F2010 BVP smoother 
            Hsill=None,     # height (in m above ocean floor) of the "sill", where Psi_Ek is tapered
            HEk=None,       # depth of surface Ekman layer
            Htapertop=None, # A quadratic tapering of the GM streamfunction at the surface
            Htaperbot=None, # A quadratic tapering of the GM streamfunction at the bottom
            smax=0.01,      # maximum slope for clipping of GM streamfunction
            Psi_Ek=None,    # eulerian overturning streamfunction (array, out) 
            Psi_GM=None,    # GM-type eddy overturning streamfunction (array, out) 
            Psi=None,       # residual overturning streamfunction (array, out)
    ):
 
        # initialize grid:
        if isinstance(z,np.ndarray):
            self.z = z
        else:
            raise TypeError('z needs to be numpy array providing grid levels') 
                      
        if isinstance(y,np.ndarray):
            self.y = y
        else:
            raise TypeError('y needs to be numpy array providing horizontal grid (or boundaries) of ACC') 
        
            
        self.b = self.make_func(self.z,b,'b')     
        self.bs = self.make_func(self.y,bs,'bs')     
        self.tau = self.make_func(self.y,tau,'tau')     
        self.f = f
        self.rho = rho
        self.L = L
        self.KGM = KGM
        self.c = c
        self.bvp_with_Ek = bvp_with_Ek
        self.Hsill = Hsill
        self.HEk = HEk
        self.Htapertop = Htapertop
        self.Htaperbot = Htaperbot
        self.smax = smax
    
    def make_func(self,xi,myst,name):
    # turn mysterious object into callable function (if needed)    
        if callable(myst):
            return myst
        elif isinstance(myst,np.ndarray):
            def funfun(x): return np.interp(x,xi,myst)
            return funfun
        elif isinstance(myst,float):
            def funfun(x): return myst +0*x
            return funfun
        else:
            raise TypeError(name,'needs to be either function, numpy array, or float') 
    
    
    def ys(self, b):
        # inverse of bs(y)
        def func(y):
          return self.bs(y) - b

        if b < np.min(self.bs(self.y)):
          # if b is smaller minimum bs, isopycnals don't outcrop and get handled separately
           return self.y[0] - 1e3
        if b > self.bs(self.y[-1]):
          # if b is larger than bs at northern end, return northernmost point:
          return self.y[-1]
        else:
          # if b in range of bs return ys(b):
          # Notice that this inversion is well defined only if bs is monotonically
          # increasing (past minind). Should probably add a check to make sure
          # this is the case...
          minind = np.argmin(self.bs(self.y))
          return optimize.brentq(func, self.y[minind], self.y[-1])
          
    def calc_N2(self):
        dz = self.z[1:] - self.z[:-1];
        N2 = np.zeros(np.size(self.z))
        b = self.b(self.z)

        N2[1:-1] = (b[2:] - b[:-2])/(dz[1:] + dz[:-1])
        N2[0] = (b[1] - b[0])/dz[0]
        N2[-1] = (b[-1] - b[-2])/dz[-1]

        return self.make_func(self.z, N2, 'N2')
    
    def calc_Ekman(self):
        # compute Ekman transport on z grid
        # based on average wind stress between outcrop and northern end of channel 
        tau_ave = 0*self.z
        for ii in range(0,np.size(self.z)):
           y0 = self.ys(self.b(self.z[ii])) # outcrop latitude
           tau_ave[ii] =  np.mean(self.tau(np.linspace(y0,self.y[-1], 100)))
        if self.Hsill is not None:
           silltaper = 1.-np.maximum(self.z[0]+self.Hsill-self.z,0.)**2./self.Hsill**2.
        else:
           silltaper = 1.
        if self.HEk is not None: 
           Ektaper = 1-np.maximum(self.z+self.HEk,0)**2./self.HEk**2.
        else:
           Ektaper = np.ones(np.size(self.z))
           Ektaper[-1] = 0.
        return tau_ave/self.f/self.rho*self.L*silltaper*Ektaper   
    
    
    def calc_GM(self):
        # compute GM ransport on z grid
        # based on average isopycnal slope
        dy_atz = 0*self.z
        eps=0.1 # minimum dy (in meters) (to avoid div. by 0)
        for ii in range(0,np.size(self.z)):
           dy_atz[ii]=max(self.y[-1]-self.ys(self.b(self.z[ii])),eps)
        if self.Htaperbot is not None:
           bottaper=1.-np.maximum(self.z[0]+self.Htaperbot-self.z,0.)**2./self.Htaperbot**2.
        else:
           bottaper=1.
        if self.Htapertop is not None:
           toptaper=1-np.maximum(self.z+self.Htapertop,0)**2./self.Htapertop**2.
        else:
           toptaper=1.
        if self.c is not None:
           temp=self.make_func(self.z,self.KGM*self.z/dy_atz*self.L*toptaper*bottaper,'psiGM')
           N2=self.calc_N2()
           if self.bvp_with_Ek:
              def bc(ya, yb):
                  return np.array([ya[0]+self.Psi_Ek[0]*1e6, yb[0]+self.Psi_Ek[-1]*1e6]) 
           else:    
              def bc(ya, yb): return np.array([ya[0], yb[0]])
           def ode(z, y): 
              return np.vstack((y[1], N2(z)/self.c**2.*(y[0] -temp(z))))
           #Solve the boundary value problem
           res = integrate.solve_bvp(ode, bc, self.z, np.zeros((2, np.size(self.z))))
           # return solution interpolated onto original grid        
           temp= res.sol(self.z)[0, :]
        else:
           temp= self.KGM*np.maximum(self.z/dy_atz,-self.smax)*self.L*toptaper*bottaper
        # limit Psi_GM to -Psi_Ek on isopycnals that don't outcrop:
        temp[dy_atz>self.y[-1]-self.y[0]]=np.maximum(
             temp[dy_atz>self.y[-1]-self.y[0]],
            -self.Psi_Ek[dy_atz>self.y[-1]-self.y[0]]*1e6)
        return temp
                        
    
    def solve(self):
        #Solve for overturning circ.
        self.Psi_Ek=self.calc_Ekman()/1e6;
        self.Psi_GM=self.calc_GM()/1e6;
        self.Psi=self.Psi_Ek + self.Psi_GM;
        # Notice that the bottom boundary is somewhat poorly defined, as it is
        # only for BC. To avoid random fluctuation in bottom Psi, we here simply 
        # set it to value of last gridpoint above
        self.Psi[0]=self.Psi[1]         
        
    
    def update(self,b=None,bs=None):
        # update buoyancy profiles
        if b is not None:
            self.b=self.make_func(self.z,b,'b')
        if bs is not None:
            self.bs=self.make_func(self.y,bs,'bs')     
          
    
        