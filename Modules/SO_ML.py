'''
Instances of this class represent 1D representations of buoyancy 
in the S.O. mixed layer. An advective-diffusive equation is solved.
Surface fluxes can be prescribed via a fixed flux and/or restoring
Northen and southern boundary conditions are determined internally based on the
overturning streamfunction and buoyancy profile in the basin
'''

import numpy as np

class SO_ML(object):
    # This module creates an advective-diffusive column
    # Notice that the column here represents a horizontal integral, rather than
    # an average, thus allowing for the area of to be a function of depth
    def __init__(
            self,
            y=None,            # grid (input)
            Ks=0.,             # hor. diffusivity (input)
            h=50.,             # ML depth (input)
            L=4e6,             # zonal width (input)
            surflux=0.,        # prescribed surface buoyancy flux (in m^2/s^3; input)
                               # Notice that the first and last gridpoints are BCs and should have no flux 
            rest_mask=0.,      # mask for surface restoring (1 where restoring is applied 0 elsewhere)
            b_rest=0.,         # surface buoyancy towards which we are restoring 
            v_pist=1.5/86400., # Piston velocity for restoring in SL (input)
            bs=0.0,            # Surface buoyancy (input, output)
            Psi_s=None         # Overturning in the ML (output)
    ):
 
        # initialize grid:
        if isinstance(y,np.ndarray):
            self.y = y
        else:
            raise TypeError('y needs to be numpy array providing (regular) grid') 
                              
        self.Ks=Ks; self.h=h; self.L=L
        self.surflux=self.make_array(surflux,'surflux')
        self.rest_mask=self.make_array(rest_mask,'rest_mask')
        self.b_rest=self.make_array(b_rest,'b_rest')
        self.v_pist=v_pist; self.Psi_s=Psi_s
        self.bs=self.make_array(bs,'bs') 
 
    
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
    

                         
    def solve_equi(self):
        #Solve for equilibrium solution given inputs
        raise TypeError('This functionality is not yet implemented') 
        
        
    
    def advdiff(self, b_basin, Psi_b, dt):
      # update surface buoyancy profile via advect. and diff
      
      # First we need to determine Psi at the surface in the ML:
      # The first line here is a hack to reduce problems with interpolation  
      # due to finite SO resolution when PsiSO goes to zero at the bottom
      #due to non-outcropping isopycnals
      # If Psi_b is zero at non-outcropping isopycnals, it can end up 
      # very near zero also at the last (non-boundary) gridpoint
      # at the surface as a result of the interpolation procedure 
      # - this is unphysical since psi by definition only vanishes on isopycnals 
      # that don't outcrop (the problem is that this vanishing here is a step
      # function which messes with the interpolation.) 
      # For the purpose of the interpolation to the surface we therefore set psi
      # on non-outcropping isopycnals to the last non-zero value above 
      ind=np.nonzero(Psi_b)[0][0]; Psi_b[:ind]=Psi_b[ind]
      self.Psi_s=np.interp(self.bs,b_basin,Psi_b)
      self.Psi_s[0]=0.# This value doesn't actually enter/matter, but zero overturning
                      # at southern boundary makes more sense for diag purposes
      
      #set boundary conditions:
      self.bs[-1]=b_basin[-1]
      if self.Psi_s[0]>0:
          # set buoyancy at southern boundary buoyancy to mean of upwelling water
          #bup=int(d_z(psi)b)dz/int(d_z(psi))dz
          dpsidz=Psi_b[2:]-Psi_b[:-2]
          bint=b_basin[1:-1];
          ind=np.logical_and(Psi_b[1:-1]>0., b_basin[1:-1]<self.bs[1]) 
          bup=np.sum(dpsidz[ind]*bint[ind])/np.sum(dpsidz[ind])
          self.bs[0]=bup
      else:
          # no-flux BC
          self.bs[0]=self.bs[1]
   
      # Compute tendency due to surface b-flux
      dbdt_flux=self.surflux/self.h + self.rest_mask*self.v_pist/self.h*(
                               self.b_rest-self.bs);
      # Compute advective tendency via upwind advection
      dy= self.y[1]-self.y[0] # Notice that the current implementation assumes an evenly spaced grid!
      dbdt_ad=0.*self.y;
      indneg=self.Psi_s[1:-1]<0.;
      indpos=self.Psi_s[1:-1]>0.;
      dbdt_ad[1:-1][indneg]=-self.Psi_s[1:-1][indneg]*1e6*(self.bs[2:][indneg]-self.bs[1:-1][indneg])/self.h/self.L/dy
      dbdt_ad[1:-1][indpos]=-self.Psi_s[1:-1][indpos]*1e6*(self.bs[1:-1][indpos]-self.bs[:-2][indpos])/self.h/self.L/dy
       
      # add tendencies from surface flux and advection
      self.bs=self.bs+dt*(dbdt_flux+dbdt_ad)       
      
      # Do implicit diffusion:
      s=self.Ks*dt/dy**2;
      U=(np.diag(-s/2.*np.ones(len(self.y)-1), -1)
       + np.diag((1+s)*np.ones(len(self.y)), 0)
       + np.diag(-s/2.*np.ones(len(self.y)-1), 1))
      U[0,0]=1;U[0,1]=0; U[-1,-2]=0;U[-1,-1]=1;
      Uinv=np.linalg.inv(U)
      V=(np.diag(s/2.*np.ones(len(self.y)-1), -1)
       + np.diag((1-s)*np.ones(len(self.y)), 0)
       + np.diag(s/2.*np.ones(len(self.y)-1), 1))
      V[0,0]=1;V[0,1]=0; V[-1,-2]=0;V[-1,-1]=1;
      self.bs=np.dot(np.dot(Uinv,V),self.bs)      
      
      
    
    
    def timestep(self, b_basin=None, Psi_b=None ,dt=1.):
        #Integrate buoyancy profile evolution for one time-step
        if not isinstance(b_basin,np.ndarray):
            raise TypeError('b_basin needs to be numpy array providing buoyancy levels in basin') 
        
        if not isinstance(Psi_b,np.ndarray):
            raise TypeError('Psi_b needs to be numpy array providing overturning at buoyancy levels given by b_basin') 
   
        self.advdiff(b_basin=b_basin,Psi_b=Psi_b,dt=dt)
        
