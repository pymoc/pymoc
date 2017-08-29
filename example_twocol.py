'''
This script shows an example of a "two column" model for the 
overturning circulation in a basin. The first column represents the 
basin, while the second column represents the northern sinking region
'''

from model_PsiNA import Model_PsiNA
from model_vertadvdiff import Model_VertAdvDiff
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs=0.03; bbot=-0.0004
# We will here assume b_N=0 (the default)

A_basin=8e13  #area of the basin
A_north=A_basin/100.  #area of northern sinking region

# The next few lines are a an example for a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
kappa_back=1e-5
kappa_s=3e-5
kappa_4k=3e-4
kappa = lambda z: (kappa_back + kappa_s*np.exp(z/100)+kappa_4k*np.exp(-z/1000 - 4))  

# create grid
z=np.asarray(np.linspace(-3500, 0, 100))
#z=np.asarray(np.linspace(-(3500.**(3./4.)), 0, 40))
#z[:-1]=-(-z[:-1])**(4./3.) # this is for an irregular grid that's fier near the surface


# create initial guess for buoyancy profile in the basin
def b_basin(z): return bs*np.exp(z/300.)+bbot

# create overturning model instance
AMOC = Model_PsiNA(z=z,b_basin=b_basin,b_N=0.)
# and solve for initial overturning streamfunction:
AMOC.solve()

# Create figure and plot initial conditions:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(AMOC.Psi_N, AMOC.z)
#ax2.plot(m.b_N, m.z)
ax2.plot(b_basin(z), z)
plt.ylim((-4e3,0))
ax1.set_xlim((-5,20))
ax2.set_xlim((-0.01,0.04))


# create adv-diff column model instance for basin
basin= Model_VertAdvDiff(z=z,kappa=kappa,b=b_basin,bs=bs,bbot=bbot)
# create adv-diff column model instance for basin
north= Model_VertAdvDiff(z=z,kappa=kappa,b=0.,bs=0.002,bbot=bbot)

# loop to iteratively find equilibrium solution
for ii in range(0, 12*100):    
   # update buoyancy profile
   wb=AMOC.Psi_N*1e6/A_basin
   wN=-AMOC.Psi_N*1e6/A_north
   basin.timestep(w=wb,dt=30*86400)
   north.timestep(w=wN,dt=30*86400,do_conv=True)
   
   # update overturning streamfunction
   #(notice that we are only adjusting b some of the way, which leads to better convergence):
   AMOC.update(b_basin=basin.b,b_N=north.b)
   AMOC.solve()
   
   if ii%120==0:
      # Plot updated results:
      ax1.plot(AMOC.Psi_N, AMOC.z, linewidth=0.5)
      ax2.plot(basin.b, basin.z, linewidth=0.5)
      

# Plot final results:
ax1.plot(AMOC.Psi_N, AMOC.z,linewidth=2)
ax2.plot(basin.b, basin.z, linewidth=2)
ax2.plot(north.b, basin.z, linewidth=2)

