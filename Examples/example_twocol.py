'''
This script shows an example of a (time-stepping) "two column" model for the 
overturning circulation in a basin. The first column represents the 
basin, while the second column represents the northern sinking region
'''
import sys
sys.path.append('../Modules')
from model_PsiNA import Model_PsiNA
from model_vertadvdiff import Model_VertAdvDiff
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs=0.03; bs_north=0.0; bzbot= 2e-7 

A_basin=8e13  #area of the basin
A_north=A_basin/100.  #area of northern sinking region

# time-stepping parameters:
dt=86400*30                                 # time-step for vert. adv. diff. calc.
MOC_up_iters=int(np.floor(2*360*86400/dt))  # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters= int(np.ceil(500*360*86400/dt))  # plotting frequency (in iterations)
total_iters=int(np.ceil(5000*360*86400/dt))  # total number of timesteps

# The next few lines define a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
kappa_back=1e-5
kappa_s=3e-5
kappa_4k=3e-4
def kappa(z):
    return (kappa_back + kappa_s*np.exp(z/100)+kappa_4k*np.exp(-z/1000 - 4))  

# create regular grid:
#z=np.asarray(np.linspace(-3500, 0, 80))
# ..or create irregular grid:
z=np.asarray(np.linspace(-(3500.**(3./4.)), 0, 40))
z[:-1]=-(-z[:-1])**(4./3.) 

# Initial conditions for buoyancy profile in the basin:
def b_basin(z): return bs*np.exp(z/300.)

# create overturning model instance
AMOC = Model_PsiNA(z=z,b_basin=b_basin,b_N=0.)
# and solve for initial overturning streamfunction:
AMOC.solve()

# Create figure:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-5,20))
ax2.set_xlim((-0.01,0.04))


# create adv-diff column model instance for basin
basin= Model_VertAdvDiff(z=z,kappa=kappa,Area=A_basin,b=b_basin,bs=bs,bzbot=bzbot)
# create adv-diff column model instance for basin
north= Model_VertAdvDiff(z=z,kappa=kappa,Area=A_north,b=0.,bs=bs_north,bzbot=bzbot)

# Main time stepping loop
for ii in range(0, total_iters):    
   # update buoyancy profile
   wAb=AMOC.Psi*1e6
   wAN=-AMOC.Psi*1e6
   basin.timestep(wA=wAb,dt=dt)
   north.timestep(wA=wAN,dt=dt,do_conv=True)
   
   if ii%MOC_up_iters==0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(b_basin=basin.b,b_N=north.b)
      AMOC.solve()
   
   if ii%plot_iters==0:
      # Plot current state:
      ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5)
      ax2.plot(basin.b, basin.z, linewidth=0.5)
      ax2.plot(north.b, north.z, linewidth=0.5)
      

# Plot final results:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(AMOC.Psi, AMOC.z,linewidth=2,color='r')
ax2.plot(basin.b, basin.z, linewidth=2,color='b')
ax2.plot(north.b, basin.z, linewidth=2,color='c')
ax1.plot(0.*AMOC.z, AMOC.z,linewidth=0.5,color='k')
ax1.set_xlim((-5,15))
ax2.set_xlim((-0.01,0.03))
plt.ylim((-4e3,0))


