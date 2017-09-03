'''
This script shows an example of a "two column" model for the 
overturning circulation in a basin connected to a channel in the south.
The first column represents the basin, while the second column represents
the northern sinking region. The overtunning circulation is computed at
the northern end of the basin (at the interface to the northern sinking region)
and at the southern end of the basin (at the interface to the channel)
The abyssal boundary condition is formulated based on a scaling for the
abyssal BL, following Mashayek et al and others.

'''
import sys
sys.path.append('../Modules')
from model_PsiNA import Model_PsiNA
from model_SO import Model_SO
from model_vertadvdiff import Model_VertAdvDiff
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs=0.03; bs_north=0.0; bbot= -0.01 

# S.O. surface boundary conditions and grid:
y=np.asarray(np.linspace(0,2.e6, 40))
tau=0.1 #*np.sin(np.pi*y/2.e6)
bs_SO=(bs-bbot)*y/y[-1]+bbot

A_basin=8e13  #area of the basin
A_north=A_basin/50.  #area of northern sinking region

# time-stepping parameters:
dt=86400*30                                 # time-step for vert. adv. diff. calc.
MOC_up_iters=int(np.floor(2.*360*86400/dt))  # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters= int(np.ceil(500*360*86400/dt))  # plotting frequency (in iterations)
total_iters=int(np.ceil(5000*360*86400/dt))  # total number of timesteps

# The next few lines define a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
kappa_back=2e-5
kappa_4k=2e-4
def kappa(z):
    return (kappa_back + kappa_4k*np.exp(-z/1000 - 4))  

# create vertical grid:
z=np.asarray(np.linspace(-3200, 0, 80))

# create initial guess for buoyancy profile in the basin
def b_basin(z): return bs*np.exp(z/300.)

# create N.A. overturning model instance
AMOC = Model_PsiNA(z=z,b_basin=b_basin,b_N=0.)
# and solve for initial overturning streamfunction:
AMOC.solve()

# create S.O. overturning model instance
SO=Model_SO(z=z,y=y,b=b_basin(z),bs=bs_SO,tau=tau,L=6e6,KGM=1200.)
SO.solve()

# create adv-diff column model instance for basin
basin= Model_VertAdvDiff(z=z,kappa=kappa,b=b_basin,bs=bs,bbot=bbot)
# create adv-diff column model instance for basin
north= Model_VertAdvDiff(z=z,kappa=kappa,b=0.,bs=bs_north,bbot=bbot)


# Create figure:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-10,15))
ax2.set_xlim((-0.02,0.03))

# loop to iteratively find equilibrium solution
for ii in range(0, total_iters):    
   # update buoyancy profile
   wb=(AMOC.Psi-SO.Psi)*1e6/A_basin
   wN=-AMOC.Psi*1e6/A_north
   # use abyssal BL scaling for bottom boundary condition
   # (activated if d_abyss non-zero):
   basin.d_abyss=-kappa(SO.z[0])*A_basin/SO.Psi[0]/1e6
   north.d_abyss=basin.d_abyss
   basin.timestep(w=wb,dt=dt)
   north.timestep(w=wN,dt=dt,do_conv=True)
   
   if ii%MOC_up_iters==0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(b_basin=basin.b,b_N=north.b)
      AMOC.solve()
      SO.update(b=basin.b)
      SO.solve()
   
   if ii%plot_iters==0:
      # Plot current state:
      ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5, color='r')
      ax1.plot(SO.Psi, SO.z, linewidth=0.5, color='m')
      ax2.plot(basin.b, basin.z, linewidth=0.5,color='b')
      ax2.plot(north.b, north.z, linewidth=0.5,color='c')
      

# Plot final results:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(AMOC.Psi, AMOC.z,linewidth=2,color='r')
ax1.plot(SO.Psi, SO.z,linewidth=2,color='m')
ax1.plot(SO.Psi_Ek, SO.z,linewidth=1,color='m',linestyle='--')
ax1.plot(SO.Psi_GM, SO.z,linewidth=1,color='m',linestyle=':')
ax2.plot(basin.b, basin.z, linewidth=2,color='b')
ax2.plot(north.b, basin.z, linewidth=2,color='c')
ax1.plot(0.*AMOC.z, AMOC.z,linewidth=0.5,color='k')
ax1.set_xlim((-10,15))
ax2.set_xlim((-0.02,0.03))
plt.ylim((-4e3,0))









