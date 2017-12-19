'''
This script shows an example of a "two column" model for the 
overturning circulation in a basin connected to a channel in the south.
The first column represents the basin, while the second column represents
the northern sinking region. The overtunning circulation is computed at
the northern end of the basin (at the interface to the northern sinking region)
and at the southern end of the basin (at the interface to the channel).
The parameters chosen here follow more or less the "control" experiment of Nikurashin
and Vallis (2012, JPO).
'''
import sys
sys.path.append('../Modules')
from model_thermwind import Model_Thermwind
from model_SO import Model_SO
from model_column import Model_Column
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs=0.03; bs_north=0.004; bmin=0.0 

# S.O. surface boundary conditions and grid:
y=np.asarray(np.linspace(0,2.e6, 40))
#Although the model can emulate the effect of a meridionally varying wind-stress,
# that feature is in "beta" and we are here for simplicity using a constant
# wind-stress (approximately the average over the channel in NV12):
tau=0.13
# from what I can infer, NV12 use an approximately quadratic surface temperature
# profile in the channel for their SAMBUCA simulations: 
bs_SO=(bs-bmin)*(y/y[-1])**2+bmin

A_basin=6e13  #area of the basin
A_north=A_basin/50.  #area of northern sinking region

# time-stepping parameters:
dt=86400*30                                  # time-step for vert. adv. diff. calc.
MOC_up_iters=int(np.floor(2.*360*86400/dt))  # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters= int(np.ceil(300*360*86400/dt))   # plotting frequency (in iterations)
total_iters=int(np.ceil(3000*360*86400/dt))  # total number of timesteps

kappa=2e-5

# create vertical grid:
z=np.asarray(np.linspace(-4000, 0, 80))

# Initial conditions for buoyancy profile in the basin
def b_basin(z): return bs*np.exp(z/300.)

# create N.A. overturning model instance
AMOC = Model_Thermwind(z=z,b1=b_basin,b2=0.,f=1e-4)
# and solve for initial overturning streamfunction:
AMOC.solve()
# evaluate overturning to isopycnal space:
[Psi_iso_b,Psi_iso_n]=AMOC.Psibz()

# create S.O. overturning model instance
SO=Model_SO(z=z,y=y,b=b_basin(z),bs=bs_SO,tau=tau,f=1e-4,L=5e6,KGM=1000.,c=0.1, bvp_with_Ek=True)
SO.solve()


# create adv-diff column model instance for basin
basin= Model_Column(z=z,kappa=kappa,Area=A_basin,b=b_basin,bs=bs,bbot=bmin)
# create adv-diff column model instance for basin
north= Model_Column(z=z,kappa=kappa,Area=A_north,b=0.,bs=bs_north,bbot=bmin)


# Create figure:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-10,15))
ax2.set_xlim((-0.02,0.03))

# Main time-stepping loop:
for ii in range(0, total_iters):    
   # update buoyancy profile
   # using z-coordinate overturning:
   #wAb=(AMOC.Psi-SO.Psi)*1e6
   #wAN=-AMOC.Psi*1e6
   # using isopycnal overturning:
   wAb=(Psi_iso_b-SO.Psi)*1e6
   wAN=-Psi_iso_n*1e6
   basin.timestep(wA=wAb,dt=dt)
   north.timestep(wA=wAN,dt=dt,do_conv=True)   
   if ii%MOC_up_iters==0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(b1=basin.b,b2=north.b)
      AMOC.solve()
      [Psi_iso_b,Psi_iso_n]=AMOC.Psibz()
      SO.update(b=basin.b)
      SO.solve()
   if ii%plot_iters==0:
      # Plot current state:
      ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5, color='r')
      ax1.plot(SO.Psi, SO.z, linewidth=0.5, color='m')
      ax2.plot(basin.b, basin.z, linewidth=0.5,color='b')
      ax2.plot(north.b, north.z, linewidth=0.5,color='c')
      plt.pause(0.01)
      

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



