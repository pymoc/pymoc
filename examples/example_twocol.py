'''
This script shows an example of a "two-column" model for the 
diffusive overturning circulation in a basin, integrated via time-stepping.
The first column represents the basin, while the second column represents the
northern sinking region. In this example we are using an isopycnal mapping of
overturning circulation, i.e. the overturning is converted into isopycnal space
and then mapped onto the different water columns. In this configuration
the model is effectively isopycnal, with the vertical coordinate representing
the mean depth of the respective isopycnal in the column.
'''

from pymoc.modules import Psi_Thermwind, Column
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs = 0.03
bs_north = 0.0
bbot = -0.003

A_basin = 8e13    #area of the basin
A_north = A_basin / 100.    #area of northern sinking region

# time-stepping parameters:
dt = 86400 * 30    # time-step for vert. adv. diff. calc.
MOC_up_iters = int(
    np.floor(2 * 360 * 86400 / dt)
)    # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters = int(
    np.ceil(500 * 360 * 86400 / dt)
)    # plotting frequency (in iterations)
total_iters = int(
    np.ceil(4000 * 360 * 86400 / dt)
)    # total number of timesteps

# The next few lines define a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
kappa_back = 1e-5
kappa_4k = 3e-4


def kappa(z):
  return (kappa_back + kappa_4k * np.exp(-z / 1000 - 4))


# create grid:
z = np.asarray(np.linspace(-4000, 0, 80))


# Initial conditions for buoyancy profile in the basin:
def b_basin(z):
  return bs * np.exp(z / 300.)
def b_north(z):
  return 1e-3*bs * np.exp(z / 300.)



# create overturning model instance
AMOC = Psi_Thermwind(z=z, b1=b_basin, b2=b_north)
# and solve for initial overturning streamfunction:
AMOC.solve()
# evaluate overturning in isopycnal space:
[Psi_iso_b, Psi_iso_n] = AMOC.Psibz()

# Create figure:
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3, 0))
ax1.set_xlim((-5, 20))
ax2.set_xlim((-0.01, 0.04))
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)

# create adv-diff column model instance for basin
basin = Column(
    z=z, kappa=kappa, Area=A_basin, b=b_basin, bs=bs, bbot=bbot
)
# create adv-diff column model instance for basin
north = Column(
    z=z, kappa=kappa, Area=A_north, b=b_north, bs=bs_north, bbot=bbot
)

# Main time stepping loop
for ii in range(0, total_iters):
  # update buoyancy profile
  wAb = Psi_iso_b * 1e6
  wAN = -Psi_iso_n * 1e6
  basin.timestep(wA=wAb, dt=dt)
  north.timestep(wA=wAN, dt=dt, do_conv=True)
  if ii % MOC_up_iters == 0:
    # update overturning streamfunction (can be done less frequently)
    AMOC.update(b1=basin.b, b2=north.b)
    AMOC.solve()
    # evaluate overturning in isopycnal space:
    [Psi_iso_b, Psi_iso_n] = AMOC.Psibz()
  if ii % plot_iters == 0:
    # Plot current state:
    ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5)
    ax2.plot(basin.b, basin.z, linewidth=0.5)
    ax2.plot(north.b, north.z, linewidth=0.5)
    plt.pause(0.01)

# Plot final results:
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(AMOC.Psi, AMOC.z, linewidth=2, color='r')
ax2.plot(basin.b, basin.z, linewidth=2, color='b')
ax2.plot(north.b, basin.z, linewidth=2, color='c')
ax1.plot(0. * AMOC.z, AMOC.z, linewidth=0.5, color='k')
ax1.set_xlim((-5, 15))
ax2.set_xlim((-0.01, 0.03))
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)
plt.ylim((-4e3, 0))
