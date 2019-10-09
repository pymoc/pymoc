'''
This script shows an example of how to use the psi_thermwind module
together with the column module to build an iterative solver
for the diffusive overturning circulation in a basin
'''

import sys
from pymoc.modules import Psi_Thermwind
from pymoc.modules import Column
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs = 0.03
bbot = -0.0004
# We will here assume b_N=0 (the default)

A_basin = 8e13    #area of the basin

# The next few lines are an example for a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
kappa_back = 1e-5
kappa_s = 3e-5
kappa_4k = 3e-4
kappa = lambda z: (
    kappa_back + kappa_s * np.exp(z / 100) + kappa_4k * np.exp(-z / 1000 - 4)
)

# create grid
z = np.asarray(np.linspace(-3500, 0, 100))


# create initial guess for buoyancy profile in the basin
def b_basin(z):
  return bs * np.exp(z / 300.) + bbot


# create overturning model instance
AMOC = Psi_Thermwind(z=z, b1=b_basin)
# and solve for initial overturning streamfunction:
AMOC.solve()

# Create figure and plot initial conditions:
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(AMOC.Psi, AMOC.z)
#ax2.plot(m.b_N, m.z)
ax2.plot(b_basin(z), z)
plt.ylim((-4e3, 0))
ax1.set_xlim((-5, 20))
ax2.set_xlim((-0.01, 0.04))
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)

# create adv-diff column model instance for basin
basin = Column(z=z, kappa=kappa, Area=A_basin, b=b_basin, bs=bs, bbot=bbot)

# loop to iteratively find equilibrium solution
for ii in range(0, 30):
  # update buoyancy profile
  wA = AMOC.Psi * 1e6
  basin.solve_equi(wA)

  # update overturning streamfunction
  #(notice that we are only adjusting b some of the way, which leads to better convergence):
  AMOC.update(b1=0.8 * AMOC.b1(z) + 0.2 * basin.b)
  AMOC.solve()

  # Plot updated results:
  ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5)
  ax2.plot(basin.b, basin.z, linewidth=0.5)
  plt.pause(0.01)

# Plot final results:
ax1.plot(AMOC.Psi, AMOC.z, linewidth=2)
ax2.plot(basin.b, basin.z, linewidth=2)
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)
