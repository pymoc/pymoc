'''
This script shows an example of a three column model for the 
overturning circulation in two basins connected to a channel in the south.
The model is identical two the one in twobasin_NadeauJansen.py,
except that we are here making use of the "model" class to couple the
various circulation modules.
'''

from pymoc import model
from pymoc.modules import Psi_Thermwind, Psi_SO, Column
from pymoc.plotting import Interpolate_channel, Interpolate_twocol
import numpy as np
from matplotlib import pyplot as plt

diag_file = None

# boundary conditions:
bs = 0.02
# the minimum surface buoyancy in the GCM in the North Atlantic is 3.5890e-04 (or 0.1795C)
bs_north = 0.00036
bAABW = -0.0011
bbot = min(bAABW, bs_north)

# S.O. surface boundary conditions and grid:
y = np.asarray(np.linspace(0, 3.e6, 51))
tau = 0.16
offset = 0.0345 * (1 - np.cos(np.pi * (5.55e5-1.e5) / 8e6))
bs_SO = (
    0.0345 * (1 - np.cos(np.pi * (y-1.e5) / 8e6)) * (y > 5.55e5) +
    (bAABW-offset) / 5.55e5 * np.maximum(0, 5.55e5 - y) + offset *
    (y < 5.55e5)
)

# time-stepping parameters:
dt = 86400. * 30.    # time-step for vert. adv. diff. calc.
MOC_up_iters = int(
    np.floor(2. * 360 * 86400 / dt)
)    # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters = int(
    np.ceil(500 * 360 * 86400 / dt)
)    # plotting frequency (in iterations)
total_iters = int(
    np.ceil(4000 * 360 * 86400 / dt)
)    # total number of timesteps


# Effective diffusivity profile
def kappaeff(z):    # effective diffusivity profile with tapering in BBL
  return 1.0 * (
      1e-4 * (
          1.1 - np.tanh(
              np.maximum(z + 2000., 0) / 1000. +
              np.minimum(z + 2000., 0) / 1300.
          )
      ) * (1. - np.maximum(-4000. - z + 600., 0.) / 600.)**2
  )


A_Atl = 7e13    #6e13
A_north = 5.5e12
A_Pac = 1.7e14    #1.44e14

Lx = 1.3e+07    #(length of the channel)
Latl = 6. / 21. * Lx
Lpac = 15. / 21. * Lx
K = 1800.
N2min = 2e-7

# create vertical grid:
z = np.asarray(np.linspace(-4000, 0, 80))


# create initial guess for buoyancy profile in the Atl
def b_Atl(z):
  return bs * np.exp(z / 300.) + z / z[0] * bbot


#def b_Atl(z): return 0.3*bs+z/z[0]*(bbot-0.3*bs)
def b_Pac(z):
  return bs * np.exp(z / 300.) + z / z[0] * bbot

# Create new Model instance
m = model.Model()

# Created and add module for N.A. overturning to model
m.new_module(Psi_Thermwind, {
    'z': z,
}, 'AMOC')

# create and add module for interbasin zonal overturning
m.new_module(Psi_Thermwind, {'z': z, 'f': 1e-4}, 'ZOC')

# create and add S.O. overturning module for Atlantic sector
m.new_module(
    Psi_SO, {
        'z': z,
        'y': y,
        'b': b_Atl(z),
        'bs': bs_SO,
        'tau': tau,
        'L': Latl,
        'KGM': K
    }, 'SO_Atl'
)

# create and add S.O. overturning mododule for Pacific sector
m.new_module(
    Psi_SO, {
        'z': z,
        'y': y,
        'b': b_Pac(z),
        'bs': bs_SO,
        'tau': tau,
        'L': Lpac,
        'KGM': K
    }, 'SO_Pac'
)

# create and add adv-diff column module for Atl. basin
m.new_module(
    Column, {
        'z': z,
        'kappa': kappaeff,
        'b': b_Atl,
        'bs': bs,
        'bbot': bbot,
        'Area': A_Atl,
        'N2min': N2min
    },
    'Atl',
    neighbors=[{
        'module': m.get_module('amoc'),
        'direction': 'right',
        'module': m.get_module('so_atl'),
        'direction': 'left',
        'module': m.get_module('zoc'),
        'direction': 'right'
    }]
)

# create and add adv-diff column module for northern sinking region
m.new_module(
    Column, {
        'z': z,
        'kappa': kappaeff,
        'b': b_Atl,
        'bs': bs_north,
        'bbot': bbot,
        'Area': A_north,
        'N2min': N2min,
        'do_conv': True
    },
    'North',
    neighbors=[{
        'module': m.get_module('amoc'),
        'direction': 'left'
    }]
)

# create and add adv-diff column module for Pac. basin
m.new_module(
    Column, {
        'z': z,
        'kappa': kappaeff,
        'b': b_Pac,
        'bs': bs,
        'bbot': bbot,
        'Area': A_Pac,
        'N2min': N2min
    },
    'Pac',
    neighbors=[{
        'module': m.get_module('zoc'),
        'direction': 'left',
        'module': m.get_module('so_pac'),
        'direction': 'left'
    }]
)

# Create figure:
fig = plt.figure(figsize=(6, 9))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3, 0))
ax1.set_xlim((-20, 30))
ax2.set_xlim((-0.02, 0.030))

for snapshot in m.run(
    basin_dt=dt,
    coupler_dt=MOC_up_iters,
    steps=total_iters,
    snapshot_start=plot_iters,
    snapshot_interval=plot_iters
):
  # Plot snapshots of current state at prescribed interval:
  ax1.plot(m.amoc.Psi, m.amoc.z, '--r', linewidth=0.5)
  ax1.plot(m.zoc.Psi, m.zoc.z, ':c', linewidth=0.5)
  ax1.plot(m.so_atl.Psi - m.zoc.Psi, m.so_atl.z, '--b', linewidth=0.5)
  ax1.plot(m.so_pac.Psi + m.zoc.Psi, m.so_pac.z, '--g', linewidth=0.5)
  ax1.plot(m.so_atl.Psi + m.so_pac.Psi, z, '--c', linewidth=0.5)
  ax2.plot(m.atl.b, m.atl.z, '-b', linewidth=0.5)
  ax2.plot(m.pac.b, m.pac.z, '-g', linewidth=0.5)
  ax2.plot(m.north.b, m.north.z, '-r', linewidth=0.5)
  plt.pause(0.01)
 
# Plot final state:
ax1.plot(m.amoc.Psi, m.amoc.z, '--r', linewidth=1.5)
ax1.plot(m.zoc.Psi, m.zoc.z, ':c', linewidth=1.5)
ax1.plot(m.so_atl.Psi - m.zoc.Psi, m.so_atl.z, '--b', linewidth=1.5)
ax1.plot(m.so_pac.Psi + m.zoc.Psi, m.so_pac.z, '--g', linewidth=1.5)
ax1.plot(m.so_atl.Psi + m.so_pac.Psi, z, '--c', linewidth=1.5)
ax2.plot(m.atl.b, m.atl.z, '-b', linewidth=1.5)
ax2.plot(m.pac.b, m.pac.z, '-g', linewidth=1.5)
ax2.plot(m.north.b, m.north.z, '-r', linewidth=1.5)

plt.show()
