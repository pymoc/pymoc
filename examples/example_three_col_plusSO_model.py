'''
This script shows an example of a "two column" model for the 
overturning circulation in a basin connected to a channel in the south.
The example is the same as "example_twocol_plusSO" but we are here
using the "model" class to patch together the different components
of the setup.
'''
from pymoc import model
from pymoc.modules import Psi_Thermwind, Psi_SO, Column, SO_ML
from pymoc.plotting import Interpolate_channel
import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

# boundary conditions:
bs = 0.02
# the minimum surface buoyancy in the GCM in the North Atlantic is 3.5890e-04 (or 0.1795C)
bs_north = 0.0004
bs_east = 0.0008
bmin = 0.000001

# S.O. surface boundary conditions and grid:
l = 2.e6
y = np.asarray(np.linspace(0, l, 40))
h = 50.    # SO ML depth. This is an important tuning parameter since it controls diffusve b exchange between fixed flux and restoring region
#Although the model can emulate the effect of a meridionally varying wind-stress,
# that feature is in "beta" and we are here for simplicity using a constant
# wind-stress (approximately the average over the channel in NV12):
tau = 0.13
# from what I can infer, NV12 use an approximately quadratic surface temperature
# profile in the channel for their SAMBUCA simulations:
bs_SO = (bs-bmin) * (y / y[-1])**2 + bmin

A_basin = 6e13    #area of the basin
A_north = A_basin / 50.    #area of northern sinking region
A_east = A_basin * 2.

# time-stepping parameters:
dt = 86400 * 30    # time-step for vert. adv. diff. calc.
MOC_up_iters = int(
    np.floor(2. * 360 * 86400 / dt)
)    # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters = int(
    np.ceil(300 * 360 * 86400 / dt)
)    # plotting frequency (in iterations)
total_iters = int(
    np.ceil(3000 * 360 * 86400 / dt)
)    # total number of timesteps

kappa = 2e-5

# create vertical grid:
z = np.asarray(np.linspace(-4000, 0, 80))

# Initial conditions for buoyancy profile in the basin and north
def b_basin(z):
  return bs * np.exp(z / 300.)
def b_north(z):
  return bs_north * np.exp(z / 300.)
def b_east(z):
  return bs_east * np.exp(z / 300.)

# Here we are creating the model object, to which we then add the
# various components ("modules")
model = model.Model()
# Add module for SO overturning streamfunction:
model.new_module(
    Psi_SO, {
        'z': z,
        'y': y,
        'b': b_basin(z),
        'bs': bs_SO,
        'tau': tau,
        'f': 1e-4,
        'L': 5e6,
        'KGM': 1000.,
        'c': 0.1,
        'bvp_with_Ek': False
    }, 'Psi SO Atl'
)
# Add column module for adv-diff "Atlantic basin"
model.new_module(
    Column,
    {
        'z': z,
        'kappa': kappa,
        'Area': A_basin,
        'b': b_basin,
        'bs': bs,
        'bbot': bmin
    },
    'Atlantic Basin',
    left_neighbors=[model.get_module('psi_so_atl')],
)
# Add module for N.A. overturning
model.new_module(
    Psi_Thermwind,
    {
        'z': z,
        'b1': b_basin,
        'b2': b_north,
        'f': 1e-4
    },
    'AMOC',
    left_neighbors=[model.get_module('atlantic_basin')],
)

# Add column module for northern deep water formation region
model.new_module(
    Column,
    {
        'z': z,
        'kappa': kappa,
        'Area': A_north,
        'b': b_north,
        'bs': bs_north,
        'bbot': bmin,
        'do_conv': True
    },
    'North Atlantic',
    left_neighbors=[model.get_module('amoc')],
)
# Add module for Pacific overturning
model.new_module(
    Psi_SO, {
        'z': z,
        'y': y,
        'b': b_east(z),
        'bs': bs_SO,
        'tau': tau,
        'f': 1e-4,
        'L': 5e6,
        'KGM': 1000.,
        'c': 0.1,
        'bvp_with_Ek': False
    }, 'Psi SO EPac'
)
model.new_module(
    Psi_Thermwind,
    {
        'z': z,
        'b1': b_east,
        'b2': b_basin,
        'f': 1e-4
    },
    'ZOC',
    right_neighbors=[model.get_module('atlantic_basin')],
)
model.new_module(
    Column,
    {
        'z': z,
        'kappa': kappa,
        'Area': A_east,
        'b': b_east,
        'bs': bs_east,
        'bbot': bmin,
        'do_conv': True
    },
    'East Pacific',
    right_neighbors=[model.get_module('zoc')],
    left_neighbors=[model.get_module('psi_so_epac')],
)

# Add column module for northern deep water formation region
# model.new_module(
#     Column,
#     {
#         'z': z,
#         'kappa': kappa,
#         'Area': A_north,
#         'b': b_north,
#         'bs': bs_north,
#         'bbot': bmin,
#         'do_conv': True
#     },
#     'North Atlantic',
#     left_neighbors=[model.get_module('amoc')],
# )

L=4e6
kaps = 400.    # horizontal ML diffusivity is SO (m^2/s)
surflux = 0. * y
Bloss = 5.9e3 / L / 2e5    #buoyancy loss around Antarctica - notice that the argument is for the globally integrated b-loss rate
surflux[1:6] = -Bloss
rest_mask = 0. * y
rest_mask[6:-1] = 1.
bminSO = 0.0 
#bs_SO_eq = 0.*y + bminSO
bs_SO_eq = bs_SO
vpist = 1.5 / 86400.    # Piston Vel for SO restoring

model.new_module(
    SO_ML,
    {
       'y': y,
       'h': h,
       'L': L,
       'Ks': kaps,
       'surflux': surflux,
       'rest_mask': rest_mask,
       'b_rest': bs_SO_eq,
       'v_pist': vpist,
       'bs': bs_SO
    },
    'Southern Ocean',
    right_neighbors=[model.get_module('psi_so_atl'), model.get_module('psi_so_epac')],
)

fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3, 0))
ax1.set_xlim((-10, 15))
ax2.set_xlim((-0.02, 0.03))
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)
model.run(basin_dt=dt, coupler_dt=MOC_up_iters, steps=total_iters)
# # Plot final results over time-iteration plot:
ax1.plot(model.amoc.Psi, model.amoc.z, linewidth=2, color='r')
ax1.plot(model.psi_so_atl.Psi, model.psi_so_atl.z, linewidth=2, color='m')
ax1.plot(
    model.psi_so_atl.Psi_Ek,
    model.psi_so_atl.z,
    linewidth=1,
    color='m',
    linestyle='--'
)
ax1.plot(
    model.psi_so_atl.Psi_GM, model.psi_so_atl.z, linewidth=1, color='m', linestyle=':'
)
ax2.plot(
    model.atlantic_basin.b, model.atlantic_basin.z, linewidth=2, color='b'
)
ax2.plot(
    model.north_atlantic.b, model.atlantic_basin.z, linewidth=2, color='c'
)
ax1.plot(0. * model.amoc.z, model.amoc.z, linewidth=0.5, color='k')

# #==============================================================================
# # Everything below is just to make fancy 2D plots:
# # This is mostly an exercise in interpolation, but useful to visualize solutions

# # first interpolate buoyancy in channel along constant-slope isopycnals:
bint = Interpolate_channel(y=y, z=z, bs=bs_SO, bn=model.atlantic_basin.b)
bsouth = bint.gridit()
# buoyancy in the basin is all the same:
ybasin = np.linspace(200., 12000., 60)
bbasin = np.tile(model.atlantic_basin.b, (len(ybasin), 1))
# finally, interpolate buoyancy in the north:
ynorth = np.linspace(12100., 13000., 10)
bnorth = np.zeros((len(ynorth), len(z)))
for iz in range(len(z)):
  bnorth[:, iz] = interp1d([11900., 12000., 12900., 13000.], [
      model.atlantic_basin.b[iz], model.atlantic_basin.b[iz],
      model.north_atlantic.b[iz], model.north_atlantic.b[iz]
  ],
                           kind='quadratic')(ynorth)
# now stick it all together:
ynew = np.concatenate(((y-l) / 1e3, ybasin, ynorth))
bnew = np.concatenate((bsouth, bbasin, bnorth))

# Evaluate isopycnal overturning streamfunction at all latitudes:
psiarray_iso = np.zeros((len(ynew), len(z)))
psiarray_z = np.zeros((len(ynew), len(z)))
for iy in range(1, len(y)):
  # in the channel, interpolate model.psi_so.Psi onto local isopycnal depth:
  psiarray_iso[iy, :] = np.interp(
      bnew[iy, :], model.atlantic_basin.b, model.psi_so_atl.Psi
  )
  psiarray_z[iy, :] = psiarray_iso[iy, :]
for iy in range(len(y), len(y) + len(ybasin)):
  # in the basin, linearly interpolate between Psi_SO and Psi_AMOC:
  psiarray_iso[iy, :] = (
      ynew[iy] * model.amoc.Psibz()[0] + (10000. - ynew[iy]) * model.psi_so_atl.Psi
  ) / 10000.
  psiarray_z[iy, :] = (
      ynew[iy] * model.amoc.Psi + (10000. - ynew[iy]) * model.psi_so_atl.Psi
  ) / 10000.
for iy in range(len(y) + len(ybasin), len(ynew)):
  # in the north, interpolate AMOC.psib to local isopycnal depth:
  psiarray_iso[iy, :] = np.interp(
      bnew[iy, :], model.amoc.bgrid, model.amoc.Psib()
  )
  psiarray_z[iy, :] = ((13000. - ynew[iy]) * model.amoc.Psi) / 1000.
psiarray_iso[-1, :] = 0.

blevs = np.arange(0.001, 0.03, 0.002)
# Notice that the contour level for Psi is 0.5 SV for negative values
# and 1 SV for positive values. This appers to be what was used in NV2012.
plevs = np.concatenate((np.arange(-15., -0.4, 0.5), np.arange(1., 16., 1.)))

# plot isopycnal overturning:
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(111)
CS = ax1.contour(
    ynew,
    z,
    bnew.transpose(),
    levels=blevs,
    colors='k',
    linewidths=1.0,
    linestyles='solid'
)
ax1.clabel(CS, fontsize=10)
ax1.contour(
    ynew,
    z,
    psiarray_iso.transpose(),
    levels=plevs,
    colors='k',
    linewidths=0.5
)
CS = ax1.contourf(
    ynew,
    z,
    psiarray_iso.transpose(),
    levels=plevs,
    cmap=plt.cm.bwr,
    vmin=-15,
    vmax=15
)
#fig.colorbar(CS, ticks=plevs[0::2], orientation='vertical')
ax1.set_xlabel('y [km]')
ax1.set_ylabel('Depth [m]')
ax1.set_title('Isopycnal Overturning')

# plot z-coord. overturning:
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(111)
CS = ax1.contour(
    ynew,
    z,
    bnew.transpose(),
    levels=blevs,
    colors='k',
    linewidths=1.0,
    linestyles='solid'
)
ax1.clabel(CS, fontsize=10)
ax1.contour(
    ynew, z, psiarray_z.transpose(), levels=plevs, colors='k', linewidths=0.5
)
CS = ax1.contourf(
    ynew,
    z,
    psiarray_z.transpose(),
    levels=plevs,
    cmap=plt.cm.bwr,
    vmin=-15,
    vmax=15
)
#fig.colorbar(CS, ticks=plevs[0::2], orientation='vertical')
ax1.set_xlabel('y [km]')
ax1.set_ylabel('Depth [m]')
ax1.set_title('Depth-averaged Overturning')
plt.savefig('test.png')
# plt.show()
