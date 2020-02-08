'''
This script performes the simulations used for the manuscript
Jansen and Nadeau (submitted to JPO).
The setup assumes a single basin, represented by two columns, and
a Southern Ocean channel region. The first column represents the basin interior,
while the second column represents the northern sinking region. An additional
1D advection-diffusion equation is solved at the surface mixed layer in the SO.
The regions are connected using the diagnostic relationships for the overturning 
circulation computed by the "psi_thermwind" and "psi_SO" modules.

This script is designed to be executed from the command line and takes various
options to set up the different experiments.
'''

import sys
from pymoc.modules import Psi_Thermwind, Psi_SO, SO_ML, Column
import numpy as np
import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--fixbSO', action='store_true')
  parser.add_argument('--fixPsiN', action='store_true')
  parser.add_argument('--fixPsiSO', action='store_true')
  parser.add_argument('--adiabatic', action='store_true')
  parser.add_argument('--db', type=float, default=0.0)
  parser.add_argument('--B', type=float, default=5.9e3)
  parser.add_argument('--pickup', default=None)
  parser.add_argument('--diagfile', default='diags.npz')
  parser.add_argument('--pickup_save_file', default=None)
  args = parser.parse_args()

  # boundary conditions:
  deltabs = args.db    # the amount of surface buoayancy change imposed at t=0
  bs = 0.02 + deltabs
  bs_north = -0.001 + deltabs
  bminSO = 0.0 + deltabs
  h = 50.    # SO ML depth. This is an important tuning parameter since it controls diffusve b exchange between fixed flux and restoring region
  L = 4e6    # zonal length of SO
  Bloss = args.B / L / 2e5    #buoyancy loss around Antarctica - notice that the argument is for the globally integrated b-loss rate
  kaps = 400.    # horizontal ML diffusivity is SO (m^2/s)
  kapGM = 800.    #GM diffusivity in SO (m^2/s)
  vpist = 1.5 / 86400.    # Piston Vel for SO restoring
  tau = 0.12
  if args.adiabatic:
    kaps = 0.25 * kaps
    L = 3. * L
    Bloss = Bloss / 3.    # We want to keep total b-loss fixed, so b-loss per unit area is reduced

  fixbSO = args.fixbSO
  fixpsiN = args.fixPsiN
  fixpsiSO = args.fixPsiSO

  if fixbSO:
    print('bSO fixed')
  if fixpsiN:
    print('psiN fixed')
  if fixpsiSO:
    print('psiSO fixed')

  if args.pickup is not None:
    pickup = np.load(args.pickup)
  else:
    pickup = None

  pickup_save_file = args.pickup_save_file

  diagfile = args.diagfile

  # S.O. surface boundary conditions and grid:
  l = 2.e6
  y = np.asarray(np.linspace(0, l, 51))
  dy = y[1]

  # target ("equilibrium") surface buoyancy profile for restoring
  bs_SO_eq = 0.*y + bminSO
  alpha = (1. - np.cos(np.pi * (l - y[5]) / 7.4e6))
  bs_SO_eq[6:] = (bs-bminSO) * (
      1. - np.cos(np.pi * (y[6:] - y[5]) / 7.4e6)
  ) / alpha + bminSO

  # fixed flux and restoring mask for SO ML:
  surflux = 0. * y
  surflux[1:6] = -Bloss
  rest_mask = 0. * y
  rest_mask[6:-1] = 1.

  A_basin = 8e13    # area of the basin
  A_north = A_basin / 50.    # area of northern sinking region

  # time-stepping parameters:
  dt = 86400. * 30    # time-step for vert. adv. diff. calc.
  MOC_up_iters = int(
      np.floor(1. * 360. * 86400. / dt)
  )    # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
  total_iters = int(
      np.ceil(12000. * 360 * 86400. / dt)
  )    # total number of timesteps
  Diag_iters = 10 * MOC_up_iters    # multiplier for Diags - needs to be multiple of MOC_up_iters

  # diffusivity profile from GCM simulations:
  kapgcm = np.array([
      1.2e-4, 0.882e-4, 0.544e-4, 0.393e-4, 0.305e-4, 0.235e-4, 0.207e-4,
      0.210e-4, 0.213e-4, 0.216e-4, 0.220e-4, 0.226e-4, 0.247e-4, 0.316e-4,
      0.377e-4, 0.407e-4, 0.389e-4, 0.407e-4, 0.454e-4, 0.517e-4, 0.633e-4,
      0.757e-4, 0.899e-4, 1.056e-4, 1.246e-4, 1.584e-4, 1.884e-4, 2.053e-4,
      2.168e-4, 2.332e-4
  ])
  zgcm = -1e3 * np.array([
      0.0, 0.0200, 0.045, 0.075, 0.110, 0.150, 0.200, 0.260, 0.330, 0.410,
      0.500, 0.600, 0.720, 0.860, 1.020, 1.200, 1.400, 1.600, 1.800, 2.000,
      2.200, 2.400, 2.600, 2.800, 3.000, 3.200, 3.400, 3.600, 3.800, 4.000
  ])
  if args.adiabatic:
    kapgcm = 0.25 * kapgcm

  # generate functions for full and effective diffusivity profile:
  def kappa(z):
    return np.interp(
        -z, -zgcm, kapgcm
    )    #notice that coordinate has to be monotonically increasing, hence we use -z

  def kappaeff(z):    # effective diffusivity profile with tapering in BBL
    return np.interp(-z, -zgcm, kapgcm) * (
        1. - np.maximum(-4000. - z + 500., 0.) / 500.
    )**2    #notice that coordinate has to be monotonically increasing, hence we use -z

  # create vertical grid:
  z = np.asarray(np.linspace(-4000., 0., 81))

  # number of vertical buoyancy levels for isopycnal computations
  nb = 500

  # Set initial conditions
  if pickup is not None:
    b_basin = 1.0 * pickup['arr_0']
    b_north = 1.0 * pickup['arr_1']
    bs_SO = 1.0 * pickup['arr_2']
  else:
    b_basin = bs * np.exp(z / 300.) + bs_north * z / z[0]
    b_north = bs_north * (z / z[0])**2.
    bs_SO = bs_SO_eq.copy()

  # create N.A. overturning model instance
  AMOC = Psi_Thermwind(z=z, b1=b_basin, b2=b_north, f=1.2e-4)
  # and solve for initial overturning streamfunction:
  AMOC.solve()

  # create S.O. overturning model instance
  PsiSO = Psi_SO(
      z=z, y=y, b=b_basin, bs=bs_SO, tau=tau, f=1.2e-4, L=L, KGM=kapGM
  )
  # and solve for initial overturning streamfunction:
  PsiSO.solve()

  # notice that the next few lines need to come after Psi_SO is computed
  # in case the latter is fixed to the pickup conditions
  bs_SO[-1] = bs
  # make sure buoyancy at northern end of channel (serves as BC for diffusion) equals surface b in basin
  if fixbSO:
    bs_SO[:-1
          ] = bs_SO[:-1] + deltabs    # instanteneously adjust surface b in SO

  # create vert. adv-diff column model instance for basin
  basin = Column(
      z=z, kappa=kappaeff, Area=A_basin, b=b_basin, bs=bs, bbot=b_basin[0]
  )
  # create adv-diff column model instance for basin
  north = Column(
      z=z,
      kappa=kappaeff,
      Area=A_north,
      b=b_north,
      bs=bs_north,
      bbot=b_north[0]
  )

  # create hor. adv-diff model instance for SO/channel ML
  channel = SO_ML(
      y=y,
      h=h,
      L=L,
      Ks=kaps,
      surflux=surflux,
      rest_mask=rest_mask,
      b_rest=bs_SO_eq,
      v_pist=vpist,
      bs=bs_SO
  )

  # Create empty arrays to save time-dependent diagnostics
  AMOC_save = np.zeros((len(z), int(total_iters / Diag_iters)))
  AMOC_b_save = np.zeros((nb, int(total_iters / Diag_iters)))
  bgrid_save = np.zeros((nb, int(total_iters / Diag_iters)))
  b_basin_save = np.zeros((len(z), int(total_iters / Diag_iters)))
  b_north_save = np.zeros((len(z), int(total_iters / Diag_iters)))
  bs_SO_save = np.zeros((len(y), int(total_iters / Diag_iters)))
  Psi_SO_save = np.zeros((len(z), int(total_iters / Diag_iters)))

  # *****************************************************************************
  for ii in range(0, total_iters):
    # Main time-stepping loop:

    if ii % MOC_up_iters == 0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(
          b1=basin.b, b2=north.b
      )    # Notice that, even with fixed PsiN,
      #buoyancy profiles need to be updated for isopycnal mapping
      if not fixpsiN:
        # update northern overturning
        AMOC.solve()
      [Psi_res_b, Psi_res_n] = AMOC.Psibz(nb=nb)    #isopycnal mapping
      PsiSO.update(b=basin.b, bs=channel.bs)
      if not fixpsiSO:
        # update SO overturning
        PsiSO.solve()
      if ii % Diag_iters == 0:
        # save diagnostics:
        AMOC_save[:, int(ii / Diag_iters)] = AMOC.Psi
        AMOC_b_save[:, int(ii / Diag_iters)] = AMOC.Psib(nb=nb)
        bgrid_save[:, int(ii / Diag_iters)] = AMOC.bgrid
        b_basin_save[:, int(ii / Diag_iters)] = basin.b
        b_north_save[:, int(ii / Diag_iters)] = north.b
        bs_SO_save[:, int(ii / Diag_iters)] = channel.bs
        Psi_SO_save[:, int(ii / Diag_iters)] = PsiSO.Psi

    # update residual vertical velocity in columns:
    wAb = (Psi_res_b - PsiSO.Psi) * 1e6
    wAN = -Psi_res_n * 1e6

    # adjust bottom boundary conditions and bbl kappa:
    if PsiSO.Psi[1] < 0:
      # bottom water coming in from the south
      basin.bbot = channel.bs[0]
      basin.kappa = kappaeff
    if Psi_res_b[1] > 0 and north.b[0] < basin.b[1] and north.b[
        0] < channel.bs[0]:
      # bottom water coming in from the north
      basin.bbot = north.b[0]
      basin.kappa = kappaeff
    elif PsiSO.Psi[
        1] >= 0:    # notice that this implies neither of the two above
      # no bottom water coming in - no flux BBC and flat isopycnals
      basin.bbot = basin.b[1]
      basin.kappa = kappa
    if Psi_res_n[1] < 0 and basin.b[0] < north.b[1]:
      # bottom water coming in from basin:
      north.bbot = basin.b[0]
      north.kappa = kappaeff
    else:
      # no bottom water coming in - no flux BBC and flat isopycnals
      north.bbot = north.b[1]
      north.kappa = kappa

    # time-step adv-diff equations in columns:
    basin.timestep(wA=wAb, dt=dt, do_conv=True)
    north.timestep(wA=wAN, dt=dt, do_conv=True)
    if not fixbSO:
      # time-step SO ML buoyancy
      channel.timestep(b_basin=basin.b, Psi_b=PsiSO.Psi, dt=dt)

  # **************** end of main time-stepping loop *************************

  # write-out pickup and diagnostics:
  if pickup_save_file is not None:
    np.savez(pickup_save_file, basin.b, north.b, channel.bs)
  if diagfile is not None:
    np.savez(
        diagfile, AMOC_save, AMOC_b_save, b_basin_save, b_north_save,
        bs_SO_save, z, bgrid_save, y, Psi_SO_save, tau, kapGM
    )
