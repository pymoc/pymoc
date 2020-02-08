'''
This script performes the simulations using a single global-scalebasin,
represented by two columns, and a Southern Ocean channel region. The setup is
similar to that of Jansen and Nadeau 2019, but parameteres have been adjusted
to be representative of the global ocean. The two columns represent
the interior of the basins and the northern North Atlannorthern sinking region.
An additional 1D advection-diffusion equation is solved at the surface mixed
layer in the SO. The regions are connected using the diagnostic relationships for
the overturning circulation computed by the "Psi_Thermwind" and "Psi_SO" modules.

This script is designed to be executed from the command line and takes various
optional inputs for parameters.
'''

from pymoc.modules import Psi_Thermwind, Psi_SO, SO_ML, Column
import numpy as np
import argparse

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--bs', type=float, default=0.025)
  parser.add_argument('--B', type=float, default=5.0e4)
  parser.add_argument('--Ks', type=float, default=1.0e3)
  parser.add_argument('--KGM', type=float, default=1.0e3)
  parser.add_argument('--tau', type=float, default=0.12)
  parser.add_argument('--kapfac', type=float, default=1.0)
  parser.add_argument('--Nyears', type=float, default=10000.0)
  parser.add_argument('--dt', type=float, default=30.0)
  parser.add_argument('--pickup', default=None)
  parser.add_argument('--diagfile', default='diags.npz')
  parser.add_argument('--pickup_save_file', default=None)
  args = parser.parse_args()


  # total run time in years:
  Nyears=args.Nyears;
  # Time step in days:
  dtdays=args.dt

  # boundary conditions:
  bs = args.bs
  bs_north = 0.0 
  bminSO = 0.0 
  h = 50.    # SO ML depth. This is an important tuning parameter since it controls diffusve b exchange between fixed flux and restoring region
  L = 2e7    # zonal length of SO
  Bloss = args.B / L / 2e5    #buoyancy loss around Antarctica - notice that the argument is for the globally integrated b-loss rate
  kaps = args.Ks    # horizontal ML diffusivity is SO (m^2/s)
  kapGM = args.KGM    #GM diffusivity in SO (m^2/s)
  vpist = 1.5 / 86400.    # Piston Vel for SO restoring
  tau = args.tau
  kapfac=args.kapfac

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

  A_basin = 3.2e14    # area of the basin (About the area of the Atlantic+Pacific+Indian)
  A_north = A_basin / 100.    # area of northern sinking region

  # time-stepping parameters:
  dt = 86400. * dtdays    # time-step for vert. adv. diff. calc.
  MOC_up_iters = int(
      np.floor(2. * 360. * 86400. / dt)
  )    # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
  total_iters = int(
      np.ceil(Nyears * 360 * 86400. / dt)
  )    # total number of timesteps
  Diag_iters = 10 * MOC_up_iters    # multiplier for Diags - needs to be multiple of MOC_up_iters

  
  # generate functions for full and effective diffusivity profile:
  def kappa(z):
    return kapfac*9e-6*np.exp(-z/1200) + 7e-5*np.exp(z/50)

  def kappaeff(z):    # effective diffusivity profile with tapering in BBL
    return (kapfac*9e-6*np.exp(-z/1200) + 7e-5*np.exp(z/50)) * (
        1. - np.maximum(-4500. - z + 500., 0.) / 500. )**2   

  # create vertical grid:
  z = np.asarray(np.linspace(-4500., 0., 46))

  # number of vertical buoyancy levels for isopycnal computations
  nb = 500

  # Set initial conditions
  if pickup is not None:
    b_basin = 1.0 * pickup['arr_0']
    b_north = 1.0 * pickup['arr_1']
    bs_SO = 1.0 * pickup['arr_2']
  else:
    b_basin = bs * np.exp(z / 400.) - 0.0001 * z / z[0]
    b_north = bs_north - 0.0001*(z / z[0])**2.
    bs_SO = bs_SO_eq.copy(); bs_SO[:6]= -0.0001

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
      # update northern overturning
      AMOC.solve()
      [Psi_res_b, Psi_res_n] = AMOC.Psibz(nb=nb)    #isopycnal mapping
      PsiSO.update(b=basin.b, bs=channel.bs)
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
