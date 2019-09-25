import numpy as np
from scipy import integrate
from pymoc.utils import make_func, make_array, check_numpy_version


class Column(object):
  r"""
  .. module:: Column

  :synopsis: Vertical Advection-Diffusion Column Model

  Instances of this class represent 1D representations of buoyancy in a water
  column governed by vertical advection and diffusion. The velocity
  profile is required as an input. The script can either compute the equilibrium
  buoyancy profile for a given vertical velocity profile and boundary conditions
  or compute the tendency and perform a time-step of given length.
  BCs have to be fixed buoyancy at the top and either fixed b or db/dz at the bottom
  The time-stepping version can also handle horizontal advection
  into the column. This is, however, not (yet) implemented for the equilibrium solver
  """

  # This module creates an advective-diffusive column
  # Notice that the column here represents a horizontal integral, rather than
  # an average, thus allowing for the area of to be a function of depth
  def __init__(
      self,
      z=None,    # grid (input)
      kappa=None,    # diffusivity profile (input)
      bs=0.025,    # surface buoyancy bound. cond (input)
      bbot=0.0,    # bottom buoyancy boundary condition (input)  
      bzbot=None,    # bottom strat. as alternative boundary condition (input) 
      b=0.0,    # Buoyancy profile (input, output)
      Area=None,    # Horizontal area (can be function of depth)
      N2min=1e-7    # Minimum strat. for conv adjustment
  ):
    """
    Parameters
    ----------

    z : ndarray; input
        Vertical depth levels of column grid. Units: m
    kappa : number, function, or ndarray; input
            Vertical diffusivity profile. Units: m\ :sup:`2`/s
    bs : number; input
         Surface level buoyancy boundary condition. Units: m/s\ :sup:`2`
    bbot : number; optional; input
           Bottom level buoyancy boundary condition. Units: m/s\ :sup:`2`
    bzbot : number; optional; input
            Bottom level buoyancy stratification. Can be used as an alternative to **bbot**. Units: s\ :sup:`-2`
    b : number, function, or ndarray; input, output
        Initial vertical buoyancy profile. Recalculated on model run. Units: m/s
    Area : number, function, or ndarray; input
           Horizontal area of basin. Units: m\ :sup:`2`
    N2min : number; optional; input
            Minimum stratification for convective adjustment. Units: s\ :sup:`-1`
    """

    # initialize grid:
    if isinstance(z, np.ndarray) and len(z) > 0:
      self.z = z
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    self.kappa = make_func(kappa, self.z, 'kappa')
    self.Area = make_func(Area, self.z, 'Area')

    self.bs = bs
    self.bbot = bbot
    self.bzbot = bzbot

    self.N2min = N2min

    self.b = make_array(b, self.z, 'b')

    if check_numpy_version():
      self.bz = np.gradient(self.b, z)
    else:
      self.bz = 0. * z    # notice that this is just for initialization of ode solver

  def Akappa(self, z):
    r"""
    Compute the area integrated diffusivity :math:`A\kappa`
    at depth(s) z.

    Parameters
    ----------

    z : number or ndarray; input
        Vertical depth level(s) at which to retrieve the integrated diffusivity.
    """

    return self.Area(z) * self.kappa(z)

  def dAkappa_dz(self, z):
    r"""
    Compute the area integrated diffusivity gradient :math:`\partial_z\left(A\kappa\right)`
    at depth(s) z.

    Parameters
    ----------

    z : number or ndarray; input
        Vertical depth level(s) at which to retrieve the integrated diffusivity gradient.
    """

    if not check_numpy_version():
      raise ImportError(
          'You need NumPy version 1.13.0 or later. Please upgrade your NumPy libary.'
      )
    return np.gradient(self.Akappa(z), z)

  def bc(self, ya, yb):
    r"""
    Calculate the residuals oof boundary conditions for the advective-diffusive
    boundary value problem.

    Parameters
    ----------

    ya : ndarray; input
         Bottom boundary condition. Units: m/s\ :sup:`2`
    yb : ndarray; input
         Surface boundary condition. Units: m/s\ :sup:`2`
    """

    #return the boundary conditions
    if self.bzbot is None:
      return np.array([ya[0] - self.bbot, yb[0] - self.bs])
    else:
      return np.array([ya[1] - self.bzbot, yb[0] - self.bs])

  def ode(self, z, y):
    r"""
    Generate the ordinary differential equation for the equilibrium buoyancy profile,
    to be solved as a boundary value problem:

    :math:`\partial_tb\left(z\right)=-w^\dagger\partial_zb+\partial_z\left(\kappa_{e\!f\!f}\partial_zb\right)`

    Parameters
    ----------

    z : ndarray; input
        Vertical depth levels of column grid on which to solve the ode. Units: m
    y : ndarray; input
        Initial values for buoyancy and buoyancy gradient profiles.
    """

    #return the equation to be solved
    return np.vstack(
        (y[1], (self.wA(z) - self.dAkappa_dz(z)) / self.Akappa(z) * y[1])
    )

  def solve_equi(self, wA):
    r"""
    Solve for the equilibrium buoyancy profile, given a specified vertical
    velocity profile, and pre-set surface and bottom boundary conditions.

    Parameters
    ----------

    wA : ndarray; input
         Area integrated velocity profile for the equilibrium solution. Units: m\ :sup:`3`/s
    """

    #Solve for equilibrium solution given vert. vel profile and BCs
    self.wA = make_func(wA, self.z, 'w')
    sol_init = np.zeros((2, np.size(self.z)))
    sol_init[0, :] = self.b
    sol_init[1, :] = self.bz
    res = integrate.solve_bvp(self.ode, self.bc, self.z, sol_init)
    # interpolate solution for b and db/dz onto original grid
    self.b = res.sol(self.z)[0, :]
    self.bz = res.sol(self.z)[1, :]

  def vertadvdiff(self, wA, dt):
    r"""
    Calculate and apply the upwind forcing from advection and diffusion on the vertical buoyancy
    profile, for the timestepping solution.

    Parameters
    ----------

    wA : number or ndarray; input
         Area integrated velocity profile for the timestepping solution. Units: m\ :sup:`3`/s
    dt : number; input
         Numerical timestep over which solution are iterated. Units: s
    """

    #upwind vert. adv. and diffusion
    wA = make_array(wA, self.z, 'wA')
    dz = self.z[1:] - self.z[:-1]

    # apply boundary conditions:
    self.b[-1] = self.bs
    self.b[0] = self.bbot if self.bzbot is None else self.b[
        1] - self.bzbot * dz[0]

    bz = (self.b[1:] - self.b[:-1]) / dz
    bz_up = bz[1:]
    bz_down = bz[:-1]
    bzz = (bz_up-bz_down) / (0.5 * (dz[1:] + dz[:-1]))

    #upwind advection:
    weff = wA - self.dAkappa_dz(self.z)
    bz = bz_down
    bz[weff[1:-1] < 0] = bz_up[weff[1:-1] < 0]

    db_dt = (
        -weff[1:-1] * bz / self.Area(self.z[1:-1]) +
        self.kappa(self.z[1:-1]) * bzz
    )
    self.b[1:-1] = self.b[1:-1] + dt*db_dt

  def convect(self):
    r"""
    Carry out downward convective adustment of the vertical buoyancy profile to
    the minimum stratification, N2min. This adjustment assumes a fixed surface 
    buoyancy boundary condition.
    """
    # do convective adjustment to minimum strat N2min
    # notice that this parameterization currently only handles convection
    # from the top down which is the only case we really encounter here...
    # it also assumes a fixed surface buoyancy BC, and hence any
    # surface heat flux that's required to adjust the buoyancy of the
    # convecting column:
    ind = self.b > self.b[-1] + self.N2min * self.z
    self.b[ind] = self.b[-1] + self.N2min * self.z[ind]
    # Below is an energy conserving version that could be used
    # for model formulations without fixed surface b. But for fixed surface
    # b, the simpler version above is preferable as it deals better with long time steps
    # (for infinitesimal time-step and fixed surface b, the two are equivalent)
    # this version does currently also not include adjustment to finite strat.
    # dz=self.z[1:]-self.z[:-1];
    # dz=np.append(dz,dz[-1]);dz=np.insert(dz,0,dz[0])
    # dz=0.5*(dz[1:]+dz[0:-1]);
    # self.b[ind]=(np.mean(self.b[ind]*dz[ind]*self.Area(self.z[ind]))
    #            /np.mean(dz[ind]*self.Area(self.z[ind])) )

  def horadv(self, vdx_in, b_in, dt):
    # upwind horizontal advection:
    # notice that vdx_in is the total transport per unit height
    # into the column (units m^2/s, sign positive for velocity into the column)
    vdx_in = make_array(vdx_in, self.z, 'vdx_in')
    b_in = make_array(b_in, self.z, 'b_in')

    adv_idx = vdx_in > 0.0
    db = b_in - self.b

    self.b[adv_idx] = self.b[adv_idx] + dt * vdx_in[adv_idx] * db[
        adv_idx] / self.Area(self.z[adv_idx])

  def timestep(self, wA=0., dt=1., do_conv=False, vdx_in=None, b_in=None):
    #Integrate buoyancy profile evolution for one time-step
    # do vertical advection and diffusion:
    self.vertadvdiff(wA=wA, dt=dt)
    if vdx_in is not None:
      # do horizontal advection: (optional)
      if b_in is not None:
        self.horadv(vdx_in=vdx_in, b_in=b_in, dt=dt)
      else:
        raise TypeError('b_in is needed if vdx_in is provided')
    if do_conv:
      # do convection: (optional)
      self.convect()
