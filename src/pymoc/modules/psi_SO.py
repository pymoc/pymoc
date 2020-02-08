# coding: utf-8

import numpy as np
from scipy import integrate, optimize
from pymoc.utils import make_func


class Psi_SO(object):
  r"""
  Southern Ocean Overturning Transport Model

  Instances of this class represent a 1D model of the overturning transport in
  the Southern Ocean interior, calculated based on a density profile in the
  adjoining basin, and local surface buoyancy and surface wind stress in the SO.

  """
  def __init__(
      self,
      z=None,    # vertical grid (array, in)
      y=None,    # horizontal grid (array, in)  
      b=None,    # buoyancy profile at northern end of ACC (function, array or float, in)
      bs=None,    # surface buoyancy (function, array or float, in)
      tau=None,    # surface wind stress (function, array or float, in)  
      f=1.2e-4,    # Coriolis parameter (in)
      rho=1030,    # Density of sea water (in)
      L=1e7,    # Zonal length of the ACC (in)
      KGM=1e3,    # GM coefficient (in)
      c=None,    # phase speed for F2010 BVP smoother of GM streamfunction
      bvp_with_Ek=False,    # if true, apply boundary condition that Psi_GM=-Psi_EK in F2010 BVP smoother 
      Hsill=None,    # height (in m above ocean floor) of the "sill", where Psi_Ek is tapered
      HEk=None,    # depth of surface Ekman layer
      Htapertop=None,    # A quadratic tapering of the GM streamfunction at the surface
      Htaperbot=None,    # A quadratic tapering of the GM streamfunction at the bottom
      smax=0.01,    # maximum slope for clipping of GM streamfunction
  ):
    r"""
    Parameters
    ----------

    z : ndarray
        Vertical depth levels of overturning grid. Units: m
    y : ndarray
        Meridional overturning grid. Units: m
    b : float, function, or ndarray
        Vertical buoyancy profile from the adjoining basin, on the north
        side of the ACC. Units: m/s\ :sup:`2`
    bs : float, function, or ndarray
         Surface level buoyancy boundary condition. Can be a constant,
         or an array or function in y. Units: m/s\ :sup:`2`
    tau : float, function, or ndarray
          Surface wind stress. Can be a constant, or an array or function
          in y. Units: N/m\ :sup:`2`
    f : float
        Coriolis parameter. Units s\ :sup:`-1`
    rho : float
          Density of sea water for Boussinesq approximation. Units: kg/m\ :sup:`3`
    L : float
        Zonal length of the modeled ACC. Units: m
    KGM : float
          Gent & McWilliams (GM) eddy diffusivity coefficient. Units: 
    c : float
        Phase speed cutoff for smoothing when solving the GM boundary value problem. Units: m/s 
    bvp_with_Ek : logical
            Whether to enforce the boundary condition that Psi_GM=-Psi_Ek at the ocean
            surface and bottom when solving the boundary value problem for the GM streamfunction.
    Hsill : float
            Height above the bottom at which the Ekman streamfunction is tapered. Units: m
    Hek : float
          Depth of the surface Ekman layer. Units: m
    Htapertop : float
                Height of the quadratic surface tapering layer for the GM streamfunction. Units: m
    Htaperbot : float
                Height of the quadratic bottom tapering layer for the GM streamfunction. Units: m
    smax : float
           Maximum slope of the GM streamfunction, above which Psi_GM is clipped. Units: m\ :sup:`-1`
    """

    # initialize grid:
    if isinstance(z, np.ndarray):
      self.z = z
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    if isinstance(y, np.ndarray):
      self.y = y
    else:
      raise TypeError(
          'y needs to be numpy array providing horizontal grid (or boundaries) of ACC'
      )

    self.b = make_func(b, self.z, 'b')
    self.bs = make_func(bs, self.y, 'bs')
    self.tau = make_func(tau, self.y, 'tau')
    self.f = f
    self.rho = rho
    self.L = L
    self.KGM = KGM
    self.c = c
    self.bvp_with_Ek = bvp_with_Ek
    self.Hsill = Hsill
    self.HEk = HEk
    self.Htapertop = Htapertop
    self.Htaperbot = Htaperbot
    self.smax = smax

  def ys(self, b):
    r"""
    Inversion function of :math:`bs\left(y\right)`. This is equivalent to the outcopping
    latitude of the isopycnal of density class :math:`b`.

    Parameters
    ----------

    b: float
       The surface buoyancy value for whose meridional location is being calculated.

    Returns
    -------

    ys : float
        The meridional location at which the surface buoyancy bs is equal to the supplied
        buoyancy value b.

    """
    def func(y):
      return self.bs(y) - b

    if b < np.min(self.bs(self.y)):
      # if b is smaller minimum bs, isopycnals don't outcrop and get handled separately
      return self.y[0] - 1e3
    if b > self.bs(self.y[-1]):
      # if b is larger than bs at northern end, return northernmost point:
      return self.y[-1]
    else:
      # if b in range of bs return ys(b):
      # Notice that this inversion is well defined only if bs is monotonically
      # increasing (past minind). Should probably add a check to make sure
      # this is the case...
      minind = np.argmin(self.bs(self.y))
      return optimize.brentq(func, self.y[minind], self.y[-1])

  def calc_N2(self):
    r"""
    Calculate the buouyancy (Brunt-Väisällä) frequency profile for the Southern Ocean

    Returns
    -------

    N2 : function
         A depth dependent function that returns the buoyancy frequency :math:`N^2` for a given depth :math:`z`.

    """

    dz = self.z[1:] - self.z[:-1]
    N2 = np.zeros(np.size(self.z))
    b = self.b(self.z)

    N2[1:-1] = (b[2:] - b[:-2]) / (dz[1:] + dz[:-1])
    N2[0] = (b[1] - b[0]) / dz[0]
    N2[-1] = (b[-1] - b[-2]) / dz[-1]

    return make_func(N2, self.z, 'N2')

  def calc_bottom_taper(self, H, z):
    r"""
    Calculate the quadratic tapering profile relative to the ocean floor.

    Parameters
    ----------

    H : float
        Height above the bottom at which the streamfunction is tapered. Units: m
    z : ndarray
        Vertical depth levels of overturning grid. Units: m

    Returns
    -------

    bottom_taper : ndarray
                   An array containing weights (0.0-1.0) corresponding to how much of the 
                   streamfunction should remain after tapering.

    """

    if H is not None:
      return 1. - np.maximum(z[0] + H - z, 0.)**2. / H**2.
    return 1.

  def calc_top_taper(self, H, z, scalar=True):
    r"""
    Calculate the quadratic tapering profile relative to the ocean surface.

    Parameters
    ----------

    H : float
        Depth from the surface at which the streamfunction is tapered. Units: m
    z : ndarray
        Vertical depth levels of overturning grid. Units: m

    Returns
    -------

    top_taper : ndarray
                An array containing weights (0.0-1.0) corresponding to how much of the 
                streamfunction should remain after tapering.

    """
    if H is not None:
      return 1 - np.maximum(z + H, 0)**2. / H**2.
    elif scalar:
      return 1.
    else:
      taper = np.ones(np.size(z))
      taper[-1] = 0.
      return taper

  def calc_Ekman(self):
    r"""
    Compute the Ekman transport from the wind stress averaged from the
    northern boundary of the domain to the latitude of the northernmost
    outcropped isopycnal.

    .. math::
      \Psi_{Ek} = -\frac{\tau L_x}{\rho_of_{SO}}

    Returns
    -------

    Psi_Ek : ndarray
             An array representing the meridional average of the Ekman transport at
             each vertical level of the Southern Ocean model.

    """

    tau_ave = 0 * self.z
    for ii in range(0, np.size(self.z)):
      y0 = self.ys(self.b(self.z[ii]))    # outcrop latitude
      tau_ave[ii] = np.mean(self.tau(np.linspace(y0, self.y[-1], 100)))

    silltaper = self.calc_bottom_taper(self.Hsill, self.z)
    Ektaper = self.calc_top_taper(self.HEk, self.z, scalar=False)
    return tau_ave / self.f / self.rho * self.L * silltaper * Ektaper

  def bc_GM(self, ya, yb):
    r"""
    Calculate the residuals of boundary conditions for the eddy driven transport
    boundary value problem.

    Parameters
    ----------

    ya : ndarray
         Bottom boundary condition. Units: Sv
    yb : ndarray
         Surface boundary condition. Units: Sv

    Returns
    -------

    bc : ndarray
         If bvp_with_Ek is false, an array containing the bottom boundary condition
         ya[0] and surface boundary condition yb[0].
         If bvp_with_Ek is true, an array containing the residuals of the supplied
         boundary conditions and the value of the Ekman transport :math:`\Psi_{Ek}`
         at those boundaries.
         
    """

    if self.bvp_with_Ek:
      return np.array([
          ya[0] + self.Psi_Ek[0] * 1e6, yb[0] + self.Psi_Ek[-1] * 1e6
      ])
    else:
      return np.array([ya[0], yb[0]])

  def calc_GM(self):
    r"""
    Compute the eddy (Gent & Mcwilliams) transport based on the meridionally
    averaged isopycnal slope.

    .. math::
      \begin{aligned}
      \Psi_{GM} &= K_{GM}\cdot s \\
      s\left(b\right) &\equiv \frac{z_B\left(b\right)}{L_y - y_{SO}\left(b\right)}
      \end{aligned}

    Where :math:`z_B\left(b\right)` is the depth of isopycnals of density class :math:`b`
    in the adjoining basin, and :math:`y_{SO}\left(b\right)` is the outcropping latitude
    of isopycnals of density class :math:`b` in the Southern Ocean, available via
    :meth:`pymoc.modules.Psi_SO.ys`.

    Returns
    -------

    Psi_GM : ndarray
             An array representing the meridional average of the eddy transport at
             each vertical level of the Southern Ocean model.

    """

    dy_atz = 0 * self.z
    eps = 0.1    # minimum dy (in meters) (to avoid div. by 0)
    for ii in range(0, np.size(self.z)):
      dy_atz[ii] = max(self.y[-1] - self.ys(self.b(self.z[ii])), eps)
    bottaper = self.calc_bottom_taper(self.Htaperbot, self.z)
    toptaper = self.calc_top_taper(self.Htapertop, self.z)
    if self.c is not None:
      temp = make_func(
          self.KGM * self.z / dy_atz * self.L * toptaper * bottaper, self.z,
          'psiGM'
      )
      N2 = self.calc_N2()

      def ode(z, y):
        return np.vstack((y[1], N2(z) / self.c**2. * (y[0] - temp(z))))

      #Solve the boundary value problem
      res = integrate.solve_bvp(
          ode, self.bc_GM, self.z, np.zeros((2, np.size(self.z)))
      )
      # return solution interpolated onto original grid
      temp = res.sol(self.z)[0, :]
    else:
      temp = self.KGM * np.maximum(
          self.z / dy_atz, -self.smax
      ) * self.L * toptaper * bottaper
    # limit Psi_GM to -Psi_Ek on isopycnals that don't outcrop:
    idx = dy_atz > self.y[-1] - self.y[0]
    temp[idx] = np.maximum(temp[idx], -self.Psi_Ek[idx] * 1e6)
    return temp

  def solve(self):
    r"""
    Compute the residual overturning transport in the Southern Ocean.

    .. math::
      \Psi_{SO} = \Psi_{Ek} + \Psi_{GM}

    Returns
    -------

    Psi : ndarray
          An array representing the meridional average of the residual overturning
          transport at each vertical level of the Southern Ocean. This is a scaled
          sum of the Ekman and eddy transports.
    """

    self.Psi_Ek = self.calc_Ekman() / 1e6
    self.Psi_GM = self.calc_GM() / 1e6
    self.Psi = self.Psi_Ek + self.Psi_GM
    # Notice that the Psi at the bottom boundary is somewhat poorly defined,
    # and only used foir plotting purposes, for which it makes sense to simply set it to zero:
    self.Psi[0] = 0.

  def update(self, b=None, bs=None):
    r"""
    Update the vertical buoyancy profile and surface buoyancy, based on changes
    in the adjoining basin and/or in the surface boundary conditions.

    Parameters
    ----------

    b : float, function, or ndarray
        Vertical buoyancy profile from the adjoining basin, on the north
        side of the ACC. Units: m/s\ :sup:`2`
    bs : float, function, or ndarray
         Surface level buoyancy boundary condition. Can be a constant,
         or an array or function in y. Units: m/s\ :sup:`2`
    """

    if b is not None:
      self.b = make_func(b, self.z, 'b')
    if bs is not None:
      self.bs = make_func(bs, self.y, 'bs')
