import numpy as np
from scipy import integrate, optimize
from pymoc.utils import make_func


class Psi_SO(object):
  r"""
  Southern Ocean Overturning Trasport Model

  A 1D model of the Southern Ocean overturning transport,
  calculated based on a density profile in the adjoining basin,
  and local surface buoyancy and surface wind stress in the SO.
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

    z : ndarray; input
        Vertical depth levels of overturning grid. Units: m
    y : ndarray; input
        Meridional overturning grid. Units: m
    b : number, function, or ndarray; input
        Vertical buoyancy profile from the adjoining basin, on the north
        side of the ACC. Units: m/s\ :sup:`2`
    bs : number, function, or ndarray; input
         Surface level buoyancy boundary condition. Can be a constant,
         or an array or function in y. Units: m/s\ :sup:`2`
    tau : number, function, or ndarray; input
          Surface wind stress. Can be a constant, or an array or function
          in y. Units: N/m\ :sup:`2`
    f : number; input
        Coriolis parameter. Units s\ :sup:`-1`
    rho : number; input
          Density of sea water for Boussinesq approximation. Units: kg/m\ :sup:`3`
    L : number; input
        Zonal length of the modeled ACC. Units: m
    KGM : number; input
          Gent & McWilliams (GM) eddy diffusivity coefficient. Units: 
    c : number; input
        Phase speed cutoff for smoothing when solving the GM boundary value problem. Units: m/s 
    bvp_with_Ek : logical; input
            Whether to enforce the boundary condition that Psi_GM=-Psi_Ek at the ocean
            surface and bottom when solving the boundary value problem for the GM streamfunction.
    Hsill : number; input
            Height above the bottom at which the Ekman streamfunction is tapered. Units: m
    Hek : number; input
          Depth of the surface Ekman layer. Units: m
    Htapertop : number; input
                Height of the quadratic surface tapering layer for the GM streamfunction. Units: m
    Htaperbot : number; input
                Height of the quadratic bottom tapering layer for the GM streamfunction. Units: m
    smax : number; input
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
    Inversion function of :math:`bs\left(y\right).

    Parameters
    ----------

    b: number; input
       The surface buoyancy value for whose meridional location is being calculated.

    Returns
    -------

    The meridional location at which the surface buoyancy bs is equal to the supplied buoyancy value b.
    """

    # inverse of bs(y)
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
    dz = self.z[1:] - self.z[:-1]
    N2 = np.zeros(np.size(self.z))
    b = self.b(self.z)

    N2[1:-1] = (b[2:] - b[:-2]) / (dz[1:] + dz[:-1])
    N2[0] = (b[1] - b[0]) / dz[0]
    N2[-1] = (b[-1] - b[-2]) / dz[-1]

    return make_func(N2, self.z, 'N2')

  def calc_bottom_taper(self, H, z):
    if H is not None:
      return 1. - np.maximum(z[0] + H - z, 0.)**2. / H**2.
    return 1.

  def calc_top_taper(self, H, z, scalar=True):
    if H is not None:
      return 1 - np.maximum(z + H, 0)**2. / H**2.
    elif scalar:
      return 1.
    else:
      taper = np.ones(np.size(z))
      taper[-1] = 0.
      return taper

  def calc_Ekman(self):
    # compute Ekman transport on z grid
    # based on average wind stress between outcrop and northern end of channel
    tau_ave = 0 * self.z
    for ii in range(0, np.size(self.z)):
      y0 = self.ys(self.b(self.z[ii]))    # outcrop latitude
      tau_ave[ii] = np.mean(self.tau(np.linspace(y0, self.y[-1], 100)))

    silltaper = self.calc_bottom_taper(self.Hsill, self.z)
    Ektaper = self.calc_top_taper(self.HEk, self.z, scalar=False)
    return tau_ave / self.f / self.rho * self.L * silltaper * Ektaper

  def bc_GM(self, ya, yb):
    if self.bvp_with_Ek:
      return np.array([
          ya[0] + self.Psi_Ek[0] * 1e6, yb[0] + self.Psi_Ek[-1] * 1e6
      ])
    else:
      return np.array([ya[0], yb[0]])

  def calc_GM(self):
    # compute GM ransport on z grid
    # based on average isopycnal slope
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
    #Solve for overturning circ.
    self.Psi_Ek = self.calc_Ekman() / 1e6
    self.Psi_GM = self.calc_GM() / 1e6
    self.Psi = self.Psi_Ek + self.Psi_GM
    # Notice that the bottom boundary is somewhat poorly defined, as it is
    # only for BC. To avoid random fluctuation in bottom Psi, we here simply
    # set it to value of last gridpoint above
    self.Psi[0] = self.Psi[1]

  def update(self, b=None, bs=None):
    # update buoyancy profiles
    if b is not None:
      self.b = make_func(b, self.z, 'b')
    if bs is not None:
      self.bs = make_func(bs, self.y, 'bs')
