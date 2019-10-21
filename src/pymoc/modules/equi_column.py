import numpy as np
from scipy import integrate
from pymoc.utils import check_numpy_version


class Equi_Column(object):
  r"""
  Equilibrium 1D Column Model

  Instances of this class represent a 1D column model for the MOC. The model is written in terms
  of a boundary value problem solving the the ODE:

  .. math::
    d_{zzzz}(\Psi_{N}) = (\kappa A)^{-1}(\Psi_{N} - \Psi_{SO} - A d_z(\kappa))d_{zzz}(\Psi_N)

  That is subject to the boundary conditions:

  .. math::
    \begin{aligned}
    (1)&\ \Psi_N(0) = 0 \\
    (2)&\ \Psi_N(-H) = 0  \\
    (3)&\ b(0)=-f \partial_{zz} \Psi_N (0) = b_s \\
    (4)&\ b(-H)=-f \partial_{zz} \Psi_N (H) = b_{bot} \\
       &\textrm{or} \\
       &\ d_z b(-H) = -f d_{zzz} \Psi_N (-H) =  (A \kappa(-H))^{-1} B_{int} \\
    \end{aligned}

  Where :math:`H` is the total depth of the upper cell, which can also be solved for with the additional BC
  that: :math:`d_z\psi_N(-H) = 0`

  The solution is found by non-dimensionalizing the equations using :math:`H` and math:`f` as length and time scales, 
  and solving between :math:`z^*=z/H=0..1`. Notice that :math:`H` then appears as a parameter in the equations.
  """
  def __init__(
      self,
      f=1.2e-4,    # Coriolis parameter (input)
      b_s=0.025,    # surface buoyancy (input)
      b_bot=None,    # bottom buoyancy (input)
      B_int=3e3,    # integrated downward buoyancy flux at the bottom of the upper cell (input)
      A=7e13,    # Area of the basin (input) 
      nz=100,    # minimum number of vert. layers for numerical solver (input)
      sol_init=None,    # Initial conditions for ODE solver (input)
      H_guess=1500.,    # guess for depth if solved for (input)
      kappa=6e-5,    # Diapycnal diffusivity (input; can be const., function, or array on grid given by z)
      dkappa_dz=None,    # Vertical derivative of diffusivity profile (input; function or nothing)
      psi_so=None,    # SO streamfunction (input; function or array on grid given by z)
      z=None,    # vertical grid for I/O (input / output) 
      H=None,    # depth of cell (input / output)
  ):
    r"""
    Parameters
    ----------

    f : float
        Coriolis parameter. Units s\ :sup:`-1`
    b_s : float
          Buoyancy at the surface of the column. Units:
    b_bot : float
            Buoyancy at the bottom of the column. Units:
    B_int : float
            Integrated downward buyancy flux at the bottom of the surface cell. Units:
    A : float
           Horizontal area of basin. Units: m\ :sup:`2`
    nz : int
         Number of levels in the non-dimensional vertical grid for the numerical solver.
    sol_init : ndarray
               Initial guess of the solution to the boundary value problem.
    H_guess : float
              Initial guess for the depth of the upper cell, if solving for it. Units: m
    kappa : float, function, or ndarray
            Vertical diffusivity profile. Units: m\ :sup:`2`/s
    dkappa_dz : float, function, or ndarray
                Vertical diffusivity gradient profile. Units: m/s
    psi_so : float, function, or ndarray
             Streamfunction in the adjoining Southern Ocean basin. Units:
    z : ndarray
        Vertical depth levels of column grid. Units: m
    H : float
        Depth of the upper cell, if specifying. Units: m
    """

    self.f = f
    self.A = A
    self.H = H
    self.H_guess = H_guess
    self.z = z
    self.zi = np.asarray(
        np.linspace(-1, 0, nz)
    )    # grid for initial conditions for solver

    if not callable(dkappa_dz) and (
        callable(kappa) or isinstance(kappa, np.ndarray)
    ) and not check_numpy_version():
      raise ImportError(
          'You need NumPy version 1.13.0 or later if you want to automatically compute dkappa_dz. Please upgrade your NumPy libary.'
      )

    self.kappa = self.init_kappa(kappa)
    self.dkappa_dz = self.init_dkappa_dz(kappa, dkappa_dz)
    self.init_psi_so(psi_so)
    self.init_b_boundaries(b_s, b_bot, B_int)
    self.sol_init = self.calc_sol_init(sol_init, nz, b_bot)

  def init_kappa(self, kappa):
    r"""
    Initialize the vertical diffusivity profile.

    Parameters
    ----------
    
    kappa : float, function, or ndarray
            Vertical diffusivity profile. Units: m\ :sup:`2`/s

    Returns
    -------
    kappa(z,H) : function
                 A function that returns the nondimensionalized vertical diffusivity profile :math:`\frac{\kappa}{H^2f}`
                 as a function of depth :math:`z` and upper cell depth :math:`H`.
    """
    if callable(kappa):
      return lambda z, H: kappa(z * H) / (
          H**2 * self.f
      )    # non-dimensionalize (incl. norm. of vertical coordinate)
    elif isinstance(kappa, np.ndarray):
      return lambda z, H: np.interp(z * H, self.z, kappa) / (H**2 * self.f)
    else:
      return lambda z, H: kappa / (H**2 * self.f)

  def init_dkappa_dz(self, kappa, dkappa_dz=None):
    r"""
    Initialize the vertical diffusivity gradient profile.

    Parameters
    ----------
    
    kappa : float, function, or ndarray
            Vertical diffusivity profile. Units: m\ :sup:`2`/s
    dkappa_dz : float, function, or ndarray
                Vertical diffusivity profile. Units: m/s

    Returns
    -------
    dkappa_dz(z,H) : function
                     A function that returns the nondimensionalized vertical diffusivity gradient profile
                     :math:`\frac{\partial_z\kappa}{H\cdot f}` as a function of depth :math:`z` and upper cell depth :math:`H`.
    """

    if callable(kappa) and callable(dkappa_dz):
      return lambda z, H: dkappa_dz(z * H) / (H * self.f)
    elif callable(kappa):
      return lambda z, H: np.gradient(kappa(z * H), z * H) / (H * self.f)
    elif isinstance(kappa, np.ndarray):
      dkappa_dz = np.gradient(kappa, self.z)
      return lambda z, H: np.interp(z * H, self.z, dkappa_dz) / (H * self.f)
    else:
      return lambda z, H: 0

  # Initialize Southern Ocean Streamfunction
  def init_psi_so(self, psi_so=None):
    r"""
    Initialize the nondimensionalized streamfunction passed from the Southern Ocean.

    Parameters
    ----------

    psi_so : float, function, or ndarray
             Streamfunction in the adjoining Southern Ocean basin. Units:

    Returns
    -------

    psi_so(z, H) : function
                   If psi_so is specified as a function or array in depth space, returns the
                   nondimensionalized Southern Ocean overturning streamfunction
                   :math:`\frac{\Psi_{SO}}{f\cdot H^3}` as a function of depth :math:`z` and
                   upper cell depth :math:`H`. Otherwise, a function that simply returns 0.

    """

    if callable(psi_so):
      self.psi_so = lambda z, H: psi_so(z * H) / (
          self.f * H**3
      )    # non-dimensionalize (incl. norm. of vertical coordinate)
    elif isinstance(psi_so, np.ndarray):
      self.psi_so = lambda z, H: np.interp(z * H, self.z, psi_so
                                           ) / (self.f * H**3)
    else:
      self.psi_so = lambda z, H: 0

  def calc_sol_init(self, sol_init, nz=None, b_bot=None):
    r"""
    Initialize the initial guess for the solutions to the system of ordinary differential equations
    defined by :meth:`pymoc.modules.Equi_Column.ode`

    Parameters
    ----------

    sol_init : ndarray
               User specified initial guess of the solution to the boundary value problem.
    nz : integer
         Number of levels in the non-dimensional vertical grid for the numerical solver.
    b_bot : float
            Buoyancy at the bottom of the column. Units:

    Returns
    -------

    sol_init : ndarray
               If a user specified sol_init is provided, return that value. Otherwise, a shape (4, nz) array,
               where the first column contains 1 at each row, the second and third columns contain 0 at each
               row, and the fourth column contains -100 at each row if b_bot is specified, otherwise the 
               stratification at the 1500m depth level provided by :meth:`pymoc.modules.Equi_Column.bz`.
    """

    if sol_init is not None:
      return sol_init

    if nz is None:
      nz = len(self.z)
    b_init = -100.0 if b_bot is not None else -self.bz(1500.)

    sol_init = np.zeros((4, nz))
    sol_init[0, :] = np.ones((nz))
    sol_init[2, :] = 0.0 * np.ones((nz))
    sol_init[3, :] = b_init * np.ones((nz))
    return sol_init

  # initialize non-dimensional surface buoyancy and bottom buoyancy
  # or abyssal buoyancy flux boundary condition
  def init_b_boundaries(self, b_s, b_bot=None, B_int=None):
    r"""
    Initialize and set the nondimensionalize surface buoyance and bottom buoyance, or the
    abyssal buoyancy flux boundary condition if no bottom buoyancy is specified.

    Parameters
    ----------

    b_s : float
          Buoyancy at the surface of the column. Units:
    b_bot : float
            Buoyancy at the bottom of the column. Units:
    B_int : float
            Integrated downward buyancy flux at the bottom of the surface cell. Units:

    """
    if b_bot is None and B_int is None:
      raise Exception(
          'You need to specify either b_bot or B_int for bottom boundary condition'
      )
    self.bs = -b_s / self.f**2
    if b_bot is not None:
      self.b_bot = -b_bot / self.f**2
    else:
      self.B_int = B_int

  def alpha(self, z, H):
    r"""
    Calculate the nondimensional factor:

    .. math::
       \alpha\left(z, H\right) = \frac{H^2}{A\cdot\kappa\left(z, H\right)}

    Parameters
    ----------

    z : ndarray
        Vertical depth levels of column grid. Units: m
    H : float
        Depth of the upper cell, if specifying. Units: m

    Returns
    -------
    alpha : float
            The value of the nondimentional factor :math:`\alpha` for the specified
            depth :math:`z` and upper cell depth :math:`H`.
    """

    return H**2 / (self.A * self.kappa(z, H))

  def bz(self, H):
    r"""
    Calculate the nondimenzionalized stratification at the bottom of the upper cell

    .. math::
       \partial_zb\approx\frac{B_{int}}{f^3\cdot H^2\cdot A \cdot\kappa\left(-1, H\right)}

    Parameters
    ----------

    H : float
        Depth of the upper cell, if specifying. Units: m

    Returns
    -------
  
    bz : float
         The value of the nondimensionalized stratification at the bottom of the cell with
         specified depth H.

    """

    return self.B_int / (self.f**3 * H**2 * self.A * self.kappa(-1, H))

  def bc(self, ya, yb, p=None):
    r"""
    Calculate the boundary conditions for the equilibrium column
    boundary value problem.

    Parameters
    ----------

    ya : ndarray
         Bottom boundary condition. Units:
    yb : ndarray
         Surface boundary condition. Units:
    p : ndarray
        Initial guess for the unknown parameter, H, if not specified at model initialization. Units: m

    Returns
    -------

    bc : ndarray
         An array containing the boundary conditions of the ODE:

         .. math::
           \begin{aligned}
            bc\left[0\right]&=\Psi_N(0) \\
            bc\left[1\right]&=\Psi_N(-H) \\
            bc\left[2\right]&=\begin{cases}
              \partial_{zz}\Psi_N\left(0\right) - \frac{b_{bot}}{H} & b_{bot} \textrm{ is specified}\\
              \partial_{zzz}\Psi_N\left(0\right) + d_zb\left(-H\right) & b_{bot} \textrm{ is unspecified}
            \end{cases}\\
            bc\left[3\right]&=\partial_{zz}\Psi_N\left(H\right)-\frac{b_s}{H} \\
           \end{aligned}

         If :math:`H` is not specified by the user, then there is an additional boundary condition  inserted:

         .. math::
           bc\left[2\right] = \partial_z\Psi_N\left(-H\right)

         and each subsequent row moves up an index.

    """
    try:
      bbot_set = getattr(self, 'b_bot', None) is not None
      y = np.array([ya[0], yb[0]])
      if self.H is None:
        y = np.append(y, ya[1])
      d = p[0] if self.H is None else self.H
      if bbot_set:
        y = np.append(y, ya[2] - self.b_bot / d)
      else:
        y = np.append(y, ya[3] + self.bz(d))
      y = np.append(y, yb[2] - self.bs / d)

      return y
    except TypeError:
      raise TypeError(
          'Must provide a p array if column does not have an H value'
      )

  def ode(self, z, y, p=None):
    r"""
    Generate the ordinary differential equation for the equilibrium column overturning streamfunction,
    to be solved as a boundary value problem:

    .. math::
      d_{zzzz}(\Psi_{N}) = (\kappa A)^{-1}(\Psi_{N} - \Psi_{SO} - A d_z(\kappa))d_{zzz}(\Psi_N)

    Parameters
    ----------

    z : ndarray
        Vertical depth levels of column grid on which to solve the ode. Units: m
    y : ndarray
        Initial values for the streamfunction and its vertical gradient.
    p : ndarray
        Initial guess for the unknown parameter, H, if not specified at model initialization. Units: m

    Returns  
    -------
    ode : ndarray
          A vertically oriented array, containing the system of linear equations:

          .. math::
            \begin{aligned}
            \partial_zy_1 &= y_2 \\
            \partial_zy_2 &= y_3 \\
            \partial_zy_3 &= y_4 \\
            \partial_zy_4 &= \alpha\left(z, H\right)\cdot y_4\left(y_1 - \Psi_{SO}\left(z, H\right)-\frac{A\cdot\partial_z\kappa\left(z, H\right)}{H^2}\right)
            \end{aligned}

    """

    if self.H is None and p is not None and len(p) > 0:
      H = p[0]
    elif self.H is not None:
      H = self.H
    else:
      raise TypeError(
          'Must provide a p array if column does not have an H value'
      )
    return np.vstack((
        y[1], y[2], y[3], self.alpha(z, H) * y[3] *
        (y[0] - self.psi_so(z, H) - self.A * self.dkappa_dz(z, H) / (H**2))
    ))

  def solve(self):
    r"""
    Solve for the thermal wind overturning streamfunction as a boundary value problem
    based on the system of equations defined in :meth:`pymoc.modules.Equi_Column.ode`.

    """

    if self.H is None:
      res = integrate.solve_bvp(
          self.ode, self.bc, self.zi, self.sol_init, p=[self.H_guess]
      )
      self.H = res.p[0]
    else:
      res = integrate.solve_bvp(
          self.ode, self.bc, self.zi, self.sol_init, p=None
      )

    # if self.z does not yet exist use mesh from solver:
    if self.z is None:
      self.z = res.x * self.H
      self.psi = res.y[0, :] * self.f * self.H**3 / 1e6
      self.b = -res.y[2, :] * self.f**2 * self.H
    # if self.z does exist, compute solution at those points (setting all points below z=-H to NaN):
    else:
      self.psi = res.sol(self.z / self.H)[0, :] * self.f * self.H**3 / 1e6
      self.psi[self.z < -self.H] = np.NaN
      self.b = -res.sol(self.z / self.H)[2, :] * self.f**2 * self.H
      self.b[self.z < -self.H] = np.NaN
