import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from pymoc.utils import make_func, make_array


class Psi_Thermwind(object):
  r"""
  Thermal Wind Closure

  Instances of this class represent the overturning circulation between two columns,
  given buoyancy profiles in those columns.

  The model assumes a thermal-wind based equation for the overturning circulation
  as in Nikurashin and Vallis (2012):

  .. math::
    d_{zz}\left(\Psi\right) = f^{-1} (b_2 - b_1)

  This equation is solved subject to the boundary conditions:

  .. math::
    \Psi(0) = \Psi(-H) = 0 

  (these BCs are different than NV2012)

  An upwind isopycnal mapping is used to compute the isopycnal overturning transport.

  """
  def __init__(
      self,
      f=1.2e-4,    # Coriolis parameter (input)
      z=None,    # grid (input)
      sol_init=None,    # Initial conditions for ODE solver (input)
      b1=None,    # Buoyancy in the basin (input, output)
      b2=0.,    # Buoyancy in the deep water formation region (input, output)
  ):
    r"""
    Parameters
    ----------

    f : float
        Coriolis parameter. Units s\ :sup:`-1`
    z : ndarray
        Vertical depth levels of overturning grid. Units: m
    sol_init : ndarray; optional
               Initial guess at the solution to the thermal wind overturning streamfunction. Units: [...]
    b1 : float, function, or ndarray; optional
         Vertical buoyancy profile from the southern basin. Units: m/s\ :sup:`2`
    b2 : float, function, or ndarray; optional
         Vertical buoyancy profile from the northern basin, representing the
         deepwater formation region. Units: m/s\ :sup:`2`
    """

    self.f = f
    # initialize grid:
    if isinstance(z, np.ndarray):
      self.z = z
      nz = np.size(z)
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    self.b1 = make_func(b1, self.z, 'b1')
    self.b2 = make_func(b2, self.z, 'b2')

    # Set initial conditions for BVP solver
    if sol_init is None:
      self.sol_init = np.zeros((2, nz))
    else:
      self.sol_init = sol_init

  def bc(self, ya, yb):
    r"""
    Calculate the residuals of boundary conditions for the thermal wind closure
    boundary value problem.

    Parameters
    ----------

    ya : ndarray
         Bottom boundary condition. Units:
    yb : ndarray
         Surface boundary condition. Units:

    Returns
    -------

    bc : ndarray
         An array containing the value of the streamfunction at the top and bottom
         levels of the vertical grid.
    """

    return np.array([ya[0], yb[0]])

  def ode(self, z, y):
    r"""
    Generate the ordinary differential equation for the thermal wind overturning streamfunction,
    to be solved as a boundary value problem:

    :math:`d_{zz}\left(\Psi\right) = f^{-1} (b_2 - b_1)`

    Parameters
    ----------

    z : ndarray
        Vertical depth levels of column grid on which to solve the ode. Units: m
    y : ndarray
        Initial values for the streamfunction and its vertical gradient.

    Returns
    -------
    ode : ndarray
          A vertically oriented array, containing the system of linear equations:

          .. math::
            \begin{aligned}
            \partial_zy_1 &= y_2 \\
            \partial_zy_2 &= \frac{b_2\left(z\right) - b_1\left(z\right)}{f}
            \end{aligned}

    """

    return np.vstack((y[1], 1. / self.f * (self.b2(z) - self.b1(z))))

  def solve(self):
    r"""
    Solve for the thermal wind overturning streamfunction as a boundary value problem
    based on the system of equations defined in :meth:`pymoc.modules.Psi_Thermwind.ode`.
    """

    # Note: The solution to this BVP is a relatively straightforward integral
    # it would probably be faster to just code it up that way.
    res = integrate.solve_bvp(self.ode, self.bc, self.z, self.sol_init)
    # interpolate solution for overturning circulation onto original grid (and change units to SV)
    self.Psi = res.sol(self.z)[0, :] / 1e6

  def Psib(self, nb=500):
    r"""
    Remap the overturning streamfunction from physical depth space, into isopycnal
    space

    .. math::
      \Psi^b\left(b\right) = \int_{-H}^0 \partial_z\Psi\left(z\right)\mathcal{H}\left[b - b_{up}\left(z\right)\right]

    by computing upwind density classes

    .. math::
      \begin{aligned}
      b_{up}\left(z\right) = 
      \begin{cases} 
        b_N\left(z\right), & \partial_z\Psi\left(z\right)  > 0 \\
        b_B\left(z\right), & \partial_z\Psi\left(z\right)  < 0
      \end{cases}
      \end{aligned}

    where :math:`b_N\left(z\right)` is the density profiles in the northern region, :math:`b_B\left(z\right)` is
    the density profile in the southern basin, and :math:`\mathcal{H}` is the Heaviside step function.

    Parameters
    ----------
    nb : int; optional
         Number of upstream density classes into which the streamfunction is to be remapped.

    Returns
    -------
    psib : ndarray
           An array representing the values of the overturning streamfunction in each upwind density class.
    """
    # map overturning into isopycnal space:
    b1 = make_array(self.b1, self.z, 'b1')
    b2 = make_array(self.b2, self.z, 'b2')
    bmin = min(np.min(b1), np.min(b2))
    bmax = max(np.max(b1), np.max(b2))
    self.bgrid = np.linspace(bmin, bmax, nb)
    udydz = -(self.Psi[1:] - self.Psi[:-1])
    psib = 0. * self.bgrid
    bup_bot = b1[:-1].copy()
    bup_top = b1[1:].copy()
    idx = udydz < 0
    bup_bot[idx] = b2[:-1][idx]
    bup_top[idx] = b2[1:][idx]
    for i in range(0, len(self.bgrid)):
      mask = np.clip((bup_top - self.bgrid[i]) / (bup_top-bup_bot), 0., 1.)
      psib[i] = np.sum(mask * udydz)
    return psib

  def Psibz(self, nb=500):
    r"""
    Remap the overturning streamfunction onto the native isopycnal-depth space
    of the columns in the northern region and southern basin.

    Parameters
    ----------
    nb : int; optional
         Number of upstream density classes into which the streamfunction is to be remapped in the intermediate :meth:`pymoc.modules.Psi_Thermwind.Psib` step.

    Returns
    -------
    psibz : ndarray
            An array where the first element is an array representing the streamfunction at each depth level in the southern basin, and the second represents the same in the northern region.
    """
    # map isopycnal overturning back into isopycnal-depth space of each column
    psib = self.Psib(nb)
    # This does a linear interploation in b:
    return [
        np.interp(self.b1(self.z), self.bgrid, psib),
        np.interp(self.b2(self.z), self.bgrid, psib)
    ]
    # Ths instead first estimates the depth levels for the bgrid and then does linear interpolation in z
    # either has pros and cons depending on the situation...
    #z1_of_bgrid=np.interp(self.bgrid,self.b1(self.z),self.z)
    #z2_of_bgrid=np.interp(self.bgrid,self.b2(self.z),self.z)
    #return [np.interp(self.z,z1_of_bgrid,psib),np.interp(self.z,z2_of_bgrid,psib)]

  def update(self, b1=None, b2=None):
    r"""
    Update the vertical buoyancy profiles from the southern basin and northern region.

    Parameters
    ----------

    b1 : float, function, or ndarray; optional
         Vertical buoyancy profile from the southern basin. Units: m/s\ :sup:`2`
    b2 : float, function, or ndarray; optional
         Vertical buoyancy profile from the northern basin, representing the
         deepwater formation region. Units: m/s\ :sup:`2`
    """
    # update buoyancy profiles
    if b1 is not None:
      self.b1 = make_func(b1, self.z, 'b1')
    if b2 is not None:
      self.b2 = make_func(b2, self.z, 'b2')
