'''
This script defines a model class that can be used to solve a 1D column model
for the MOC.
 
The model is written in terms of a boundary value problem
 
In dimensional units, the ODE is
d_{zzzz}(\Psi_{N}) = (\kappa A)^{-1}(\Psi_{N} - \Psi_{SO} - A d_z(\kappa))d_{zzz}(\Psi_N)

This equation is solved subject to the boundary conditions:
(1) \Psi_N(0) = 0
(2) \Psi_N(-H) = 0 
(3) b(0)=-f \partial_{zz} \Psi_N (0) = b_s
(4) either  b(-H)=-f \partial_{zz} \Psi_N (H) = b_bot
    or d_z b(-H) = -f d_{zzz} \Psi_N (-H) =  (A \kappa(-H))^{-1} B_{int}
Where H is the total depth of the upper cell, which can also be solved for
with the additional BC that:
(5) d_z\psi_N(-H) = 0

The solution is found by non-dimensionalizing the equations using H and f as length and time scales
The model is then solved between z^*=z/H=0..1
Notice that H then appears as a parameter in the equations.
5 boundary conditions are needed if we also want to solve for the parameter H
'''

import numpy as np
from scipy import integrate
from pymoc.utils import check_numpy_version


class Equi_Column(object):
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
      psi=None,    # streamfunction (output)
      b=None,    # streamfunction (output)
      H=None,    # depth of cell (input / output)
  ):

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
    self.sol_init = sol_init if sol_init is not None else self.calc_sol_init(
        sol_init, nz, b_bot
    )

  def init_kappa(self, kappa):
    # Initialize vertical diffusivity profile:
    if callable(kappa):
      return lambda z, H: kappa(z * H) / (
          H**2 * self.f
      )    # non-dimensionalize (incl. norm. of vertical coordinate)
    elif isinstance(kappa, np.ndarray):
      return lambda z, H: np.interp(z * H, self.z, kappa) / (H**2 * self.f)
    else:
      return lambda z, H: kappa / (H**2 * self.f)

  def init_dkappa_dz(self, kappa, dkappa_dz=None):
    if callable(kappa) and callable(dkappa_dz):
      return lambda z, H: dkappa_dz(z * H) / (H * self.f)
    elif callable(kappa):
      return lambda z, H: np.gradient(kappa(z * H), z * H) / (H * self.f)
    elif isinstance(kappa, np.ndarray):
      dkappa_dz = np.gradient(kappa, self.z)
      return lambda z, H: np.interp(z * H, self.z, dkappa_dz) / (H * self.f)
    else:
      return lambda z, H: 0

  def calc_sol_init(self, sol_init, nz=None, b_bot=None):
    # Set initial conditions for ODE solver
    if nz is None:
      nz = len(self.z)
    b_init = -100.0 if b_bot is not None else -self.bz(1500.)

    sol_init = np.zeros((4, nz))
    sol_init[0, :] = np.ones((nz))
    sol_init[2, :] = 0.0 * np.ones((nz))
    sol_init[3, :] = b_init * np.ones((nz))
    return sol_init

  # Initialize Southern Ocean Streamfunction
  def init_psi_so(self, psi_so=None):
    if callable(psi_so):
      self.psi_so = lambda z, H: psi_so(z * H) / (
          self.f * H**3
      )    # non-dimensionalize (incl. norm. of vertical coordinate)
    elif isinstance(psi_so, np.ndarray):
      self.psi_so = lambda z, H: np.interp(z * H, self.z, psi_so
                                           ) / (self.f * H**3)
    else:
      self.psi_so = lambda z, H: 0

  # initialize non-dimensional surface buoyancy and bottom buoyancy
  # or abyssal buoyancy flux boundary condition
  def init_b_boundaries(self, b_s, b_bot=None, B_int=None):
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
    #return factor on the RHS of ODE
    return H**2 / (self.A * self.kappa(z, H))

  def bz(self, H):
    # return the properly non-dimensionalized stratification at the bottom of the cell
    return self.B_int / (self.f**3 * H**2 * self.A * self.kappa(-1, H))

  def bc(self, ya, yb, p=None):
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
    #return the ODE to be solved
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
    #Solve the boundary value problem
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
