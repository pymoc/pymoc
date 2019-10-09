'''
This script defines a model class that can be used to compute the
overturning circulation between two columns, given buoyancy profiles
in the columns.
 
The model assumes a thermal-wind based equation for the overturning circulation
as in Nikurashin and Vallis (2012):
d_{zz}(\Psi) = f^{-1} (b_{2} - b_{1})

This equation is solved subject to the boundary conditions that
\Psi(0) = \Psi(-H) = 0 
(these BCs are different to NV2012)

An upwind isopycnal mapping is used to compute the isopycnal overturning transport
'''

import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from pymoc.utils import make_func, make_array


class Psi_Thermwind(object):
  def __init__(
      self,
      f=1.2e-4,    # Coriolis parameter (input)
      z=None,    # grid (input)
      sol_init=None,    # Initial conditions for ODE solver (input)
      b1=None,    # Buoyancy in the basin (input, output)
      b2=0.,    # Buoyancy in the deep water formation region (input, output)
      Psi=None,    # Streamfunction (output) 
  ):

    self.f = f
    # initialize grid:
    if isinstance(z, np.ndarray):
      self.z = z
      nz = np.size(z)
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    self.b1 = self.make_func(b1, 'b1', self.z)
    self.b2 = self.make_func(b2, 'b2', self.z)

    # Set initial conditions for BVP solver
    if sol_init is None:
      self.sol_init = np.zeros((2, nz))
    else:
      self.sol_init = sol_init

  # end of init

  def make_func(self, myst, name, zin):
    return make_func(myst, zin, name)

  def make_array(self, myst, name):
    return make_array(myst, self.z, name)

  def bc(self, ya, yb):
    #return the boundary conditions
    return np.array([ya[0], yb[0]])

  def ode(self, z, y):
    #return the equation to be solved
    return np.vstack((y[1], 1. / self.f * (self.b2(z) - self.b1(z))))

  def solve(self):
    #Solve the boundary value problem
    # Note: The solution to this BVP is a relatively straightforward integral
    # it would probably be faster to just code it up that way.
    res = integrate.solve_bvp(self.ode, self.bc, self.z, self.sol_init)
    # interpolate solution for overturning circulation onto original grid (and change units to SV)
    self.Psi = res.sol(self.z)[0, :] / 1e6

  def Psib(self, nb=500):
    # map overturning into isopycnal space:
    b1 = self.make_array(self.b1, 'b1')
    b2 = self.make_array(self.b2, 'b2')
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
    # update buoyancy profiles
    if b1 is not None:
      self.b1 = self.make_func(b1, 'b1', self.z)
    if b2 is not None:
      self.b2 = self.make_func(b2, 'b2', self.z)
