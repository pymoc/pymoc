'''
A class to interpolate the buoyancy field in the North along lines
with constant slope - only to make fancy plots
'''

import sys
import numpy as np
from scipy.optimize import brenth
from pymoc.utils import gridit, make_func


class Interpolate_twocol(object):
  def __init__(
      self,
      y=None,    # y-grid
      z=None,    # z-grid
      bs=None,    # buoyancy profile in the south
      bn=None,    # buoyancy profile in the north
  ):

    # initialize grid:
    if isinstance(y, np.ndarray):
      self.y = y
    else:
      raise TypeError('y needs to be numpy array providing grid levels')
    if isinstance(z, np.ndarray):
      self.z = z
    else:
      raise TypeError('z needs to be numpy array providing grid levels')

    self.bs = self.make_func(bs, 'bs', self.z)
    self.bn = self.make_func(bn, 'bn', self.z)

  def make_func(self, myst, name, xin):
    return make_func(myst, xin, name)

  def __call__(self, y, z):
    l = self.y[-1]
    bsurf = self.make_func(
        self.y / l * self.bn(0) + (1 - self.y / l) * self.bs(0), 'bsurf',
        self.y
    )
    if z == 0 and y == 0:
      # slope ill defined at (0,0); evaluate infinitesimally below the surface:
      z = -0.01
    if z == self.z[0]:
      # slope also potentially ill defined at z= -H; evaluatejust above the bottom:
      z = 0.9999 * self.z[0]

    def fint(x):
      # function to help determine slope at bottom of vent. region
      return self.bn(0) - self.bs(-x * l)

    def fup(x):
      # function to help determine slope in vent. region
      return self.bs(z - x*y) - bsurf(y - z/x)

    def fdeep(x):
      # function to help determine slope below vent. region
      return self.bs(z - x*y) - self.bn(z + x * (l-y))

    # first determine slope at bottom of vent. region here
    sbot = brenth(fint, 0., 1.)
    # than set slope for stuff above and below...
    if z > -sbot * (l-y):
      s = brenth(fup, 1e-10, 1.0)
    else:
      s = brenth(fdeep, -1.0, 1.0)
    return self.bs(z - s*y)

  def gridit(self):
    return gridit(self.y, self.z, self)
