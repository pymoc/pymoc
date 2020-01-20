'''
An class to interpolate the buoyancy field in the SO along lines
with constant slope - only to mke fancy looking plots
'''

import sys
import numpy as np
from scipy.optimize import brenth
from pymoc.utils import gridit, make_func


class Interpolate_channel(object):
  def __init__(
      self,
      y=None,    # y-grid
      z=None,    # z-grid
      bs=None,    # surface buoyancy
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

    self.bs = self.make_func(bs, 'bs', self.y)
    self.bn = self.make_func(bn, 'bn', self.z)

  def make_func(self, myst, name, xin): # Seems unecessary to define a method that already exists identically as a function, no?
    return make_func(myst, xin, name)

  def __call__(self, y, z):
      l=self.y[-1]
      if y==l:
          # slope can be ill defined at y=l, and solution is trivial, so it makes sense to treat this separately
          return self.bn(z)
      else:
         def f2(x):
             # function to help determine slope at bottom of vent. region
             return self.bn(x)-self.bs(0)
         
         def f(x):
             # function to help determine slope in vent. region
             return self.bn(z - x * (l-y)) - self.bs(y + z/x)
         
         # first determine slope at bottom of vent. region here
         sbot=-brenth(f2, self.z[0],0.)/l
         # than set slope for stuff above and below...
         if -z>sbot*y:
           s=sbot
         else:
           s=brenth(f, 1.e-12,1.0)
         return self.bn(z - s * (l-y))


  def gridit(self):
    return gridit(self.y, self.z, self)
