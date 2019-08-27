import numpy as np


def make_func(myst, axis, name):
  # turn mysterious object into callable function (if needed)
  if callable(myst):
    return myst
  elif isinstance(myst, np.ndarray):

    def funfun(x):
      return np.interp(x, axis, myst)

    return funfun
  elif isinstance(myst, float):

    def funfun(x):
      return myst + 0*x

    return funfun
  else:
    raise TypeError(name, 'needs to be either function, numpy array, or float')
