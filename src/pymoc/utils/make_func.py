import numpy as np


def make_func(myst, axis, name):
  r"""
  Make an argument of unknown type into a function if needed, in the dimension of the specified array.

  Parameters
  ----------

  myst : float, function, or ndarray
         The argument to be transformed. Either a single numerical value, a function in a single dimension, or an array.
  axis : ndarray
         Grid points in the direction along which the function is to operate. If myst is an array, axis should be of the same length.
  name : string
         Name of the variable being transformed into a function.

  Returns
  -------

  made_func : ndarray
               An function representing the data in myst, along axis:
               
               - If myst is a funcion, simply returns myst.
               - If myst is a float, returns a function that returns the value of myst.
               - If myst is an array, returns an function that operates along the same dimension as axis, which returns the value :func:`made_func(axis[i])=myst[i]`, and interpolates values between points in axis.

  """

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
