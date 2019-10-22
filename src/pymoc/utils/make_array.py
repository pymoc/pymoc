import numpy as np


def make_array(myst, axis, name):
  r"""
  Make an argument of unknown type into an array if needed, along the specified array.

  Parameters
  ----------

  myst : float, function, or ndarray
         The argument to be transformed. Either a single numerical value, a function in a single dimension, or an array.
  axis : ndarray
         Grid points in the direction along which the array is to be created.
  name : string
         Name of the variable being transformed into an array.

  Returns
  -------

  made_array : ndarray
               An array representing the data in myst, along axis:
               
               - If myst is an ndarray, simply returns myst.
               - If myst is a float, returns an array of length :func:`len(axis)` where each row contains the value of myst.
               - If myst is a function, returns an array of length :func:`len(axis)` where each row contains the value :func:`made_array[i]=myst(axis[i])`.

  """

  if isinstance(myst, np.ndarray):
    return myst
  elif callable(myst):
    return myst(axis)
  elif isinstance(myst, float):
    return myst + 0*axis
  else:
    raise TypeError(name, 'needs to be either function, numpy array, or float')
