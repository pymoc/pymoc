import numpy as np


def make_array(myst, axis, name):
  # turn mysterious object into array(if needed)
  if isinstance(myst, np.ndarray):
    return myst
  elif callable(myst):
    return myst(axis)
  elif isinstance(myst, float):
    return myst + 0*axis
  else:
    raise TypeError(name, 'needs to be either function, numpy array, or float')
