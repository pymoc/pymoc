import numpy as np


def gridit(x1, x2, f):
  r"""
  Generate a gridded dataset, based on a function in x1 and x2.

  Parameters
  ----------

  x1 : ndarray
       Grid points in the first dimension along which values are to be calculated.
  x2 : ndarray
       Grid points in the second dimension along which values are to be calculated.
  f : function
      A function in dimensions x1 and x2, that returns a single value for each pair of coordinate values.

  Returns
  -------
  gridded: ndarray
           A 2D array of shape :func:`(len(x1), len(x2))`, where each point :func:`gridded[i, j] = f(x1[i], x2[j])`.

  """
  n1 = len(x1)
  n2 = len(x2)
  array = np.zeros((n1, n2))
  for i in range(0, n1):
    for j in range(0, n2):
      array[i, j] = f(x1[i], x2[j])
  return array
