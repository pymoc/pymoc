import numpy as np


def check_numpy_version():
  r"""
  Check whether your system is using a numpy version of at least 1.13,
  required for some modules for gradient calculations.

  Returns
  -------

  is_sufficient : logical
                  True if numpy version >= 1.13, False otherwise.

  """

  # check numpy version (version >= 1.13 needed to automatically compute db/dz)
  v = [int(i) for i in np.version.version.split('.')]
  if v[0] <= 1 and v[1] < 13:
    return False
  return True
