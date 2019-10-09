import numpy as np


def check_numpy_version():
  # check numpy version (version >= 1.13 needed to automatically compute db/dz)
  v = [int(i) for i in np.version.version.split('.')]
  if v[0] <= 1 and v[1] < 13:
    return False
  return True
