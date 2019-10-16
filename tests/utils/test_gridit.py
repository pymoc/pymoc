import sys
import pytest
import numpy as np
sys.path.append('/pymoc/src/pymoc/utils')
from gridit import gridit


class TestGridit(object):
  def test_gridit(self):
    f = lambda x, y: x ^ 2 + y
    nx = 3
    ny = 2
    test = np.asarray([[2, 3], [3, 2], [0, 1]])
    sol = gridit(range(nx), range(ny), f)

    for i in range(nx):
      for j in range(ny):
        assert test[i, j] == sol[i, j]
