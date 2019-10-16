import sys
import pytest
import numpy as np
sys.path.append('/pymoc/src/pymoc/utils')
from make_func import make_func


class TestMakeFunc(object):
  def test_make_func(self):
    zlevels = np.asarray(np.linspace(-4000, 0, 80))
    myst = lambda: 42
    assert make_func(myst, zlevels, 'myst')() == myst()
    myst = np.arange(0.0, 8.0, 0.1)
    for z in zlevels:
      assert make_func(myst, zlevels, 'myst')(z) == np.interp(z, zlevels, myst)
    myst = 6.0
    for z in zlevels:
      assert make_func(myst, zlevels, 'myst')(z) == myst
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      make_func(myst, zlevels, 'myst')
    assert (
        str(
            mystinfo.value
        ) == "('myst', 'needs to be either function, numpy array, or float')"
    )
