import sys
import pytest
import numpy as np
sys.path.append('/pymoc/src/pymoc/utils')
from make_array import make_array


class TestMakeArray(object):
  def test_make_array(self):
    zlevels = np.asarray(np.linspace(-4000, 0, 80))
    myst = np.arange(0.0, 8.0, 0.1)
    assert all(make_array(myst, zlevels, 'myst') == myst)
    myst = lambda n: 42 + n
    assert all(make_array(myst, zlevels, 'myst') == myst(zlevels))
    myst = 5.0
    assert all(
        make_array(myst, zlevels, 'myst') == 5.0 * np.ones((len(zlevels)))
    )
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      make_array(myst, zlevels, 'myst')
    assert (
        str(
            mystinfo.value
        ) == "('myst', 'needs to be either function, numpy array, or float')"
    )
