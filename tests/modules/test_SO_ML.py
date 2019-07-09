import pytest
import inspect
import numpy as np
from scipy import integrate
from pymoc.SO_ML import SO_ML

@pytest.fixture(scope="module")
def so_ml(request):
  return SO_ML(**{ 'y': np.asarray(np.linspace(-70.0, -55.0, 60)) }) 

class TestSO_ML(object):
  def test_make_array(self, so_ml):
    myst = np.arange(0.0, 8.0, 0.1)
    assert all(so_ml.make_array(myst, 'myst') == myst)
    myst = lambda n: 42 + n
    assert all(so_ml.make_array(myst, 'myst') == myst(so_ml.y))
    myst = 5.0
    assert all(so_ml.make_array(myst, 'myst') == 5.0 * np.ones((len(so_ml.y))))
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      so_ml.make_array(myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  def test_solve_equi(self, so_ml):
    with pytest.raises(TypeError) as info:
      so_ml.solve_equi()
    assert(str(info.value) == "This functionality is not yet implemented")

  def test_timestep(self):
    