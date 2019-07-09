import pytest
import inspect
import numpy as np
from scipy import integrate
from pymoc.SO_ML import SO_ML

@pytest.fixture(scope="module", params=[
  {
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'Ks': 100,
    'h': 50,
    'L': 4e6,
    'surflux': 5.9e3,
    'rest_mask': 0.0,
    'b_rest': 0.0,
    'v_pist': 2.0/86400.0,
    'bs': 0.02
  },
  {
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
  },
  {
    'y': None,
  },
  {
    'y': 1e6
  }
])
def so_ml_config(request):
  return request.param

@pytest.fixture(scope="module")
def so_ml(request):
  return SO_ML(y=np.asarray(np.linspace(0, 2.0e6, 51))) 

class TestSO_ML(object):
  def test_so_ml_init(self, so_ml_config):
    if not isinstance(so_ml_config['y'], np.ndarray) or not len(so_ml_config['y']):
      with pytest.raises(TypeError) as yinfo:
        SO_ML(**so_ml_config)
      assert(str(yinfo.value) == "y needs to be numpy array providing (regular) grid")
      return

    so_ml = SO_ML(**so_ml_config)
    for k in ['y', 'Ks', 'h', 'L', 'surflux', 'rest_mask', 'b_rest', 'v_pist', 'Psi_s', 'bs']:
      assert hasattr(so_ml, k)

    so_ml_signature = inspect.signature(SO_ML)

    for k in ['Ks', 'h', 'L', 'v_pist', 'Psi_s']:
      assert getattr(so_ml, k)  == (so_ml_config[k] if k in so_ml_config and so_ml_config[k] else so_ml_signature.parameters[k].default)

    for k in ['surflux', 'rest_mask', 'b_rest', 'bs']:
      assert all(getattr(so_ml, k)  == so_ml.make_array((so_ml_config[k] if k in so_ml_config and so_ml_config[k] else so_ml_signature.parameters[k].default), k))

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
