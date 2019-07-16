import pytest
import inspect
import sys
import numpy as np
sys.path.append('/pymoc/src/modules')
from psi_SO import Psi_SO

@pytest.fixture(scope="module", params=[
  {
  },
  {
    'z': -2000
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': 1e6
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
    'bs': 0.05,
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
    'bs': 0.05,
    'tau': 0.12
  },
  # {
  #   'z': np.asarray(np.linspace(-4000, 0, 80)),
  #   'y':
  #   'b':
  #   'bs':
  #   'tau':
  #   'f':
  #   'rho':
  #   'L':
  #   'KGM':
  #   'c':
  #   'bvp_with_EK':
  #   'Hsill':
  #   'HEk':
  #   'Htapertop':
  #   'Htaperbot':
  #   'Psi_Ek':
  #   'Psi_GM':
  #   'Psi':
  #   'smax':
  # }
])
def psi_so_config(request):
  return request.param

class TestPsi_SO(object):
  def test_psi_so_init(self, psi_so_config):
    if not 'z' in psi_so_config or not isinstance(psi_so_config['z'], np.ndarray) or not len(psi_so_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Psi_SO(**psi_so_config)
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return
    if not 'y' in psi_so_config or not isinstance(psi_so_config['y'], np.ndarray) or not len(psi_so_config['y']):
      with pytest.raises(TypeError) as yinfo:
        Psi_SO(**psi_so_config)
      assert(str(yinfo.value) == "y needs to be numpy array providing horizontal grid (or boundaries) of ACC")
      return

    for k in ['b', 'bs', 'tau']:
      if not k in psi_so_config:
        with pytest.raises(TypeError) as kinfo:
          Psi_SO(**psi_so_config)
        assert(str(kinfo.value) == "('" + k + "', 'needs to be either function, numpy array, or float')")
        return


    psi_so = Psi_SO(**psi_so_config)
    psi_so_signature = inspect.signature(Psi_SO)
    assert 1