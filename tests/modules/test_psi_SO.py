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
  #   'bvp_with_EskK':
  #   'Hsill':
  #   'HEk':
  #   'Htapertop':
  #   'Htaperbot':
  #   'smax':
  # }
])
def psi_so_config(request):
  return request.param

@pytest.fixture(scope="module")
def psi_so(request):
  return Psi_SO(**{
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
    'bs': 0.05,
    'tau': 0.12
  })

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

    # Implicit exceptions thrown if these params are missing, should change these to explicit checks
    for k in ['b', 'bs', 'tau']:
      if not k in psi_so_config:
        with pytest.raises(TypeError) as kinfo:
          Psi_SO(**psi_so_config)
        assert(str(kinfo.value) == "('" + k + "', 'needs to be either function, numpy array, or float')")
        return

    psi_so = Psi_SO(**psi_so_config)
    for k in ['z', 'y', 'b', 'bs', 'tau', 'f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek', 'Hsill', 'HEk', 'Htapertop', 'Htaperbot', 'smax']:
      assert hasattr(psi_so, k)

    psi_so_signature = inspect.signature(Psi_SO)
    
    for k in ['f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek', 'Hsill', 'HEk', 'Htapertop', 'Htaperbot', 'smax']:
      assert getattr(psi_so, k)  == (psi_so_config[k] if k in psi_so_config and psi_so_config[k] else psi_so_signature.parameters[k].default)

    for k in ['b']:
      f = getattr(psi_so, k)
      ft = psi_so.make_func(psi_so_config['z'], psi_so_config[k], k)
      for z in psi_so.z:
        assert f(z) == ft(z)

    for k in ['bs', 'tau']:
      f = getattr(psi_so, k)
      ft = psi_so.make_func(psi_so_config['y'], psi_so_config[k], k)
      for y in psi_so.y:
        assert f(y) == ft(y)

  def test_make_func(self, psi_so):
    myst = lambda: 42
    assert psi_so.make_func(psi_so.z, myst, 'myst')() == myst()
    myst = np.arange(0.0, 8.0, 0.1)
    for z in psi_so.z:
      assert psi_so.make_func(psi_so.z, myst, 'myst')(z) == np.interp(z, psi_so.z, myst)
    myst = 6.0
    for y in psi_so.y:
      assert psi_so.make_func(psi_so.y, myst, 'myst')(y) == myst
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      psi_so.make_func(psi_so.z, myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")