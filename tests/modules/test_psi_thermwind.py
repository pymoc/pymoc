import pytest
import sys
import inspect
import numpy as np
from numpy import testing
sys.path.append('/pymoc/src/modules')
from psi_thermwind import Psi_Thermwind

@pytest.fixture(scope="module", params=[
  {
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'b1': np.linspace(0.03, -0.001, 80),
    'b2': np.linspace(0.02, 0.0, 80),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'b1': np.linspace(0.03, -0.001, 80),
    'b2': np.linspace(0.02, 0.0, 80),
    'sol_init': np.ones((2, 80))
  },
])
def psi_config(request):
  return request.param

class TestPsi_Thermwind(object):
  def test_psi_thermwind_init(self, psi_config):
    if not 'z' in psi_config or not isinstance(psi_config['z'], np.ndarray) or not len(psi_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Psi_Thermwind(**psi_config)
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return

    psi = Psi_Thermwind(**psi_config)
    for k in ['z', 'f', 'sol_init', 'b1', 'b2']:
      assert hasattr(psi, k)

    psi_signature = inspect.signature(Psi_Thermwind)

    for k in ['f']:
      assert getattr(psi, k)  == (psi_config[k] if k in psi_config and psi_config[k] else psi_signature.parameters[k].default)

    if not 'sol_init' in psi_config:
      testing.assert_array_equal(psi.sol_init, np.zeros((2, len(psi.z))) )

    for k in ['b1', 'b2']:
      f = getattr(psi, k)
      ft = psi.make_func(psi_config[k], k, psi_config['z'])
      for z in psi.z:
        assert f(z) == ft(z)
    