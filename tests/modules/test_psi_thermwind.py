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

@pytest.fixture(scope="module")
def psi(request):
  return Psi_Thermwind(**{
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'b1': np.linspace(0.03, -0.01, 80),
    'b2': np.linspace(0.02, 0.0, 80),
    'sol_init': np.ones((2, 80))
  })

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

  def test_make_func(self, psi):
    myst = lambda: 42
    assert psi.make_func(myst, 'myst', psi.z)() == myst()
    myst = np.arange(0.0, 8.0, 0.1)
    for z in psi.z:
      assert psi.make_func(myst, 'myst', psi.z)(z) == np.interp(z, psi.z, myst)
    myst = 6.0
    for z in psi.z:
      assert psi.make_func(myst, 'myst', psi.z)(z) == myst
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      psi.make_func(myst, 'myst', psi.z)
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  def test_make_array(self, psi):
    myst = np.arange(0.0, 8.0, 0.1)
    assert all(psi.make_array(myst, 'myst') == myst)
    myst = lambda n: 42 + n
    assert all(psi.make_array(myst, 'myst') == myst(psi.z))
    myst = 5.0
    assert all(psi.make_array(myst, 'myst') == 5.0 * np.ones((len(psi.z))))
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      psi.make_array(myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  def test_bc(self, psi):
    testing.assert_array_equal(psi.bc([1, 2], [3, 4]), np.array([1, 3]))

  def test_ode(self, psi):
    ode = psi.ode(-1e3, [0, 1])
    b1 = 0.0
    b2 = 0.005
    f = 1.2e-4
    assert ode[0] == 1
    testing.assert_approx_equal(ode[1][0], 1.0/f*(b2 - b1))

  def test_solv(self, psi):
    return

  def test_Psib(self, psi):
    return

  def test_Psibz(self, psi):
    return
  
  def test_update(self, psi):
    return