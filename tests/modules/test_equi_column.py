import pytest
import sys
import inspect
import numpy as np
from numpy import testing
sys.path.append('/pymoc/src/modules')
from equi_column import Equi_Column

@pytest.fixture(scope="module", params=[
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'B_int': 3e3,
    'A': 2.0e14,
    'kappa': 3e-5,
    'H': 500.0
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'B_int': 3e3,
    'A': 2.0e14,
    'kappa': np.asarray(np.linspace(3e-5, 1e-5, 80)),
    'H': 500.0
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'B_int': 3e3,
    'A': 2.0e14,
    'kappa': lambda z: -3e-5 * z / 4e3,
    'H': 500.0
  }
])
def column_config(request):
  return request.param

@pytest.fixture(scope="module")
def column(request):
  return Equi_Column(**{
    'B_int': 3e3,
    'A': 2.0e14,
    'kappa': 3e-5
  })

class TestEqui_Column(object):
  def test_Equi_Column_init(self, column_config):
    column_signature = inspect.signature(Equi_Column)

    column = Equi_Column(**column_config)

    for k in ['f', 'A', 'H', 'H_guess']:
      assert getattr(column, k)  == (column_config[k] if k in column_config and column_config[k] else column_signature.parameters[k].default)
    for k in ['z']:
      testing.assert_array_equal(getattr(column, k), (column_config[k] if k in column_config and not column_config[k] is None else column_signature.parameters[k].default))

    nz = column_config['nz'] if 'nz' in column_config and not column_config['nz'] is None else 100
    for i in range(nz):
      testing.assert_approx_equal(column.zi[i], i/(nz - 1) - 1)

    assert callable(column.kappa)
    assert callable(column.dkappa_dz)
    if 'kappa' in column_config and not column_config['kappa'] is None:
      if callable(column_config['kappa']):
        for z in column_config['z']:
          assert column.kappa(z, 1.0) == column_config['kappa'](z) / (column.f)
      elif isinstance(column_config['kappa'], np.ndarray):
        for i in range(len(column_config['z'])):
          assert column.kappa(column_config['z'][i], 1.0) == column_config['kappa'][i] / (column.f)
      else:
        for z in column_config['z']:
          assert column.kappa(z, 200.0) == column_config['kappa'] / (4e4*column.f)
          assert column.dkappa_dz(z, 200.0) == 0



  def test_alpha(self):
    return 

  def test_bz(self):
    return 

  def test_bc(self):
    return 

  def test_ode(self):
    return 

  def test_solve(self):
    return 