import sys
import os
import funcsigs
import numpy as np
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc/modules')
from pymoc.modules import ModuleWrapper, Column

@pytest.fixture(scope="module")
def column(request):
  return Column(
      **{
          'Area': 6e13,
          'z': np.asarray(np.linspace(-4000, 0, 80)),
          'kappa': 2e-5,
          'bs': 0.05,
          'bbot': 0.02,
          'bzbot': 0.01,
          'b': 0.03,
          'N2min': 2e-7
      }
  )
@pytest.fixture(
  scope="module",
  params=[
    {
      'module': 'column',
      'name': 'Atlantic Ocean',
      'left_neighbors': None,
      'right_neighbors': None,
    },
  ]
)
def module_wrapper_config(request, column):
  param = request.param
  if param['module'] == 'column':
    param['module'] = column
  return param

class TestModuleWrapper(object):
  def test_module_wrapper_init(self, module_wrapper_config):
    print(module_wrapper_config)
    module_wrapper = ModuleWrapper(**module_wrapper_config)
    for k in ['module', 'name', 'key', 'do_psi_bz', 'b_type', 'psi', 'left_neighbors', 'right_neighbors']:
      assert hasattr(module_wrapper, k)

    module_wrapper_signature = funcsigs.signature(ModuleWrapper)
    # The constructor initializes all scalar properties
    # Uses explicit property if present, or the default
    for k in ['module', 'name']:
      assert getattr(module_wrapper, k) == (
          module_wrapper_config[k] if k in module_wrapper_config and module_wrapper_config[k] else
          module_wrapper_signature.parameters[k].default
      )
    #   def test_module_type(self):

    #   def test_timestep_basin(self):

    #   def test_update_coupler(self):

    #   def test_b(self):

    #   def test_validate_neighbors(self):

    #   def test_backlink_neighbor_lins(self):
