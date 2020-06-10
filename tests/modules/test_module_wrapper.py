import sys
import os
import funcsigs
import numpy as np
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc/modules')
from pymoc.modules import ModuleWrapper, Column, Psi_Thermwind

@pytest.fixture(
  scope="module",
  params=[
    {
      'Area': 6e13,
      'z': np.asarray(np.linspace(-4000, 0, 80)),
      'kappa': 2e-5,
      'bs': 0.05,
      'bbot': 0.02,
      'bzbot': 0.01,
      'b': 0.03,
      'N2min': 2e-7
    },
  ]
)
def column(request):
  return Column(**request.param)

@pytest.fixture(
  scope="module",
  params=[
    {
      'z': np.asarray(np.linspace(-4000, 0, 80)),
      'b1': np.linspace(0.03, -0.001, 80),
      'b2': np.linspace(0.02, 0.0, 80),
      'sol_init': np.ones((2, 80))
    },
  ]
)
def psi_thermwind(request):
  return Psi_Thermwind(**request.param)

@pytest.fixture(
  scope="module",
  params=[
    {
      'module': 'column',
      'name': 'Atlantic Ocean',
      'left_neighbors': None,
      'right_neighbors': None,
    },
    {
      'module': 'psi_thermwind',
      'name': 'AMOC',
      'left_neighbors': None,
      'right_neighbors': None,
    },
  ]
)
def module_wrapper_config(request, column, psi_thermwind):
  param = request.param
  if param['module'] == 'column':
    param['module'] = column
    param['left_neighbors'] = [ModuleWrapper(
      name='AMOC',
      module=psi_thermwind,
    )]
  elif param['module'] == 'psi_thermwind':
    param['module'] = psi_thermwind
    param['right_neighbors'] = [ModuleWrapper(
      name='Atlantic Ocean',
      module=column,
    )]
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

    assert module_wrapper.psi == [0, 0]

    if module_wrapper.name == 'Atlantic Ocean':
      assert module_wrapper.key == 'atlantic_ocean'
      assert module_wrapper.do_psi_bz == False
      assert module_wrapper.b_type == 'b'
      assert len(module_wrapper.left_neighbors) == 1
      assert len(module_wrapper.right_neighbors) == 0
    elif module_wrapper.name == 'AMOC':
      assert module_wrapper.key == 'amoc'
      assert module_wrapper.do_psi_bz == True
      assert module_wrapper.b_type == 'bs'
      assert len(module_wrapper.left_neighbors) == 0
      assert len(module_wrapper.right_neighbors) == 1

    #   def test_module_type(self):

    #   def test_timestep_basin(self):

    #   def test_update_coupler(self):

    #   def test_b(self):

    #   def test_validate_neighbors(self):

    #   def test_backlink_neighbor_lins(self):
