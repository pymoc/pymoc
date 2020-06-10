import sys
import os
import copy
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
      # 'b': 0.03,
      'b': np.linspace(0.03, -0.001, 80),
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
  elif param['module'] == 'psi_thermwind':
    param['module'] = psi_thermwind
  return param

@pytest.fixture(scope='module')
def module_wrapper(request, module_wrapper_config):
  return ModuleWrapper(**module_wrapper_config)

class TestModuleWrapper(object):
  def test_module_wrapper_init(self, module_wrapper_config, column, psi_thermwind):
    module_wrapper_config = copy.deepcopy(module_wrapper_config)
    if module_wrapper_config['name'] == 'Atlantic Ocean':
      module_wrapper_config['left_neighbors'] = [ModuleWrapper(name='MOC', module=copy.deepcopy(psi_thermwind))]
    elif module_wrapper_config['name'] == 'AMOC':
      module_wrapper_config['right_neighbors'] = [ModuleWrapper(name='Pacific Ocean', module=copy.deepcopy(column))]

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

  def test_module_type(self, module_wrapper):
    if module_wrapper.name == 'Atlantic Ocean':
      assert module_wrapper.module_type == 'basin'
    elif module_wrapper.name == 'AMOC':
      assert module_wrapper.module_type == 'coupler'

  def test_timestep_basin(self, mocker, module_wrapper, psi_thermwind):
    dt = 0.1
    if module_wrapper.name == 'Atlantic Ocean':
      psi_wrapper = ModuleWrapper(name='MOC', module=copy.deepcopy(psi_thermwind))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_left_neighbor(psi_wrapper)
      psi_wrapper.update_coupler()
      wA = -psi_wrapper.psi[-1]*1e6
      spy = mocker.spy(wrapper.module, 'timestep')
      wrapper.timestep_basin(dt=dt)

      spy.assert_called_once()
      assert spy.call_args[1]['dt'] == dt
      assert all(spy.call_args[1]['wA'] == wA)

      psi_wrapper = ModuleWrapper(name='MOC', module=copy.deepcopy(psi_thermwind))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_right_neighbor(psi_wrapper)
      psi_wrapper.update_coupler()
      wA = psi_wrapper.psi[0]*1e6
      spy = mocker.spy(wrapper.module, 'timestep')
      wrapper.timestep_basin(dt=dt)

      spy.assert_called_once()
      assert spy.call_args[1]['dt'] == dt
      assert all(spy.call_args[1]['wA'] == wA)

    #   def test_update_coupler(self):

    #   def test_b(self):

    #   def test_validate_neighbors(self):

    #   def test_backlink_neighbor_lins(self):
