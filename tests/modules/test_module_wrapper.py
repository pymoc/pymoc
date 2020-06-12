import sys
import os
import copy
import funcsigs
import numpy as np
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc/modules')
from pymoc.modules import ModuleWrapper, Column, Psi_Thermwind, Psi_SO

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
      'z': np.asarray(np.linspace(-4000, 0, 80)),
      'y': np.asarray(np.linspace(0, 2.0e6, 51)),
      'b': np.linspace(0.03, -0.001, 80),
      'bs': 0.05,
      'tau': 0.12
    }
  ]
)
def psi_so(request):
  return Psi_SO(**request.param)

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
    {
      'module': 'psi_so',
      'name': 'Southern Ocean',
      'left_neighbors': None,
      'right_neighbors': None,
    },
  ]
)
def module_wrapper_config(request, column, psi_thermwind, psi_so):
  param = request.param
  if param['module'] == 'column':
    param['module'] = column
  elif param['module'] == 'psi_thermwind':
    param['module'] = psi_thermwind
  elif param['module'] == 'psi_so':
    param['module'] = psi_so
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
    elif module_wrapper_config['name'] == 'Southern Ocean':
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
      assert module_wrapper.b_type is None
      assert len(module_wrapper.left_neighbors) == 0
      assert len(module_wrapper.right_neighbors) == 1
    elif module_wrapper.name == 'Southern Ocean':
      assert module_wrapper.key == 'southern_ocean'
      assert module_wrapper.do_psi_bz == False
      assert module_wrapper.b_type == 'bs'
      assert len(module_wrapper.left_neighbors) == 0
      assert len(module_wrapper.right_neighbors) == 1

  def test_module_type(self, module_wrapper):
    if module_wrapper.name == 'Atlantic Ocean':
      assert module_wrapper.module_type == 'basin'
    elif module_wrapper.name == 'AMOC':
      assert module_wrapper.module_type == 'coupler'

  def test_timestep_basin(self, mocker, module_wrapper, psi_thermwind, psi_so):
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

      psi_wrapper = ModuleWrapper(name='MOC', module=copy.deepcopy(psi_thermwind))
      psi_so_wrapper = ModuleWrapper(name='SO', module=copy.deepcopy(psi_so))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_right_neighbor(psi_wrapper)
      wrapper.add_left_neighbor(psi_so_wrapper)
      psi_wrapper.update_coupler()
      wA = psi_wrapper.psi[0]*1e6 - psi_so_wrapper.psi[-1]*1e6
      spy = mocker.spy(wrapper.module, 'timestep')
      wrapper.timestep_basin(dt=dt)
      spy.assert_called_once()
      assert spy.call_args[1]['dt'] == dt
      assert all(spy.call_args[1]['wA'] == wA)
    else:
      with pytest.raises(TypeError) as uinfo:
        module_wrapper.timestep_basin(dt=dt)
      assert (
          str(uinfo.value) ==
          "Cannot use timestep_basin on non-basin modules."
      )

  def test_update_coupler(self, mocker, module_wrapper, column):
    if module_wrapper.name == 'Atlantic Ocean':
      with pytest.raises(TypeError) as uinfo:
        module_wrapper.update_coupler()
      assert (
          str(uinfo.value) ==
          "Cannot use update_coupler on non-coupler modules."
      )
    else:
      col_wrapper = ModuleWrapper(name='Atlantic Ocean', module=copy.deepcopy(column))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_left_neighbor(col_wrapper)
      update_spy = mocker.spy(wrapper.module, 'update')
      solve_spy = mocker.spy(wrapper.module, 'solve')
      wrapper.update_coupler()
      update_spy.assert_called_once()
      solve_spy.assert_called_once()
      if wrapper.name == 'AMOC':
        assert all(update_spy.call_args[1]['b1'] == col_wrapper.b)
        assert update_spy.call_args[1]['b2'] == None
        assert all(wrapper.psi[0] == wrapper.module.Psibz()[0])
      elif wrapper.name == 'Southern Ocean':
        assert update_spy.call_args[1]['b'] == None
        assert wrapper.psi == [wrapper.module.Psi]

      col_wrapper = ModuleWrapper(name='Atlantic Ocean', module=copy.deepcopy(column))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_right_neighbor(col_wrapper)
      update_spy = mocker.spy(wrapper.module, 'update')
      solve_spy = mocker.spy(wrapper.module, 'solve')
      wrapper.update_coupler()
      update_spy.assert_called_once()
      solve_spy.assert_called_once()
      if wrapper.name == 'AMOC':
        assert all(update_spy.call_args[1]['b2'] == col_wrapper.b)
        assert update_spy.call_args[1]['b1'] == None
        assert all(wrapper.psi[0] == wrapper.module.Psibz()[0])
      elif wrapper.name == 'Southern Ocean':
        assert all(update_spy.call_args[1]['b'] == col_wrapper.b)
        assert wrapper.psi == [wrapper.module.Psi]

      col_wrapper_left = ModuleWrapper(name='Atlantic Ocean', module=copy.deepcopy(column))
      col_wrapper_right = ModuleWrapper(name='Pacific Ocean', module=copy.deepcopy(column))
      wrapper = copy.deepcopy(module_wrapper)
      wrapper.add_left_neighbor(col_wrapper_left)
      wrapper.add_right_neighbor(col_wrapper_right)
      update_spy = mocker.spy(wrapper.module, 'update')
      solve_spy = mocker.spy(wrapper.module, 'solve')
      wrapper.update_coupler()
      update_spy.assert_called_once()
      solve_spy.assert_called_once()
      if wrapper.name == 'AMOC':
        assert all(update_spy.call_args[1]['b1'] == col_wrapper.b)
        assert all(update_spy.call_args[1]['b2'] == col_wrapper.b)
        assert all(wrapper.psi[0] == wrapper.module.Psibz()[0])
      elif wrapper.name == 'Southern Ocean':
        assert all(update_spy.call_args[1]['b'] == col_wrapper.b)
        assert wrapper.psi == [wrapper.module.Psi] 

  def test_b(self, module_wrapper):
    if module_wrapper.name == "Atlantic Ocean":
      assert all(module_wrapper.b == module_wrapper.module.b)
    elif module_wrapper.name == "Southern Ocean":
      assert module_wrapper.b == module_wrapper.module.bs
    else:
      assert module_wrapper.b is None
    
  def test_neighbors(self, column, psi_thermwind, psi_so):
    amoc_wrapper=ModuleWrapper(name='AMOC', module=psi_thermwind)
    so_wrapper=ModuleWrapper(name='SO', module=psi_so)
    wrapper = ModuleWrapper(
      name='Atlantic Ocean',
      module=column,
      left_neighbors=[amoc_wrapper],
      right_neighbors=[so_wrapper],
    )
    assert wrapper.neighbors == [amoc_wrapper, so_wrapper]
    
  def test_add_left_neighbor(self, mocker, column, psi_thermwind):
    amoc_wrapper=ModuleWrapper(name='AMOC', module=psi_thermwind)
    wrapper = ModuleWrapper(name='Atlantic Ocean', module=column)
    unique_spy = mocker.spy(wrapper, 'validate_neighbor_uniqueness')
    direction_spy = mocker.spy(wrapper, 'validate_coupler_neighbor_direction')
    backlink_spy = mocker.spy(wrapper, 'backlink_neighbor')
    wrapper.add_left_neighbor(amoc_wrapper)
    unique_spy.assert_called_once_with(amoc_wrapper)
    direction_spy.assert_called_once_with('left')
    backlink_spy.assert_called_once_with(amoc_wrapper)
    assert wrapper.left_neighbors == [amoc_wrapper]
    assert amoc_wrapper.right_neighbors == [wrapper]

    amoc_wrapper=ModuleWrapper(name='AMOC', module=psi_thermwind)
    wrapper = ModuleWrapper(name='Atlantic Ocean', module=column)
    unique_spy = mocker.spy(wrapper, 'validate_neighbor_uniqueness')
    direction_spy = mocker.spy(wrapper, 'validate_coupler_neighbor_direction')
    backlink_spy = mocker.spy(wrapper, 'backlink_neighbor')
    wrapper.add_left_neighbor(amoc_wrapper, backlinking=True)
    unique_spy.assert_called_once_with(amoc_wrapper)
    direction_spy.assert_called_once_with('left')
    backlink_spy.assert_not_called()
    assert wrapper.left_neighbors == [amoc_wrapper]
    assert amoc_wrapper.right_neighbors == []

  def test_add_right_neighbor(self, mocker, column, psi_thermwind):
    amoc_wrapper=ModuleWrapper(name='AMOC', module=psi_thermwind)
    wrapper = ModuleWrapper(name='Atlantic Ocean', module=column)
    unique_spy = mocker.spy(wrapper, 'validate_neighbor_uniqueness')
    direction_spy = mocker.spy(wrapper, 'validate_coupler_neighbor_direction')
    backlink_spy = mocker.spy(wrapper, 'backlink_neighbor')
    wrapper.add_right_neighbor(amoc_wrapper)
    unique_spy.assert_called_once_with(amoc_wrapper)
    direction_spy.assert_called_once_with('right')
    backlink_spy.assert_called_once_with(amoc_wrapper)
    assert wrapper.right_neighbors == [amoc_wrapper]
    assert amoc_wrapper.left_neighbors == [wrapper]

    amoc_wrapper=ModuleWrapper(name='AMOC', module=psi_thermwind)
    wrapper = ModuleWrapper(name='Atlantic Ocean', module=column)
    unique_spy = mocker.spy(wrapper, 'validate_neighbor_uniqueness')
    direction_spy = mocker.spy(wrapper, 'validate_coupler_neighbor_direction')
    backlink_spy = mocker.spy(wrapper, 'backlink_neighbor')
    wrapper.add_right_neighbor(amoc_wrapper, backlinking=True)
    unique_spy.assert_called_once_with(amoc_wrapper)
    direction_spy.assert_called_once_with('right')
    backlink_spy.assert_not_called()
    assert wrapper.right_neighbors == [amoc_wrapper]
    assert amoc_wrapper.left_neighbors == []

  def test_add_neighbors(self, mocker, column, psi_thermwind):
    wrapper = ModuleWrapper(name="Atlantic Ocean", module=column)
    left_spy = mocker.spy(wrapper, 'add_left_neighbor')
    right_spy = mocker.spy(wrapper, 'add_right_neighbor')
    wrapper.add_neighbors()
    left_spy.assert_not_called()
    right_spy.assert_not_called()

    wrapper = ModuleWrapper(name="Atlantic Ocean", module=column)
    left_spy = mocker.spy(wrapper, 'add_left_neighbor')
    right_spy = mocker.spy(wrapper, 'add_right_neighbor')
    wrapper.add_neighbors(right_neighbors=[], left_neighbors=[])
    left_spy.assert_not_called()
    right_spy.assert_not_called()

    wrapper = ModuleWrapper(name="Atlantic Ocean", module=column)
    left_spy = mocker.spy(wrapper, 'add_left_neighbor')
    right_spy = mocker.spy(wrapper, 'add_right_neighbor')
    left_neighbors=[
      ModuleWrapper(name="MOC1", module=psi_thermwind),
      ModuleWrapper(name="MOC2", module=psi_thermwind),
      ModuleWrapper(name="MOC3", module=psi_thermwind),
    ]
    right_neighbors=[
      ModuleWrapper(name="ZMOCO1", module=psi_thermwind),
      ModuleWrapper(name="ZMOCO2", module=psi_thermwind),
    ]
    wrapper.add_neighbors(right_neighbors=right_neighbors, left_neighbors=left_neighbors)
    assert left_spy.call_count == 3
    assert right_spy.call_count == 2

  def test_validate_neighbor_uniqueness(self, mocker, module_wrapper, column):
    col_wrapper = ModuleWrapper(name="Pacific Ocean", module=column)
    spy = mocker.spy(module_wrapper, 'validate_neighbor_uniqueness')
    module_wrapper.validate_neighbor_uniqueness(col_wrapper)
    assert spy.spy_exception is None

    module_wrapper.add_left_neighbor(col_wrapper)
    dup_wrapper = ModuleWrapper(name="Indian Ocean", module=column)
    module_wrapper.validate_neighbor_uniqueness(dup_wrapper)
    assert spy.spy_exception is None

    with pytest.raises(KeyError) as valinfo:
      module_wrapper.validate_neighbor_uniqueness(col_wrapper)
    assert (
        str(
            valinfo.value
        ) == "'Cannot add module " + col_wrapper.name + " as a neighbor of " + module_wrapper.name + " because they are already coupled.'"
    )
