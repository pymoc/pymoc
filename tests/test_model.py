import sys
import os
import collections
import numpy as np
from numpy import testing
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc')
from pymoc.model import Model
from pymoc.modules import Column, Psi_Thermwind, ModuleWrapper

COLUMN_PARAMS = {
      'Area': 6e13,
      'z': np.asarray(np.linspace(-4000, 0, 80)),
      'kappa': 2e-5,
      'bs': 0.05,
      'bbot': 0.02,
      'bzbot': 0.01,
      'b': np.linspace(0.03, -0.001, 80),
      'N2min': 2e-7
    }
PSI_PARAMS = {
      'z': np.asarray(np.linspace(-4000, 0, 80)),
      'b1': np.linspace(0.03, -0.001, 80),
      'b2': np.linspace(0.02, 0.0, 80),
      'sol_init': np.ones((2, 80))
    }

@pytest.fixture(
  scope="module",
  params=[COLUMN_PARAMS,]
)
def column(request):
  return Column(**request.param)

@pytest.fixture(
  scope="module",
  params=[PSI_PARAMS,]
)
def psi_thermwind(request):
  return Psi_Thermwind(**request.param)

class TestModel(object):
  def test_model_init(self):
    model = Model()
    assert model.basins == []
    assert model.couplers == []
    assert model._modules == {}
  # def __get_item__(self):

  def test_keys(self, column):
    model = Model()
    model.add_module(column, 'Atlantic')
    assert list(model.keys) == ['atlantic']

  def test_get_module(self, column):
    model = Model()
    model.add_module(column, 'Atlantic')
    assert model.get_module('atlantic').module == column
    with pytest.raises(KeyError) as pinfo:
      test_column = model.get_module('pacific')
      assert test_column == None
    assert str(pinfo.value) == "'pacific'"


  def test_validate_neighbors_input(self, column):
    model = Model()
    atlantic = ModuleWrapper(module=column, name='Atlantic')
    pacific = ModuleWrapper(module=column, name='Pacific')
    # Validation passes if no exception raised
    model.validate_neighbors_input([atlantic, pacific])
    with pytest.raises(ValueError) as ninfo:
      model.validate_neighbors_input([atlantic, atlantic])
    assert (
        str(ninfo.value) ==
        'Cannot link basins multiple times. Please check your configuration.'
    )

  def test_validate_new_module_key(self, column):
    model = Model()
    atlantic = ModuleWrapper(module=column, name='Atlantic')
    pacific = ModuleWrapper(module=column, name='Pacific')
    model.validate_new_module_key(atlantic)
    model.basins.append(atlantic)
    model._modules[atlantic.key] = atlantic
    setattr(model, atlantic.key, atlantic.module)
    model.validate_new_module_key(pacific)
    with pytest.raises(NameError) as ninfo:
      model.validate_new_module_key(atlantic)
    assert (
        str(ninfo.value) ==
        'Cannot use module name Atlantic because it would overwrite an existing key or model property.'
    )
    model.validate_new_module_key(pacific)
    keys = ModuleWrapper(module=column, name='Keys')
    with pytest.raises(NameError) as ninfo:
      model.validate_new_module_key(keys)
    assert (
        str(ninfo.value) ==
        'Cannot use module name Keys because it would overwrite an existing key or model property.'
    )


  def test_add_module(self, mocker, column, psi_thermwind):
    atlantic = ModuleWrapper(module=column, name='Atlantic')
    pacific = ModuleWrapper(module=column, name='Pacific')
    amoc = ModuleWrapper(module=psi_thermwind, name='AMOC')
    zoc = ModuleWrapper(module=psi_thermwind, name='ZOC')

    model = Model()
    neighbor_input_spy = mocker.spy(model, 'validate_neighbors_input')
    new_module_key_spy = mocker.spy(model, 'validate_new_module_key')
    model.add_module(column, 'Atlantic', left_neighbors=[amoc], right_neighbors=[zoc]) 
    neighbor_input_spy.assert_called_once()
    assert neighbor_input_spy.call_args[0] == ([amoc, zoc],)
    new_module_key_spy.assert_called_once()
    assert new_module_key_spy.call_args[0][0].module == atlantic.module
    assert new_module_key_spy.call_args[0][0].name == atlantic.name
    assert new_module_key_spy.call_args[0][0].key == atlantic.key
    assert model['atlantic'].module == atlantic.module
    assert model.atlantic == atlantic.module
    assert atlantic.module in [b.module for b in model.basins]

    model = Model()
    neighbor_input_spy = mocker.spy(model, 'validate_neighbors_input')
    new_module_key_spy = mocker.spy(model, 'validate_new_module_key')
    model.add_module(psi_thermwind, 'AMOC', left_neighbors=[atlantic], right_neighbors=[pacific]) 
    neighbor_input_spy.assert_called_once()
    assert neighbor_input_spy.call_args[0] == ([atlantic, pacific],)
    new_module_key_spy.assert_called_once()
    assert new_module_key_spy.call_args[0][0].module == amoc.module
    assert new_module_key_spy.call_args[0][0].name == amoc.name
    assert new_module_key_spy.call_args[0][0].key == amoc.key
    assert model['amoc'].module == amoc.module
    assert model.amoc == amoc.module
    assert amoc.module in [c.module for c in model.couplers]


  def test_new_module(self, mocker, column, psi_thermwind):
    atlantic = ModuleWrapper(module=column, name='Atlantic')
    pacific = ModuleWrapper(module=column, name='Pacific')
    amoc = ModuleWrapper(module=psi_thermwind, name='AMOC')
    zoc = ModuleWrapper(module=psi_thermwind, name='ZOC')

    model = Model()
    add_module_spy = mocker.spy(model, 'add_module')
    model.new_module(Column, COLUMN_PARAMS, 'Atlantic', left_neighbors=[amoc], right_neighbors=[zoc])
    add_module_spy.assert_called_once() 
    assert isinstance(add_module_spy.call_args[0][0], Column)
    for key in add_module_spy.call_args[0][0].__dict__.keys():
      if callable(getattr(add_module_spy.call_args[0][0], key)):
        test = getattr(add_module_spy.call_args[0][0], key)(1) == getattr(atlantic.module, key)(1)
      else:
        test = getattr(add_module_spy.call_args[0][0], key) == getattr(atlantic.module, key)
      if isinstance(test, collections.Iterable):
        assert all(test)
      else:
        assert test
    assert add_module_spy.call_args[1]['left_neighbors'] ==  [amoc]
    assert add_module_spy.call_args[1]['right_neighbors'] == [zoc]

    model = Model()
    add_module_spy = mocker.spy(model, 'add_module')
    model.new_module(Psi_Thermwind, PSI_PARAMS, 'AMOC', left_neighbors=[atlantic], right_neighbors=[pacific])
    add_module_spy.assert_called_once() 
    assert isinstance(add_module_spy.call_args[0][0], Psi_Thermwind)
    for key in add_module_spy.call_args[0][0].__dict__.keys():
      if callable(getattr(add_module_spy.call_args[0][0], key)):
        test = getattr(add_module_spy.call_args[0][0], key)(1) == getattr(amoc.module, key)(1)
      else:
        test = getattr(add_module_spy.call_args[0][0], key) == getattr(amoc.module, key)
      while isinstance(test, collections.Iterable) and len(np.shape(test)) > 1:
        test = [t for subtest in test for t in subtest]
      if isinstance(test, collections.Iterable):
        assert all(test)
      else:
        assert test
    assert add_module_spy.call_args[1]['left_neighbors'] ==  [atlantic]
    assert add_module_spy.call_args[1]['right_neighbors'] == [pacific]

  def test_run(self, mocker, column, psi_thermwind):
    model = Model()
    steps = 100
    dt = 1
    coupler_dt = 10

    model.new_module(Psi_Thermwind, PSI_PARAMS, 'AMOC')
    model.new_module(Psi_Thermwind, PSI_PARAMS, 'ZOC')
    model.add_module(column, 'Atlantic', left_neighbors=[model.get_module('amoc')], right_neighbors=[model.get_module('zoc')]) 

    timestep_spy = mocker.spy(model, 'timestep')
    amoc_update_spy = mocker.spy(model.amoc, 'update')
    zoc_update_spy = mocker.spy(model.zoc, 'update')

    model.run(steps, dt, coupler_dt=coupler_dt)

    assert timestep_spy.call_count == steps
    assert amoc_update_spy.call_count == steps * dt / coupler_dt + 1
    assert zoc_update_spy.call_count == steps * dt / coupler_dt + 1

  def test_run_with_snapshots(self, mocker, column, psi_thermwind):
    model = Model()
    steps = 100
    dt = 1
    coupler_dt = 10
    snapshot_start = 5
    snapshot_interval = 5

    model.new_module(Psi_Thermwind, PSI_PARAMS, 'AMOC')
    model.new_module(Psi_Thermwind, PSI_PARAMS, 'ZOC')
    model.add_module(column, 'Atlantic', left_neighbors=[model.get_module('amoc')], right_neighbors=[model.get_module('zoc')]) 

    timestep_spy = mocker.spy(model, 'timestep')
    amoc_update_spy = mocker.spy(model.amoc, 'update')
    zoc_update_spy = mocker.spy(model.zoc, 'update')

    snapshots = [s for s in model.run_with_snapshots(steps, dt, coupler_dt=coupler_dt, snapshot_start=snapshot_start, snapshot_interval=snapshot_interval)]
    snapshot_test = list(np.arange(snapshot_start, steps, snapshot_interval))
    if not snapshot_test[-1] == steps - 1:
      snapshot_test.append(steps - 1)

    assert snapshots == snapshot_test
    assert timestep_spy.call_count == steps
    assert amoc_update_spy.call_count == steps * dt / coupler_dt + 1
    assert zoc_update_spy.call_count == steps * dt / coupler_dt + 1
    
  def test_timestep(self, mocker, column):
    model = Model()
    steps = 100
    dt = 1
    coupler_dt = 10

    model.new_module(Psi_Thermwind, PSI_PARAMS, 'AMOC')
    model.new_module(Psi_Thermwind, PSI_PARAMS, 'ZOC')
    model.add_module(column, 'Atlantic', left_neighbors=[model.get_module('amoc')], right_neighbors=[model.get_module('zoc')]) 

    amoc_update_spy = mocker.spy(model.get_module('amoc'), 'update_coupler')
    zoc_update_spy = mocker.spy(model.get_module('zoc'), 'update_coupler')
    atlantic_timestep_spy = mocker.spy(model.get_module('atlantic'), 'timestep_basin')

    model.timestep(5, dt, coupler_dt=coupler_dt)

    assert atlantic_timestep_spy.call_count == 1
    assert amoc_update_spy.call_count == 0
    assert zoc_update_spy.call_count == 0

    model.timestep(20, dt, coupler_dt=coupler_dt)

    assert atlantic_timestep_spy.call_count == 2
    assert amoc_update_spy.call_count == 1
    assert zoc_update_spy.call_count == 1

  def test_snapshot(self):
    model = Model()

    assert not model.snapshot(10, None, None)
    assert not model.snapshot(10, 1, None)
    assert not model.snapshot(1, 10, 10)
    assert model.snapshot(10, 10, 10)
    assert model.snapshot(50, 10, 10)
