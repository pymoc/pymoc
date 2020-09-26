import sys
import os
import funcsigs
import numpy as np
from numpy import testing
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc')
from pymoc.model import Model
from pymoc.modules import Column, Psi_Thermwind, ModuleWrapper

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


    # def test_new_module(self):

    # def test_run(self):

    # def test_timestep(self):

    # def test_snapshot(self):
