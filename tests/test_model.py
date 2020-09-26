import sys
import os
import funcsigs
import numpy as np
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc')
from pymoc.model import Model
from pymoc.modules import Column, ModuleWrapper

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
    with pytest.raises(ValueError) as ninfo:
      # model.validate_neighbors_input([atlantic, pacific])
      model.validate_neighbors_input([atlantic, atlantic])
    assert (
        str(ninfo.value) ==
        'Cannot link basins multiple times. Please check your configuration.'
    )

    # def test_add_module(self):

    # def test_new_module(self):

    # def test_run(self):

    # def test_timestep(self):

    # def test_snapshot(self):
