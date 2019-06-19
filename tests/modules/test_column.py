import sys
import os
import inspect
import pprint as pp
import numpy as np
import pytest
from collections.abc import Iterable
from pymoc.column import Column

@pytest.fixture(scope="module", params=[
  { 'A': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 },
  {'A': None, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 },
  { 'A': 6e13, 'z': None, 'kappa': 2e-5 },
  { 'A': 6e13, 'z': 50, 'kappa': 2e-5 },
  { 'A': 6e13, 'z': np.array([]), 'kappa': 2e-5 },
  { 'A': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': None },
]
)
def column_config(request):
  return request.param
  # return Column(z=request.param['z'], kappa=Grequest.param['kappa'], Area=request.param['A'])

class TestColumn(object):
  def test_column_init(self, column_config):
    init_col = lambda: Column(z=column_config['z'], kappa=column_config['kappa'], Area=column_config['A'])
    if not column_config['A']:
      with pytest.raises(TypeError) as areainfo:
        init_col()
      assert(str(areainfo.value) == "('Area', 'needs to be either function, numpy array, or float')")
      return

    if not column_config['kappa']:
      with pytest.raises(TypeError) as kappainfo:
        init_col()
      assert(str(kappainfo.value) == "('kappa', 'needs to be either function, numpy array, or float')")
      return

    if not isinstance(column_config['z'], np.ndarray) or not len(column_config['z']):
      with pytest.raises(TypeError) as zinfo:
        init_col()
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return

    column = init_col()

    assert hasattr(column, 'z')
    assert hasattr(column, 'kappa')
    assert hasattr(column, 'Area')
    assert hasattr(column, 'b')
    assert hasattr(column, 'bs')
    assert hasattr(column, 'bbot')
    assert hasattr(column, 'bzbot')
    assert hasattr(column, 'N2min')

    column_signature = inspect.signature(Column)
    for k in ['bs', 'bbot', 'bzbot', 'N2min']:
      assert getattr(column, k)  == (column_config[k] if k in column_config and column_config[k] else column_signature.parameters[k].default)
  