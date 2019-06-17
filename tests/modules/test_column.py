import sys
import os
import numpy as np
import pytest
src_path = os.path.realpath(os.path.dirname(__file__))
src_path = '/'.join(src_path.split('/')[0:-2]) + '/src/modules'
sys.path.append(src_path)
from column import Column

@pytest.fixture(scope="module", params=[
  { 'A': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 },
  {'A': None, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 }]
)
def column_config(request):
  return request.param
  # return Column(z=request.param['z'], kappa=Grequest.param['kappa'], Area=request.param['A'])

class TestColumn(object):
  def test_column_init(self, column_config):
    if not column_config['A']:
      with pytest.raises(TypeError) as areainfo:
        print('bar')
        column = Column(
          z=column_config['z'],
          kappa=column_config['kappa'],
          Area=column_config['A']
        )
      assert(str(areainfo.value) == "('Area', 'needs to be either function, numpy array, or float')")
      return

    column = Column(
      z=column_config['z'],
      kappa=column_config['kappa'],
      Area=column_config['A']
    )

    assert hasattr(column, 'z')
    assert hasattr(column, 'kappa')
    assert hasattr(column, 'Area')
    assert hasattr(column, 'b')
    assert hasattr(column, 'bs')
    assert hasattr(column, 'bbot')
    assert hasattr(column, 'bzbot')
    assert hasattr(column, 'N2min')
  