import sys
import os
import inspect
import pprint as pp
import numpy as np
from scipy import integrate
import pytest
from collections.abc import Iterable
from pymoc.column import Column

@pytest.fixture(scope="module", params=[
  { 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 },
  { 'Area': None, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5 },
  { 'Area': 6e13, 'z': None, 'kappa': 2e-5 },
  { 'Area': 6e13, 'z': 50, 'kappa': 2e-5 },
  { 'Area': 6e13, 'z': np.array([]), 'kappa': 2e-5 },
  { 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': None },
  { 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'bzbot': 0.01, 'b': 0.03, 'N2min': 2e-7 },
  { 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'b': np.arange(0.0, 8.0, 0.1) },
])
def column_config(request):
  return request.param

@pytest.fixture(scope="module")
def column(request):
  return Column(**{ 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'bzbot': 0.01, 'b': 0.03, 'N2min': 2e-7 }) 

class TestColumn(object):
  def test_column_init(self, column_config):
    if not column_config['Area']:
      with pytest.raises(TypeError) as areainfo:
        Column(**column_config)
      assert(str(areainfo.value) == "('Area', 'needs to be either function, numpy array, or float')")
      return

    if not column_config['kappa']:
      with pytest.raises(TypeError) as kappainfo:
        Column(**column_config)
      assert(str(kappainfo.value) == "('kappa', 'needs to be either function, numpy array, or float')")
      return

    if not isinstance(column_config['z'], np.ndarray) or not len(column_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Column(**column_config)
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return

    column = Column(**column_config)

    # The constructor assigns all expected properties
    assert hasattr(column, 'z')
    assert hasattr(column, 'kappa')
    assert hasattr(column, 'Area')
    assert hasattr(column, 'b')
    assert hasattr(column, 'bs')
    assert hasattr(column, 'bbot')
    assert hasattr(column, 'bzbot')
    assert hasattr(column, 'N2min')

    column_signature = inspect.signature(Column)

    # The constructor initializes all scalar properties
    # Uses explicit property if present, or the default
    for k in ['bs', 'bbot', 'bzbot', 'N2min']:
      assert getattr(column, k)  == (column_config[k] if k in column_config and column_config[k] else column_signature.parameters[k].default)

    # The constructor initializes all vector properties
    # Uses explicit property if present, or the default
    assert all(column.z == column_config['z'])
    assert all(column.b  == column.make_array(column_config['b'] if 'b' in column_config and not column_config['b'] is None else column_signature.parameters['b'].default, 'b'))

    # The constructor initializes all z-dependent callable properties
    # Uses explicit property if present, or the default
    for k in ['Area', 'kappa']:
      f = getattr(column, k)
      ft = column.make_func(column_config[k], k)
      for z in column.z:
        assert f(z) == ft(z)

  def test_make_func(self, column):
    myst = lambda: 42
    assert column.make_func(myst, 'myst')() == myst()
    myst = np.arange(0.0, 8.0, 0.1)
    for z in column.z:
      assert column.make_func(myst, 'myst')(z) == np.interp(z, column.z, myst)
    myst = 6.0
    for z in column.z:
      assert column.make_func(myst, 'myst')(z) == myst
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      column.make_func(myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")
  
  def test_make_array(self, column):
    myst = np.arange(0.0, 8.0, 0.1)
    assert all(column.make_array(myst, 'myst') == myst)
    myst = lambda n: 42 + n
    assert all(column.make_array(myst, 'myst') == myst(column.z))
    myst = 5.0
    assert all(column.make_array(myst, 'myst') == 5.0 * np.ones((len(column.z))))
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      column.make_array(myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  def test_Akappa(self, column):
    # Constant area and kappa
    for z in column.z:
      assert column.Akappa(z) == column.Area(0) * column.kappa(0)
    # variable area and kappa
    z = column.z
    Area = lambda l: float(column.Area(0)) / float(l + 1)
    kappa = lambda l: float(column.kappa(0)) * float(l + 1)
    col = Column(z=z, Area=Area, kappa=kappa)
    for z in column.z:
      assert col.Akappa(z) == Area(z) * kappa(z)

  def test_dAkappa_dz(self, column):
    # Constant area and kappa
    for z in column.z:
      assert column.dAkappa_dz(z) == np.gradient(column.Area(0) * column.kappa(0), z)
    # variable area and kappa
    z = column.z
    Area = lambda l: float(column.Area(0)) / float(l + 1)
    kappa = lambda l: float(column.kappa(0)) * float(l + 1)
    col = Column(z=z, Area=Area, kappa=kappa)
    for z in column.z:
      assert col.Akappa(z) == Area(z) * kappa(z)

  def test_bc(self):
    column = Column(**{ 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'bzbot': 0.01, 'b': 0.03, 'N2min': 2e-7 })
    ya = [column.b[0], column.bz[0]]
    yb = [column.b[-1], column.bz[-1]]
    bc = column.bc(ya, yb) 
    assert np.round(bc[0], decimals=2)  == -0.01
    assert np.round(bc[1], decimals=2)  == -0.02

    column = Column(**{ 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'b': 0.03, 'N2min': 2e-7 })
    ya = [column.b[0], column.bz[0]]
    yb = [column.b[-1], column.bz[-1]]
    bc = column.bc(ya, yb) 
    assert np.round(bc[0], decimals=2)  == 0.01
    assert np.round(bc[1], decimals=2)  == -0.02
      
  def test_ode(self, column):
    column.wA = np.sin
    assert (column.ode(column.z, [column.b, column.bz]) == np.vstack((column.bz, (np.sin(column.z) - column.dAkappa_dz(column.z)) / column.Akappa(column.z) * column.bz))).all()

  def test_solve_equi(self):
    column = Column(**{ 'Area': 6e13, 'z': np.asarray(np.linspace(-4000, 0, 80)), 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'bzbot': 0.01, 'b': 0.03, 'N2min': 2e-7 })
    column.wA = np.sin
    sol = integrate.solve_bvp(column.ode, column.bc, column.z, [column.b, column.bz])
    sol_values = sol.sol(column.z)
    column.solve_equi(column.wA)
    assert all(np.around(column.b, decimals=2) == np.around(np.asarray(np.linspace(-39.95, 0.05, 80)), decimals=2))
    assert all(np.around(column.bz, decimals=2) == np.around(sol_values[1,:], decimals=2))

  def test_vertadvdiff(self):
    z = np.asarray(np.linspace(-4000, 0, 80))
    # For constant kappa, bottom stratification provided, negative weff
    column = Column(**{ 'Area': 6e13, 'z': z, 'kappa': 2e-5, 'bs': 0.05, 'bbot': 0.02, 'bzbot': 0.01, 'b': 0.03, 'N2min': 2e-7 })
    wA = np.sin(z)
    b = column.b.copy()
    dt = 0.1
    dz = 50.63291139 * np.ones(len(b) - 1)
    b[-1] = 0.05
    b[0] = b[1] - 0.01 * 50.63291139
    bz = (b[1:] - b[:-1]) / dz
    bz_up = bz[1:]
    bz_down = bz[:-1]
    bzz = (bz_up - bz_down) / (0.5*(dz[1:] + dz[:-1]))
    weff = wA - column.dAkappa_dz(column.z)
    bz = bz_down
    bz[weff[1:-1] < 0] = bz_up[weff[1:-1] < 0]      
    db_dt = (
      -weff[1:-1] * bz / 6e13
      + 2e-5 * bzz
    )
    b[1:-1] = b[1:-1]+dt*db_dt
    column.vertadvdiff(wA, dt)
    assert all(np.around(column.b, decimals=4) == np.around(b, decimals=4))
