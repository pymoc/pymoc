import sys
import os
import funcsigs
import numpy as np
from scipy import integrate
import pytest
from pymoc.utils import make_func, make_array
sys.path.append('/pymoc/src/pymoc/modules')
from column import Column


@pytest.fixture(
    scope="module",
    params=[
        {
            'Area': 6e13,
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'kappa': 2e-5
        },
        {
            'Area': None,
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'kappa': 2e-5
        },
        {
            'Area': 6e13,
            'z': None,
            'kappa': 2e-5
        },
        {
            'Area': 6e13,
            'z': 50,
            'kappa': 2e-5
        },
        {
            'Area': 6e13,
            'z': np.array([]),
            'kappa': 2e-5
        },
        {
            'Area': 6e13,
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'kappa': None
        },
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
        {
            'Area': 6e13,
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'kappa': 2e-5,
            'b': np.linspace(0.03, -0.001, 80)
        },
    ]
)
def column_config(request):
  return request.param


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


class TestColumn(object):
  def test_column_init(self, column_config):
    if not column_config['Area']:
      with pytest.raises(TypeError) as areainfo:
        Column(**column_config)
      assert (
          str(
              areainfo.value
          ) == "('Area', 'needs to be either function, numpy array, or float')"
      )
      return

    if not column_config['kappa']:
      with pytest.raises(TypeError) as kappainfo:
        Column(**column_config)
      assert (
          str(kappainfo.value) ==
          "('kappa', 'needs to be either function, numpy array, or float')"
      )
      return

    if not isinstance(column_config['z'],
                      np.ndarray) or not len(column_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Column(**column_config)
      assert (
          str(zinfo.value) == "z needs to be numpy array providing grid levels"
      )
      return

    column = Column(**column_config)

    # The constructor assigns all expected properties
    for k in ['z', 'kappa', 'Area', 'b', 'bs', 'bbot', 'bzbot', 'N2min']:
      assert hasattr(column, k)

    column_signature = funcsigs.signature(Column)

    # The constructor initializes all scalar properties
    # Uses explicit property if present, or the default
    for k in ['bs', 'bbot', 'bzbot', 'N2min']:
      assert getattr(column, k) == (
          column_config[k] if k in column_config and column_config[k] else
          column_signature.parameters[k].default
      )

    # The constructor initializes all vector properties
    # Uses explicit property if present, or the default
    assert all(column.z == column_config['z'])
    assert all(
        column.b == make_array(
            column_config['b'] if 'b' in column_config
            and not column_config['b'] is None else column_signature.
            parameters['b'].default, column.z, 'b'
        )
    )

    # The constructor initializes all z-dependent callable properties
    # Uses explicit property if present, or the default
    for k in ['Area', 'kappa']:
      f = getattr(column, k)
      ft = make_func(column_config[k], column.z, k)
      for z in column.z:
        assert f(z) == ft(z)

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
      assert column.dAkappa_dz(z) == np.gradient(
          column.Area(0) * column.kappa(0), z
      )
    # variable area and kappa
    z = column.z
    Area = lambda l: float(column.Area(0)) / float(l + 1)
    kappa = lambda l: float(column.kappa(0)) * float(l + 1)
    col = Column(z=z, Area=Area, kappa=kappa)
    for z in column.z:
      assert col.Akappa(z) == Area(z) * kappa(z)

  def test_bc(self):
    column = Column(
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
    ya = [column.b[0], column.bz[0]]
    yb = [column.b[-1], column.bz[-1]]
    bc = column.bc(ya, yb)
    assert np.round(bc[0], decimals=2) == -0.01
    assert np.round(bc[1], decimals=2) == -0.02

    column = Column(
        **{
            'Area': 6e13,
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'kappa': 2e-5,
            'bs': 0.05,
            'bbot': 0.02,
            'b': 0.03,
            'N2min': 2e-7
        }
    )
    ya = [column.b[0], column.bz[0]]
    yb = [column.b[-1], column.bz[-1]]
    bc = column.bc(ya, yb)
    assert np.round(bc[0], decimals=2) == 0.01
    assert np.round(bc[1], decimals=2) == -0.02

  def test_ode(self, column):
    column.wA = np.sin
    assert (
        column.ode(column.z, [column.b, column.bz]) == np.vstack((
            column.bz, (np.sin(column.z) - column.dAkappa_dz(column.z)) /
            column.Akappa(column.z) * column.bz
        ))
    ).all()

  def test_solve_equi(self):
    column = Column(
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
    column.wA = np.sin
    sol = integrate.solve_bvp(
        column.ode, column.bc, column.z, [column.b, column.bz]
    )
    sol_values = sol.sol(column.z)
    column.solve_equi(column.wA)
    assert all(
        np.around(column.b, decimals=2) ==
        np.around(np.asarray(np.linspace(-39.95, 0.05, 80)), decimals=2)
    )
    assert all(
        np.around(column.bz, decimals=2) ==
        np.around(sol_values[1, :], decimals=2)
    )

  def test_vertadvdiff(self):
    dt = 60 * 86400

    Area = 6e13
    z = np.asarray(np.linspace(-4000, 0, 80))
    kappa = 2e-5
    b = np.linspace(0.03, -0.002, 80)
    bs = 0.03
    bbot = -0.002
    wA = Area * np.sin(z)
    db_dz = 0.0004 / 50.0
    d2b_dz2 = 0
    dkappa_dz = 0
    dArea_dz = 0
    db_dt1 = (db_dz/
              Area) * (-wA + Area*dkappa_dz + kappa*dArea_dz) + kappa*d2b_dz2

    column = Column(z=z, Area=Area, kappa=kappa, b=b.copy(), bbot=bbot, bs=bs)
    column.vertadvdiff(wA, dt, do_conv=False)

    assert all(
        np.around(column.b[2:-2], decimals=3) ==
        np.around(b[2:-2] - dt * db_dt1[2:-2], decimals=3)
    )

  def test_convect(self):
    N2min = 1.5e-7
    z = np.asarray([-4000.0, -1000.0, -100.0, 0.0])
    b = np.asarray([-0.03, -0.02, 0.01, 0.01])
    bs=0.0
    
    column = Column(z=z, b=b.copy(), bs=bs, N2min=N2min, kappa=2e-5, Area=6e13)
    b[2:] = bs + N2min * (z[2:]-z[1])

    column.convect()
    assert (all(column.b == b))

  def test_horadv(self):
    z = np.asarray([-4000.0, -1000.0, -100.0, 0.0])
    b = np.asarray([-0.03, 0.01, -0.0025, -0.002])
    column = Column(z=z, b=b.copy(), kappa=2e-5, Area=6e13)

    vdx_in = np.asarray([2e8, 2.5e8, 0.0, 0.0])
    b_in = np.asarray([-0.02, 0.01, -0.001, 0.001])
    dt = 60 * 86400

    b[0] = -0.03 + dt*2e6/6e13
    column.horadv(vdx_in, b_in, dt)

    assert (all(np.around(column.b, decimals=4) == np.around(b, decimals=4)))

  def test_timestep(self):
    Area = 6e13
    z = np.asarray(np.linspace(-4000, 0, 80))
    b = np.linspace(-np.sqrt(0.04), 0.0, 80)**2.
    vdx_in = np.asarray([2e4 for n in z])
    b_in = np.asarray([-0.02 for n in z])
    wA = np.sin(z)/Area
    dt = 30 * 86400

    column1 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column2 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column1.timestep(wA=wA, dt=dt)
    column2.vertadvdiff(wA=wA, dt=dt)
    assert (all(column1.b == column2.b))
    column2.horadv(vdx_in=vdx_in, b_in=b_in, dt=dt)
    column2.convect()
    assert (any(column1.b != column2.b))

    column1 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column2 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column1.timestep(wA=wA, dt=dt, b_in=b_in, vdx_in=vdx_in)
    column2.vertadvdiff(wA=wA, dt=dt)
    column2.horadv(vdx_in=vdx_in, b_in=b_in, dt=dt)
    assert (all(column1.b == column2.b))
    column2.convect()
    assert (any(column1.b != column2.b))

    column1 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column2 = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    column1.timestep(wA=wA, dt=dt, b_in=b_in, vdx_in=vdx_in, do_conv=True)
    column2.convect()
    column2.vertadvdiff(wA=wA, dt=dt)
    column2.horadv(vdx_in=vdx_in, b_in=b_in, dt=dt)
    assert (all(column1.b == column2.b))

    column = Column(z=z, b=b.copy(), bs=-0.0, bbot=-0.04, kappa=2e-5, Area=Area)
    with pytest.raises(TypeError) as binfo:
      column.timestep(wA=wA, dt=dt, vdx_in=vdx_in)
    assert (str(binfo.value) == "b_in is needed if vdx_in is provided")
