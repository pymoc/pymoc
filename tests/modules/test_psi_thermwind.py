import pytest
import sys
import funcsigs
import numpy as np
from numpy import testing
from matplotlib import pyplot as plt
sys.path.append('/pymoc/src/pymoc/modules')
from psi_thermwind import Psi_Thermwind
from pymoc.utils import make_func


@pytest.fixture(
    scope="module",
    params=[
        {},
        {
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'b1': np.linspace(0.03, -0.001, 80),
            'b2': np.linspace(0.02, 0.0, 80),
        },
        {
            'z': np.asarray(np.linspace(-4000, 0, 80)),
            'b1': np.linspace(0.03, -0.001, 80),
            'b2': np.linspace(0.02, 0.0, 80),
            'sol_init': np.ones((2, 80))
        },
    ]
)
def psi_config(request):
  return request.param


@pytest.fixture(scope="module")
def psi(request):
  return Psi_Thermwind(
      **{
          'z': np.asarray(np.linspace(-4000, 0, 80)),
          'b1': np.linspace(0.03, -0.01, 80),
          'b2': np.linspace(0.02, 0.0, 80),
          'sol_init': np.ones((2, 80))
      }
  )


class TestPsi_Thermwind(object):
  def test_psi_thermwind_init(self, psi_config):
    if not 'z' in psi_config or not isinstance(psi_config['z'], np.ndarray
                                               ) or not len(psi_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Psi_Thermwind(**psi_config)
      assert (
          str(zinfo.value) == "z needs to be numpy array providing grid levels"
      )
      return

    psi = Psi_Thermwind(**psi_config)
    for k in ['z', 'f', 'sol_init', 'b1', 'b2']:
      assert hasattr(psi, k)

    psi_signature = funcsigs.signature(Psi_Thermwind)

    for k in ['f']:
      assert getattr(psi, k) == (
          psi_config[k] if k in psi_config and psi_config[k] else
          psi_signature.parameters[k].default
      )

    if not 'sol_init' in psi_config:
      testing.assert_array_equal(psi.sol_init, np.zeros((2, len(psi.z))))

    for k in ['b1', 'b2']:
      f = getattr(psi, k)
      ft = make_func(psi_config[k], psi_config['z'], k)
      for z in psi.z:
        assert f(z) == ft(z)

  def test_bc(self, psi):
    testing.assert_array_equal(psi.bc([1, 2], [3, 4]), np.array([1, 3]))

  def test_ode(self, psi):
    ode = psi.ode(-1e3, [0, 1])
    b1 = 0.0
    b2 = 0.005
    f = 1.2e-4
    assert ode[0] == 1
    testing.assert_approx_equal(ode[1][0], 1.0 / f * (b2-b1))

  def test_solve(self, psi):
    psi.solve()
    d = -4e3
    p = 1e-2 * ((2.5/3.0) * psi.z**3 + 5000.0 * psi.z**2 - d *
                ((2.5/3.0) * d + 5000.0) * psi.z)
    assert all([
        np.abs(p[i] - psi.Psi[i] * 1e6) / np.abs(p[i]) < 0.25
        for i in range(1,
                       len(psi.z) - 1)
    ])

  def test_Psib(self):
    psi = Psi_Thermwind(
        **{
            'z': np.asarray(np.linspace(-4000, 0, 81)),
            'b1': np.linspace(0.03, -0.01, 81),
            'b2': np.linspace(0.02, 0.0, 81),
            'sol_init': np.ones((2, 81))
        }
    )

    d = -4e3
    psi.solve()
    dp_dz = np.gradient(psi.Psi, psi.z[1] - psi.z[0])
    [zom, zop] = psi.z[np.where(np.abs(dp_dz) == np.min(np.abs(dp_dz)))]
    izom = np.where(psi.z == zom)[0][0]
    izop = np.where(psi.z == zop)[0][0]
    zom = int(izom)
    zop = int(izop)
    b_up = np.concatenate(
        (psi.b2(psi.z)[:izom], psi.b1(psi.z)[izom:izop], psi.b2(psi.z)[izop:])
    )

    bgrid = np.linspace(-0.01, 0.03, 500)
    psib = np.zeros((500))
    for i in range(500):
      b = bgrid[i]
      for j in range(len(psi.z)):
        if b < b_up[j]:
          psib[i] += dp_dz[j]
        elif b == b_up[j]:
          psib[i] += dp_dz[j] / 2.0

    psib *= 50.0
    bo = np.where(np.abs(bgrid) == np.min(np.abs(bgrid)))[0][0]
    bf = np.where(np.abs(bgrid - 0.02) == np.min(np.abs(bgrid - 0.02)))[0][0]

    f = np.ones(15) / 15
    psib = np.convolve(psib, f, 'same')
    psib -= 1.5

    assert np.sqrt((psib[bo:bf] - psi.Psib()[bo:bf])**2).mean() < 0.25

  def test_Psibz(self, psi):
    psi = Psi_Thermwind(
        **{
            'z': np.asarray(np.linspace(-4000, 0, 81)),
            'b1': np.linspace(0.03, -0.01, 81),
            'b2': np.linspace(0.02, 0.0, 81),
            'sol_init': np.ones((2, 81))
        }
    )

    d = -4e3
    psi.solve()
    dp_dz = np.gradient(psi.Psi, psi.z[1] - psi.z[0])
    [zom, zop] = psi.z[np.where(np.abs(dp_dz) == np.min(np.abs(dp_dz)))]
    izom = np.where(psi.z == zom)[0][0]
    izop = np.where(psi.z == zop)[0][0]
    zom = int(izom)
    zop = int(izop)
    b_up = np.concatenate(
        (psi.b2(psi.z)[:izom], psi.b1(psi.z)[izom:izop], psi.b2(psi.z)[izop:])
    )

    bgrid = psi.b1(psi.z)
    psib = np.zeros((len(psi.z)))
    for i in range(len(psi.z)):
      b = bgrid[i]
      for j in range(len(psi.z)):
        if b < b_up[j]:
          psib[i] += dp_dz[j]
        elif b == b_up[j]:
          psib[i] += dp_dz[j] / 2.0

    psib *= 50.0
    bo = np.where(np.abs(bgrid - 0.02) == np.min(np.abs(bgrid - 0.02)))[0][0]
    bf = np.where(np.abs(bgrid) == np.min(np.abs(bgrid)))[0][0]
    psib -= 1.5
    assert np.sqrt((psib[bo:bf] - psi.Psibz()[0][bo:bf])**2).mean() < 0.3

    bgrid = psi.b2(psi.z)
    psib = np.zeros((len(psi.z)))
    for i in range(len(psi.z)):
      b = bgrid[i]
      for j in range(len(psi.z)):
        if b < b_up[j]:
          psib[i] += dp_dz[j]
        elif b == b_up[j]:
          psib[i] += dp_dz[j] / 2.0

    psib *= 50.0
    bo = np.where(np.abs(bgrid - 0.02) == np.min(np.abs(bgrid - 0.02)))[0][0]
    bf = np.where(np.abs(bgrid) == np.min(np.abs(bgrid)))[0][0]
    psib -= 1.5
    assert np.sqrt((psib[bo:bf] - psi.Psibz()[1][bo:bf])**2).mean() < 0.3

  def test_update(self, psi):
    b1 = 10.0
    b2 = 50.0
    z = psi.z
    psi.update(b1, b2)

    b1_func = make_func(b1, z, 'b1')
    b2_func = make_func(b2, z, 'b2')

    assert (all(psi.b1(z) == b1_func(z)))
    assert (all(psi.b2(z) == b2_func(z)))
