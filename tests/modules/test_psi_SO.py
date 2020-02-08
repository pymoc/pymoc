import pytest
import funcsigs
import sys
import numpy as np
sys.path.append('/pymoc/src/pymoc/modules')
from psi_SO import Psi_SO
from pymoc.utils import make_func
from matplotlib import pyplot as plt


@pytest.fixture(
    scope="module",
    params=[{}, {
        'z': -2000
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
        'y': 1e6
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
        'b': np.linspace(0.03, -0.001, 80),
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
        'b': np.linspace(0.03, -0.001, 80),
        'bs': 0.05,
    }, {
        'z': np.asarray(np.linspace(-4000, 0, 80)),
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
        'b': np.linspace(0.03, -0.001, 80),
        'bs': 0.05,
        'tau': 0.12
    }]
)
def psi_so_config(request):
  return request.param


@pytest.fixture(scope="module")
def psi_so(request):
  return Psi_SO(
      **{
          'z': np.asarray(np.linspace(-4000, 0, 81)),
          'y': np.asarray(np.linspace(0, 2.0e6, 51)),
          'b': np.linspace(0.03, -0.001, 81),
          'bs': np.linspace(0.05, 0.10, 51),
          'tau': 0.12
      }
  )


class TestPsi_SO(object):
  def test_psi_so_init(self, psi_so_config):
    if not 'z' in psi_so_config or not isinstance(
        psi_so_config['z'], np.ndarray
    ) or not len(psi_so_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Psi_SO(**psi_so_config)
      assert (
          str(zinfo.value) == "z needs to be numpy array providing grid levels"
      )
      return
    if not 'y' in psi_so_config or not isinstance(
        psi_so_config['y'], np.ndarray
    ) or not len(psi_so_config['y']):
      with pytest.raises(TypeError) as yinfo:
        Psi_SO(**psi_so_config)
      assert (
          str(yinfo.value) ==
          "y needs to be numpy array providing horizontal grid (or boundaries) of ACC"
      )
      return

    # Implicit exceptions thrown if these params are missing, should change these to explicit checks
    for k in ['b', 'bs', 'tau']:
      if not k in psi_so_config:
        with pytest.raises(TypeError) as kinfo:
          Psi_SO(**psi_so_config)
        assert (
            str(kinfo.value) == "('" + k +
            "', 'needs to be either function, numpy array, or float')"
        )
        return

    psi_so = Psi_SO(**psi_so_config)
    for k in [
        'z', 'y', 'b', 'bs', 'tau', 'f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek',
        'Hsill', 'HEk', 'Htapertop', 'Htaperbot', 'smax'
    ]:
      assert hasattr(psi_so, k)

    psi_so_signature = funcsigs.signature(Psi_SO)

    for k in [
        'f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek', 'Hsill', 'HEk',
        'Htapertop', 'Htaperbot', 'smax'
    ]:
      assert getattr(psi_so, k) == (
          psi_so_config[k] if k in psi_so_config and psi_so_config[k] else
          psi_so_signature.parameters[k].default
      )

    for k in ['b']:
      f = getattr(psi_so, k)
      ft = make_func(psi_so_config[k], psi_so_config['z'], k)
      for z in psi_so.z:
        assert f(z) == ft(z)

    for k in ['bs', 'tau']:
      f = getattr(psi_so, k)
      ft = make_func(psi_so_config[k], psi_so_config['y'], k)
      for y in psi_so.y:
        assert f(y) == ft(y)

  def test_ys(self, psi_so):
    assert (
        np.round(psi_so.ys(0.02),
                 decimals=3) == np.round(psi_so.y[0] - 1e3, decimals=3)
    )
    assert (
        np.round(psi_so.ys(0.2),
                 decimals=3) == np.round(psi_so.y[-1], decimals=3)
    )
    for i in range(len(psi_so.y)):
      assert (
          np.round(psi_so.ys(psi_so.bs(psi_so.y[i])),
                   decimals=3) == np.round(psi_so.y[i], decimals=3)
      )

  def test_calc_N2(self, psi_so):
    N2 = (psi_so.b(psi_so.z[1]) -
          psi_so.b(psi_so.z[0])) / (psi_so.z[1] - psi_so.z[0])
    for z in psi_so.z:
      assert (
          np.round(psi_so.calc_N2()(z),
                   decimals=10) == np.round(N2, decimals=10)
      )

  def test_Ekman(self, psi_so):
    # Constant wind stress
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[-1] = 0
    assert (
        all(
            np.around(ekman, decimals=3) ==
            np.around(psi_so.calc_Ekman(), decimals=3)
        )
    )

    # Variable wind stress
    b = np.asarray(np.linspace(0.03, 0.01, 21))
    bs = np.asarray(np.linspace(0.01, 0.02, 51))
    y = np.asarray(np.linspace(0, 2e6, 51))
    z = np.asarray(np.linspace(-4000, 0, 21))
    tau = np.asarray(np.linspace(0.2, 0.12, 51))
    psi = Psi_SO(y=y, z=z, b=b, bs=bs, tau=tau)
    tau_ave = np.zeros((len(z)))
    tau_ave[0:11] = tau[-1]
    tau_ave[11:-1] = [np.mean(tau[-i * 5 - 1:]) for i in range(1, 10)]
    ekman = psi.L * tau_ave / (psi.f * psi.rho)
    assert (
        all(
            np.around(psi.calc_Ekman(), decimals=3) ==
            np.around(ekman, decimals=3)
        )
    )

    # With sill tapering
    psi_so.Hsill = 1000.0
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[:20] = [
        ekman[-1] * (1 - (-3000 - z)**2 / 1.0e6) for z in psi_so.z[:20]
    ]
    ekman[-1] = 0
    assert (
        all(
            np.around(ekman, decimals=3) ==
            np.around(psi_so.calc_Ekman(), decimals=3)
        )
    )
    psi_so.Hsill = None

    # With constrained ekman layer
    psi_so.HEk = 200
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[77:] = [
        ekman[0] * (1 - (z + 200.0)**2 / 4.0e4) for z in psi_so.z[77:]
    ]
    assert (
        all(
            np.around(ekman, decimals=3) ==
            np.around(psi_so.calc_Ekman(), decimals=3)
        )
    )
    psi_so.HEk = None

  def test_calc_GM_no_bvp(self, psi_so):
    # Test constant isopycnal sloping
    dy_atz = 2001000.0
    GM = [psi_so.L * psi_so.KGM * z / dy_atz for z in psi_so.z]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert (
        all(
            np.around(GM, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )

    # Test minimum dy limit
    dy_atz = 0.1
    psi_so.b = make_func(np.linspace(0.3, 0.2, 81), psi_so.z, 'b')
    GM = [-psi_so.L * psi_so.KGM * psi_so.smax for z in psi_so.z]
    GM[-1] = 0
    assert (
        all(
            np.around(GM, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )

    # Test maximum isopycnal slope limit
    dy_atz = 2001000.0
    psi_so.smax = -0.01
    psi_so.b = make_func(np.linspace(0.03, -0.001, 81), psi_so.z, 'b')
    GM = [psi_so.L * psi_so.KGM * 0.01 for _ in psi_so.z]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert (
        all(
            np.around(GM, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )

    # Test quadratic bottom tapering of streamfunction
    psi_so.smax = 0.01
    psi_so.Htaperbot = 800.0
    dy_atz = 2001000.0
    GM = [psi_so.L * psi_so.KGM * z / dy_atz for z in psi_so.z]
    GM[0:17] = [(1.0 - (-3200.0 - psi_so.z[i])**2 / 6.4e5) * GM[i]
                for i in range(0, 17)]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert (
        all(
            np.around(GM, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )
    psi_so.Htaperbot = None

    # Test quadratic top tapering of streamfunction
    psi_so.smax = 0.01
    psi_so.Htapertop = 500.0
    dy_atz = 2001000.0
    GM = [psi_so.L * psi_so.KGM * z / dy_atz for z in psi_so.z]
    GM[-11:] = [(1.0 - (500.0 + psi_so.z[i])**2 / 2.5e5) * GM[i]
                for i in range(70, 81)]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert (
        all(
            np.around(GM, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )
    psi_so.Htapertop = None

    # Test ekman streamfunction limiting for non-outcropping isopycnals
    psi_so = Psi_SO(
        **{
            'z': np.asarray(np.linspace(-4000, 0, 81)),
            'y': np.asarray(np.linspace(0, 2.0e3, 51)),
            'b': np.linspace(0.03, -0.001, 81),
            'bs': np.linspace(0.05, 0.10, 51),
            'tau': 0.12
        }
    )
    dy_atz = 2001000.0
    GM = [-psi_so.L * psi_so.KGM * 0.01 for _ in psi_so.z]
    GM[-1] = 0
    psi_ek = np.asarray([gm + np.abs(gm / 2.0) for gm in GM])
    psi_so.Psi_Ek = psi_ek
    assert (
        all(
            np.around(-psi_ek * 1e6, decimals=3) ==
            np.around(psi_so.calc_GM(), decimals=3)
        )
    )

  def test_calc_GM_bvp_with_Ek(self, psi_so):
    c = 1e3
    h = -4000
    Ly = 2.0e6
    Lx = 5e6
    K = 1.0e3
    yso = -1.0e3
    z = np.asarray(np.linspace(h, 0, 80))
    psi_so = Psi_SO(
        **{
            'z': z,
            'y': np.asarray(np.linspace(0, Ly, 51)),
            'b': np.linspace(0.03, -0.001, 80),
            'bs': np.linspace(0.05, 0.04, 51),
            'tau': 0.12,
            'c': c,
            'L': Lx,
            'KGM': K,
            'bvp_with_Ek': True
        }
    )

    psi_Ek = psi_so.calc_Ekman()
    psi_so.Psi_Ek = psi_Ek / 1e6
    psi_Ek = psi_Ek[0]
    GM = psi_so.calc_GM()
    N2 = psi_so.calc_N2()(z)
    N = np.sqrt(np.abs(N2))
    psio = (psi_Ek * np.exp(-N * h / c) + K * Lx * h /
            (Ly-yso)) / (np.exp(-N * h / c) - np.exp(N * h / c))
    psi = psio * (np.exp(N * z / c) -
                  np.exp(-N * z / c)) + K * Lx * z / (Ly-yso)
    assert (
        all([
            np.abs((GM[i] - psi[i]) / GM[i]) < 0.025
            for i in range(1,
                           len(GM) - 1)
        ])
    )

  def test_calc_GM_bvp(self, psi_so):
    c = 1e3
    h = -4000
    Ly = 2.0e6
    Lx = 5e6
    K = 1.0e3
    yso = -1.0e3
    z = np.asarray(np.linspace(h, 0, 80))
    psi_so = Psi_SO(
        **{
            'z': z,
            'y': np.asarray(np.linspace(0, Ly, 51)),
            'b': np.linspace(0.03, -0.001, 80),
            'bs': np.linspace(0.05, 0.04, 51),
            'tau': 0.12,
            'c': c,
            'L': Lx,
            'KGM': K
        }
    )

    psi_so.Psi_Ek = psi_so.calc_Ekman() / 1e6
    GM = psi_so.calc_GM()
    N2 = psi_so.calc_N2()(z)
    N = np.sqrt(np.abs(N2))
    psio = K * Lx * h / ((Ly-yso) * (np.exp(-N * h / c) - np.exp(N * h / c)))
    psi = psio * (np.exp(N * z / c) -
                  np.exp(-N * z / c)) + K * Lx * z / (Ly-yso)
    assert (
        all(
            np.around(GM[1:-1], decimals=1) ==
            np.around(-psi[1:-1], decimals=1)
        )
    )

  def test_solve(self):
    c = 1e3
    h = -4000
    Ly = 2.0e6
    Lx = 5e6
    K = 1.0e3
    yso = -1.0e3
    z = np.asarray(np.linspace(h, 0, 80))
    psi_so_1 = Psi_SO(
        **{
            'z': z,
            'y': np.asarray(np.linspace(0, Ly, 51)),
            'b': np.linspace(0.03, -0.001, 80),
            'bs': np.linspace(0.05, 0.04, 51),
            'tau': 0.12,
            'c': c,
            'L': Lx,
            'KGM': K,
            'bvp_with_Ek': True
        }
    )
    psi_so_2 = Psi_SO(
        **{
            'z': z,
            'y': np.asarray(np.linspace(0, Ly, 51)),
            'b': np.linspace(0.03, -0.001, 80),
            'bs': np.linspace(0.05, 0.04, 51),
            'tau': 0.12,
            'c': c,
            'L': Lx,
            'KGM': K,
            'bvp_with_Ek': True
        }
    )

    Psi_Ek = psi_so_1.calc_Ekman() / 1e6
    psi_so_1.Psi_Ek = Psi_Ek
    Psi_GM = psi_so_1.calc_GM() / 1e6
    Psi = Psi_Ek + Psi_GM
    Psi[0] = 0.
    psi_so_2.solve()
    assert (all(psi_so_2.Psi == Psi))

  def test_update(self, psi_so):
    b = 10.0
    bs = 50.0
    z = psi_so.z
    y = psi_so.y
    psi_so.update(b, bs)

    b_func = make_func(b, z, 'b')
    bs_func = make_func(bs, y, 'bs')

    assert (all(psi_so.b(z) == b_func(z)))
    assert (all(psi_so.bs(y) == bs_func(y)))
