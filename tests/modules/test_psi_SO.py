import pytest
import inspect
import sys
import numpy as np
sys.path.append('/pymoc/src/modules')
from psi_SO import Psi_SO

@pytest.fixture(scope="module", params=[
  {
  },
  {
    'z': -2000
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': 1e6
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
    'bs': 0.05,
  },
  {
    'z': np.asarray(np.linspace(-4000, 0, 80)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 80),
    'bs': 0.05,
    'tau': 0.12
  },
  # {
  #   'z': np.asarray(np.linspace(-4000, 0, 80)),
  #   'y':
  #   'b':
  #   'bs':
  #   'tau':
  #   'f':
  #   'rho':
  #   'L':
  #   'KGM':
  #   'c':
  #   'bvp_with_EskK':
  #   'Hsill':
  #   'HEk':
  #   'Htapertop':
  #   'Htaperbot':
  #   'smax':
  # }
])
def psi_so_config(request):
  return request.param

@pytest.fixture(scope="module")
def psi_so(request):
  return Psi_SO(**{
    'z': np.asarray(np.linspace(-4000, 0, 81)),
    'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    'b': np.linspace(0.03, -0.001, 81),
    'bs': np.linspace(0.05, 0.10, 51),
    'tau': 0.12
  })

class TestPsi_SO(object):
  def test_psi_so_init(self, psi_so_config):
    if not 'z' in psi_so_config or not isinstance(psi_so_config['z'], np.ndarray) or not len(psi_so_config['z']):
      with pytest.raises(TypeError) as zinfo:
        Psi_SO(**psi_so_config)
      assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
      return
    if not 'y' in psi_so_config or not isinstance(psi_so_config['y'], np.ndarray) or not len(psi_so_config['y']):
      with pytest.raises(TypeError) as yinfo:
        Psi_SO(**psi_so_config)
      assert(str(yinfo.value) == "y needs to be numpy array providing horizontal grid (or boundaries) of ACC")
      return

    # Implicit exceptions thrown if these params are missing, should change these to explicit checks
    for k in ['b', 'bs', 'tau']:
      if not k in psi_so_config:
        with pytest.raises(TypeError) as kinfo:
          Psi_SO(**psi_so_config)
        assert(str(kinfo.value) == "('" + k + "', 'needs to be either function, numpy array, or float')")
        return

    psi_so = Psi_SO(**psi_so_config)
    for k in ['z', 'y', 'b', 'bs', 'tau', 'f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek', 'Hsill', 'HEk', 'Htapertop', 'Htaperbot', 'smax']:
      assert hasattr(psi_so, k)

    psi_so_signature = inspect.signature(Psi_SO)
    
    for k in ['f', 'rho', 'L', 'KGM', 'c', 'bvp_with_Ek', 'Hsill', 'HEk', 'Htapertop', 'Htaperbot', 'smax']:
      assert getattr(psi_so, k)  == (psi_so_config[k] if k in psi_so_config and psi_so_config[k] else psi_so_signature.parameters[k].default)

    for k in ['b']:
      f = getattr(psi_so, k)
      ft = psi_so.make_func(psi_so_config['z'], psi_so_config[k], k)
      for z in psi_so.z:
        assert f(z) == ft(z)

    for k in ['bs', 'tau']:
      f = getattr(psi_so, k)
      ft = psi_so.make_func(psi_so_config['y'], psi_so_config[k], k)
      for y in psi_so.y:
        assert f(y) == ft(y)

  def test_make_func(self, psi_so):
    myst = lambda: 42
    assert psi_so.make_func(psi_so.z, myst, 'myst')() == myst()
    myst = np.arange(0.0, 8.1, 0.1)
    for z in psi_so.z:
      assert psi_so.make_func(psi_so.z, myst, 'myst')(z) == np.interp(z, psi_so.z, myst)
    myst = 6.0
    for y in psi_so.y:
      assert psi_so.make_func(psi_so.y, myst, 'myst')(y) == myst
    myst = 1
    with pytest.raises(TypeError) as mystinfo:
      psi_so.make_func(psi_so.z, myst, 'myst')
    assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  def test_ys(self, psi_so):
    assert(np.round(psi_so.ys(0.02), decimals=3) == np.round(psi_so.y[0] - 1e3, decimals=3))
    assert(np.round(psi_so.ys(0.2), decimals=3) == np.round(psi_so.y[-1], decimals=3))
    for i in range(len(psi_so.y)):
      assert(np.round(psi_so.ys(psi_so.bs(psi_so.y[i])), decimals=3) == np.round(psi_so.y[i], decimals=3))

  def test_calc_N2(self, psi_so):
    N2 = (psi_so.b(psi_so.z[1]) - psi_so.b(psi_so.z[0])) / (psi_so.z[1] - psi_so.z[0])
    for z in psi_so.z:
      assert(np.round(psi_so.calc_N2()(z), decimals=10) == np.round(N2, decimals=10))

  def test_Ekman(self, psi_so):
    # Constant wind stress
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[-1] = 0
    assert(all(np.around(ekman, decimals=3) == np.around(psi_so.calc_Ekman(), decimals=3)))

    # Variable wind stress
    b = np.asarray(np.linspace(0.03, 0.01, 21))
    bs = np.asarray(np.linspace(0.01, 0.02, 51))
    y = np.asarray(np.linspace(0, 2e6, 51))
    z = np.asarray(np.linspace(-4000, 0, 21))
    tau = np.asarray(np.linspace(0.2, 0.12, 51))
    psi = Psi_SO(y=y, z=z, b=b, bs=bs, tau=tau)
    tau_ave = np.zeros((len(z)))
    tau_ave[0:11] = tau[-1]
    tau_ave[11:-1] = [np.mean(tau[-i*5 - 1:]) for i in range(1, 10)]
    ekman = psi.L*tau_ave/(psi.f*psi.rho)
    assert(all(np.around(psi.calc_Ekman(), decimals=3) == np.around(ekman, decimals=3)))

    # With sill tapering
    psi_so.Hsill = 1000.0
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[:20] = [ekman[-1] * (1 - (-3000 - z)**2/1.0e6) for z in psi_so.z[:20]]
    ekman[-1] = 0
    assert(all(np.around(ekman, decimals=3) == np.around(psi_so.calc_Ekman(), decimals=3)))
    psi_so.Hsill = None

    # With constrained ekman layer
    psi_so.HEk = 200
    ekman = (psi_so.L * 0.12) / (psi_so.f * psi_so.rho)
    ekman = [ekman for _ in range(len(psi_so.z))]
    ekman[77:] = [ekman[0] * (1 - (z + 200.0)**2/4.0e4) for z in psi_so.z[77:]]
    assert(all(np.around(ekman, decimals=3) == np.around(psi_so.calc_Ekman(), decimals=3)))
    psi_so.HEk = None

  def test_calc_GM_no_bvp(self, psi_so):
    # Test constant isopycnal sloping
    dy_atz = 2001000.0
    GM = [psi_so.L*psi_so.KGM*z/dy_atz for z in psi_so.z]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert(all(np.around(GM, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))

    # Test minimum dy limit
    dy_atz = 0.1
    psi_so.b = psi_so.make_func(psi_so.z, np.linspace(0.3, 0.2, 81), 'b')
    GM = [-psi_so.L*psi_so.KGM*psi_so.smax for z in psi_so.z]
    GM[-1] = 0
    assert(all(np.around(GM, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))

    # Test maximum isopycnal slope limit
    dy_atz = 2001000.0
    psi_so.smax = -0.01
    psi_so.b = psi_so.make_func(psi_so.z, np.linspace(0.03, -0.001, 81), 'b')
    GM = [psi_so.L*psi_so.KGM*0.01 for _ in psi_so.z]
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert(all(np.around(GM, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))

    # Test quadratic bottom tapering of streamfunction
    psi_so.smax = 0.01
    psi_so.Htaperbot = 800.0
    dy_atz = 2001000.0
    GM = [psi_so.L*psi_so.KGM*z/dy_atz for z in psi_so.z]
    GM[0:17] = [(1.0 - (-3200.0 - psi_so.z[i])**2/6.4e5)*GM[i] for i in range(0,17)] 
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert(all(np.around(GM, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))
    psi_so.Htaperbot = None

    # Test quadratic top tapering of streamfunction
    psi_so.smax = 0.01
    psi_so.Htapertop = 500.0
    dy_atz = 2001000.0
    GM = [psi_so.L*psi_so.KGM*z/dy_atz for z in psi_so.z]
    GM[-11:] = [(1.0 - (500.0 + psi_so.z[i])**2/2.5e5)*GM[i] for i in range(70,81)] 
    psi_so.Psi_Ek = psi_so.calc_Ekman()
    assert(all(np.around(GM, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))
    psi_so.Htapertop = None

    # Test ekman streamfunction limiting for non-outcropping isopycnals
    psi_so = Psi_SO(**{
      'z': np.asarray(np.linspace(-4000, 0, 81)),
      'y': np.asarray(np.linspace(0, 2.0e3, 51)),
      'b': np.linspace(0.03, -0.001, 81),
      'bs': np.linspace(0.05, 0.10, 51),
      'tau': 0.12
    })
    dy_atz = 2001000.0
    GM = [-psi_so.L*psi_so.KGM*0.01 for _ in psi_so.z]
    GM[-1] = 0
    psi_ek = np.asarray([gm + np.abs(gm/2.0) for gm in GM])
    psi_so.Psi_Ek = psi_ek
    assert(all(np.around(-psi_ek*1e6, decimals=3) == np.around(psi_so.calc_GM(), decimals=3)))


    # eps = 0.1 # minimum dy (in meters) (to avoid div. by 0)
    # for ii in range(0, np.size(psi_so.z)):
    #   dy_atz[ii] = max(psi_so.y[-1] - psi_so.ys(psi_so.b(psi_so.z[ii])), eps)
    # print(dy_atz)
    # GM[dy_atz > psi_so.y[-1] - psi_so.y[0]]jk = np.maximum(
    #       GM[dy_atz > psi_so.y[-1] - psi_so.y[0]],
    #     -psi_so.Psi_Ek[dy_atz > psi_so.y[-1] - psi_so.y[0]]*1e6)
    # GM = np.asarray(GM)
    # print(GM)
    # psi_so.Psi_Ek = psi_so.calc_Ekman()
    # print(psi_so.Psi_Ek)
    # print(psi_so.calc_GM())
