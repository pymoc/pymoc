import pytest
import sys
import funcsigs
import numpy as np
from scipy import integrate
sys.path.append('/pymoc/src/pymoc/modules')
from SO_ML import SO_ML
from pymoc.utils import make_array


@pytest.fixture(
    scope="module",
    params=[{
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
        'Ks': 100,
        'h': 50,
        'L': 4e6,
        'surflux': 5.9e3,
        'rest_mask': 0.0,
        'b_rest': 0.0,
        'v_pist': 2.0 / 86400.0,
        'bs': 0.02
    }, {
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
    }, {
        'y': None,
    }, {
        'y': 1e6
    }]
)
def so_ml_config(request):
  return request.param


@pytest.fixture(scope="module")
def so_ml(request):
  return SO_ML(y=np.asarray(np.linspace(0, 2.0e6, 51)))


class TestSO_ML(object):
  def test_so_ml_init(self, so_ml_config):
    if not isinstance(so_ml_config['y'],
                      np.ndarray) or not len(so_ml_config['y']):
      with pytest.raises(TypeError) as yinfo:
        SO_ML(**so_ml_config)
      assert (
          str(yinfo.value
              ) == "y needs to be numpy array providing (regular) grid"
      )
      return

    so_ml = SO_ML(**so_ml_config)
    for k in [
        'y', 'Ks', 'h', 'L', 'surflux', 'rest_mask', 'b_rest', 'v_pist',
        'Psi_s', 'bs'
    ]:
      assert hasattr(so_ml, k)

    so_ml_signature = funcsigs.signature(SO_ML)

    for k in ['Ks', 'h', 'L', 'v_pist', 'Psi_s']:
      assert getattr(so_ml, k) == (
          so_ml_config[k] if k in so_ml_config and so_ml_config[k] else
          so_ml_signature.parameters[k].default
      )

    for k in ['surflux', 'rest_mask', 'b_rest', 'bs']:
      assert all(
          getattr(so_ml, k) == make_array((
              so_ml_config[k] if k in so_ml_config and so_ml_config[k] else
              so_ml_signature.parameters[k].default
          ), so_ml.y, k)
      )

  # def test_make_array(self, so_ml):
  #   myst = np.arange(0.0, 8.0, 0.1)
  #   assert all(so_ml.make_array(myst, 'myst') == myst)
  #   myst = lambda n: 42 + n
  #   assert all(so_ml.make_array(myst, 'myst') == myst(so_ml.y))
  #   myst = 5.0
  #   assert all(so_ml.make_array(myst, 'myst') == 5.0 * np.ones((len(so_ml.y))))
  #   myst = 1
  #   with pytest.raises(TypeError) as mystinfo:
  #     so_ml.make_array(myst, 'myst')
  #   assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

  # def test_solve_equi(self, so_ml):
  #   with pytest.raises(TypeError) as info:
  #     so_ml.solve_equi()
  #   assert (str(info.value) == "This functionality is not yet implemented")

  def test_timestep(self, so_ml_config):
    dt = 60 * 86400
    conf = {
        'y': np.asarray(np.linspace(0, 2.0e6, 51)),
        'Ks': 100,
        'h': 50,
        'L': 4e6,
        'surflux': 5.9e3,
        'rest_mask': 0.0,
        'b_rest': 0.0,
        'v_pist': 2.0 / 86400.0,
        'bs': 0.02
    }
    b_basin = np.linspace(0.03, -0.002, 80)
    Psi_b = np.linspace(4.0e6, 0, 80)
    so_ml1 = SO_ML(**conf)
    so_ml2 = SO_ML(**conf)

    with pytest.raises(TypeError) as info:
      so_ml1.timestep(dt=dt, Psi_b=Psi_b)
    assert (
        str(info.value) ==
        'b_basin needs to be numpy array providing buoyancy levels in basin'
    )

    with pytest.raises(TypeError) as info:
      so_ml1.timestep(dt=dt, b_basin=b_basin)
    assert (
        str(info.value) ==
        'Psi_b needs to be numpy array providing overturning at buoyancy levels given by b_basin'
    )

    so_ml1.timestep(dt=dt, b_basin=b_basin, Psi_b=Psi_b)
    so_ml2.advdiff(b_basin=b_basin, Psi_b=Psi_b, dt=dt)
    assert (all(so_ml1.bs == so_ml2.bs))
    assert (all(so_ml1.Psi_s == so_ml2.Psi_s))

  def test_advdiff(self):
    dt = 60 * 86400
    # dt=60
    y = np.asarray(np.linspace(0, 2.0e6, 51))
    z = np.asarray(np.linspace(-4000, 0, 80))
    tau = 0.12
    Ks = 100
    L = 4e6
    h = 50
    surflux = 5.9e3
    dtheta_dy = 2.0 * np.pi / 2.0e6
    b_basin = np.asarray([0.02 * (n / 2.0e6)**2 for n in y])
    bs = np.asarray([b_basin[-1] * np.cos(n * dtheta_dy) for n in y])
    conf = {
        'y': y,
        'Ks': Ks,
        'h': h,
        'L': L,
        'surflux': surflux,
        'rest_mask': 0.0,
        'b_rest': 0.0,
        'v_pist': 2.0 / 86400.0,
        'bs': bs
    }

    Psi_b = np.asarray(np.linspace(1e4, 2.0e4, 51))
    so_ml = SO_ML(**conf)

    # Explicity calculate the analytical solution fo the above setup
    dbs_dy = np.asarray([
        -dtheta_dy * b_basin[-1] * np.sin(n * dtheta_dy) for n in y
    ])
    d2bs_dy2 = np.asarray([
        -dtheta_dy**2 * b_basin[-1] * np.cos(n * dtheta_dy) for n in y
    ])
    db = -((Psi_b / (h*L)) * dbs_dy + Ks*d2bs_dy2 + surflux/h) * dt

    b = -(so_ml.bs.copy() + db)
    so_ml.advdiff(b_basin, Psi_b, dt)
    assert (
        all([np.abs(b[i] - so_ml.bs[i]) / b[i] < 0.05 for i in range(len(b))])
    )
