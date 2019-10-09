# import sys
# import numpy as np
# sys.path.append('/pymoc/src/pymoc/modules')
# from interp_channel import Interpolate_channel
# from matplotlib import pyplot as plt

# @pytest.fixture(scope="module", params=[
#   {
#   },
#   {
#     'z': np.asarray(np.linspace(-4000, 0, 81)),
#   },
#   {
#     'y': np.asarray(np.linspace(0, 2.0e6, 51)),
#   },
#   {
#     'z': np.asarray(np.linspace(-4000, 0, 81)),
#     'y': 1e6
#   },
#   {
#     'z': np.asarray(np.linspace(-4000, 0, 81)),
#     'y': np.asarray(np.linspace(0, 2.0e6, 51)),
#   },
#   {
#     'y': np.asarray(np.linspace(0, 2.0e6, 51)),
#     'z': np.asarray(np.linspace(-4.0e3, 0, 81)),
#     'bn': np.linspace(0.03, -0.01, 81),
#     'bs': np.linspace(0.02, 0.01, 51)
#   },
# ])
# def interp_channel_config(request):
#   return request.param

# @pytest.fixture(scope="module")
# def interp_channel(request):
#   return Interpolate_channel(**{
#     'y': np.asarray(np.linspace(0, 2.0e6, 51)),
#     'z': np.asarray(np.linspace(-4.0e3, 0, 81)),
#     'bn': np.linspace(0.03, 0.01, 81),
#     'bs': np.linspace(0.02, 0.01, 51)
#   })

# class TestInterpolate_channel(object):
#   def test_interp_channel_init(self, interp_channel_config):
#     if not 'y' in interp_channel_config or not isinstance(interp_channel_config['y'], np.ndarray) or not len(interp_channel_config['y']):
#       with pytest.raises(TypeError) as yinfo:
#         Interpolate_channel(**interp_channel_config)
#       assert(str(yinfo.value) == "y needs to be numpy array providing grid levels")
#       return
#     if not 'z' in interp_channel_config or not isinstance(interp_channel_config['z'], np.ndarray) or not len(interp_channel_config['z']):
#       with pytest.raises(TypeError) as zinfo:
#         Interpolate_channel(**interp_channel_config)
#       assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
#       return
#     if not 'bs' in interp_channel_config or not isinstance(interp_channel_config['bs'], np.ndarray):
#       with pytest.raises(TypeError) as bsinfo:
#         Interpolate_channel(**interp_channel_config)
#       assert(str(bsinfo.value) == "('bs', 'needs to be either function, numpy array, or float')")
#       return
#     if not 'bn' in interp_channel_config or not isinstance(interp_channel_config['bn'], np.ndarray):
#       with pytest.raises(TypeError) as bninfo:
#         Interpolate_channel(**interp_channel_config)
#       assert(str(bninfo.value) == "('bn', 'needs to be either function, numpy array, or float')")
#       return

#     interp_channel = Interpolate_channel(**interp_channel_config)
#     for k in ['z', 'y']:
#       assert hasattr(interp_channel, k)

#     bs = interp_channel.make_func(interp_channel_config['bs'], 'bs', interp_channel.y)
#     for y in interp_channel.y:
#       assert(interp_channel.bs(y) == bs(y))

#     bn = interp_channel.make_func(interp_channel_config['bn'], 'bn', interp_channel.z)
#     for z in interp_channel.z:
#       assert(interp_channel.bn(z) == bn(z))

#   def test_make_func(self, interp_channel):
#     myst = lambda: 42
#     assert interp_channel.make_func(myst, 'myst', interp_channel.z)() == myst()
#     myst = np.arange(0.0, 8.1, 0.1)
#     for z in interp_channel.z:
#       assert interp_channel.make_func(myst, 'myst', interp_channel.z)(z) == np.interp(z, interp_channel.z, myst)
#     myst = 6.0
#     for y in interp_channel.y:
#       assert interp_channel.make_func(myst, 'myst', interp_channel.y)(y) == myst
#     myst = 1
#     with pytest.raises(TypeError) as mystinfo:
#       interp_channel.make_func(myst, 'myst', interp_channel.z)
#     assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

#   def test_call_and_gridit(self, interp_channel):
#     l = interp_channel.y[-1]
#     d = interp_channel.z[0]
#     ny = len(interp_channel.y)
#     nz = len(interp_channel.z)
#     barray = interp_channel.gridit()
#     sy = -0.01 / l
#     sz = 0.02 / d
#     for iy in range(0,ny):
#       for iz in range(0,nz):
#         z = interp_channel.z[iz]
#         y = interp_channel.y[iy]
#         b = np.amin([0.03, sz * z + sy * y + 0.02])
#         assert np.abs((b - interp_channel(y, z))/b) < 0.01
#         assert np.abs((barray[iy, iz] - interp_channel(y, z))/b) < 0.01
