# import pytest
# import sys
# import numpy as np
# sys.path.append('/pymoc/src/pymoc/modules')
# from interp_twocol import Interpolate_twocol
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
#     'y': np.asarray(np.linspace(0, 2.0e6, 81)),
#     'z': np.asarray(np.linspace(-4.0e3, 0, 81)),
#     'bn': np.linspace(0.03, -0.01, 81),
#     'bs': np.linspace(0.02, 0.01, 81)
#   },
# ])
# def interp_twocol_config(request):
#   return request.param

# @pytest.fixture(scope="module")
# def interp_twocol(request):
#   return Interpolate_twocol(**{
#     'y': np.asarray(np.linspace(0, 2.0e6, 81)),
#     'z': np.asarray(np.linspace(-4.0e3, 0, 81)),
#     'bn': np.linspace(0.03, 0.01, 81),
#     'bs': np.linspace(0.02, -0.01, 81)
#   })

# class TestInterpolate_twocol(object):
#   def test_interp_twocol_init(self, interp_twocol_config):
#     if not 'y' in interp_twocol_config or not isinstance(interp_twocol_config['y'], np.ndarray) or not len(interp_twocol_config['y']):
#       with pytest.raises(TypeError) as yinfo:
#         Interpolate_twocol(**interp_twocol_config)
#       assert(str(yinfo.value) == "y needs to be numpy array providing grid levels")
#       return
#     if not 'z' in interp_twocol_config or not isinstance(interp_twocol_config['z'], np.ndarray) or not len(interp_twocol_config['z']):
#       with pytest.raises(TypeError) as zinfo:
#         Interpolate_twocol(**interp_twocol_config)
#       assert(str(zinfo.value) == "z needs to be numpy array providing grid levels")
#       return
#     if not 'bs' in interp_twocol_config or not isinstance(interp_twocol_config['bs'], np.ndarray):
#       with pytest.raises(TypeError) as bsinfo:
#         Interpolate_twocol(**interp_twocol_config)
#       assert(str(bsinfo.value) == "('bs', 'needs to be either function, numpy array, or float')")
#       return
#     if not 'bn' in interp_twocol_config or not isinstance(interp_twocol_config['bn'], np.ndarray):
#       with pytest.raises(TypeError) as bninfo:
#         Interpolate_twocol(**interp_twocol_config)
#       assert(str(bninfo.value) == "('bn', 'needs to be either function, numpy array, or float')")
#       return

#     interp_twocol = Interpolate_twocol(**interp_twocol_config)
#     for k in ['z', 'y']:
#       assert hasattr(interp_twocol, k)

#     bs = interp_twocol.make_func(interp_twocol_config['bs'], 'bs', interp_twocol.z)
#     for z in interp_twocol.z:
#       assert(interp_twocol.bs(z) == bs(z))

#     bn = interp_twocol.make_func(interp_twocol_config['bn'], 'bn', interp_twocol.z)
#     for z in interp_twocol.z:
#       assert(interp_twocol.bn(z) == bn(z))

#   def test_make_func(self, interp_twocol):
#     myst = lambda: 42
#     assert interp_twocol.make_func(myst, 'myst', interp_twocol.z)() == myst()
#     myst = np.arange(0.0, 8.1, 0.1)
#     for z in interp_twocol.z:
#       assert interp_twocol.make_func(myst, 'myst', interp_twocol.z)(z) == np.interp(z, interp_twocol.z, myst)
#     myst = 6.0
#     for y in interp_twocol.y:
#       assert interp_twocol.make_func(myst, 'myst', interp_twocol.y)(y) == myst
#     myst = 1
#     with pytest.raises(TypeError) as mystinfo:
#       interp_twocol.make_func(myst, 'myst', interp_twocol.z)
#     assert(str(mystinfo.value) == "('myst', 'needs to be either function, numpy array, or float')")

#   def test_call_and_gridit(self, interp_twocol):
#     l = interp_twocol.y[-1]
#     d = interp_twocol.z[0]
#     ny = len(interp_twocol.y)
#     nz = len(interp_twocol.z)
#     # btarray = interp_twocol.gridit()
#     barray=np.zeros((ny,nz))
#     btarray=np.zeros((ny,nz))
#     s = 0.001
#     for iy in range(0, ny):
#       for iz in range(0, nz):
#         print("iy="+str(iy)+", iz=", str(iz))
#         zi = interp_twocol.z[nz - iz - 1]
#         yi = interp_twocol.y[nz - iy - 1]
#         z = interp_twocol.z[iz]
#         y = interp_twocol.y[iy]
#         b = interp_twocol.bn(d - zi - s*(l - yi))
#         barray[iy, iz] = b
#         btarray[iy, iz] = interp_twocol(y, z)
#         # assert np.abs((b - interp_twocol(y, z))/interp_twocol(y,z)) <= 0.5
#         # assert np.abs((barray[iy, iz] - interp_twocol(y, z))/b) <= 0.6

#     plt.pcolormesh(interp_twocol.y, interp_twocol.z, barray)
#     plt.colorbar()
#     plt.savefig('btest.png')
#     plt.close()

#     plt.pcolormesh(interp_twocol.y, interp_twocol.z, btarray)
#     plt.colorbar()
#     plt.savefig('test.png')
#     plt.close()

#     plt.pcolormesh(interp_twocol.y, interp_twocol.z, barray - btarray)
#     plt.colorbar()
#     plt.savefig('difftest.png')
#     plt.close()

#     plt.contour(interp_twocol.y, interp_twocol.z, barray, colors='k')
#     plt.contour(interp_twocol.y, interp_twocol.z, btarray, colors='b')
#     plt.savefig('ctest.png')
#     plt.close()
