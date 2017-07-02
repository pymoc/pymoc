import math
import numpy as np
import scipy
from scipy import integrate
from matplotlib import pyplot as plt

f = 1.2e-4 # coriolis frequency
b_s = 0.025 # surface buoyancy in the basin
a = 6.37e6 # planetary radius (used in basin area calculation)
A = 2*np.pi*a**2*59/360*(np.sin(math.radians(69)) - np.sin(math.radians(-48))) # area of upwelling
kappa = 6e-5 # constant diffusivity
kappa_back = 1e-5
kappa_s = 3e-5
kappa_4k = 3e-4
H_max_so = [2000, 2000, 1500, 1500]
B_int=[3e3, 1.2e4, 3e3, 1.2e4]
zi = np.asarray(np.linspace(-1,0,100)) # vertical grid
psi_so_max = 4.
N = 4

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

def psi_so(z, H_so, H):
    pre = psi_so_max * 1e6 / (f * H**3)
    return pre * np.sin([-np.pi * max(H * x, -H_so) / H_so for x in z])**2

def ode_fun(z, y, p, H_so, alpha, kap, dkap_dz):
    H = p[0]
    return np.vstack((y[1], y[2], y[3], alpha(z,H) * y[3] * (y[0] - psi_so(z, H_so, H) - A * dkap_dz(z, H)/(H**2))))

for i in range(0, N):
    H_so = H_max_so[i]
    # kap = lambda z, H: kappa / (f * H**2)
    # dkap_dz = lambda z, H: 0
    kap = lambda z, H: (kappa_back + kappa_s*np.exp(z*H/100)+kappa_4k*np.exp(-z*H/1000 - 4)) / (H**2 * f) # nondim kappa
    dkap_dz=lambda z, H: (kappa_s/100*np.exp(z*H/100)-kappa_4k/1000*np.exp(-z*H/1000-4)) / (H * f) #nondimensional d(kappa)/dz  
    alpha = lambda z, H: H**2 / (A * kap(0, H))
    b = -b_s/(f**2)
    bz = lambda H: B_int[i] / (f**3 * H**2 * A * kap(0, H))
    ode = lambda z, y, p: ode_fun(z, y, p, H_so, alpha, kap, dkap_dz)
    bc = lambda ya, yb, p: np.array([ ya[0], yb[0], ya[1], ya[3] + bz(p[0]), yb[2] - b/p[0]])
    sol_init = np.zeros((4, zi.size))
    sol_init[0,:] = np.ones((zi.size))
    sol_init[2,:] = -0.1*np.ones((zi.size))
    sol_init[3,:] = -bz(H_max_so[i] + 200) * np.ones((zi.size))
    res = integrate.solve_bvp(ode, bc, zi, sol_init, p=[H_max_so[i] - 200])
    Hf = res.p[0]
    print(res.status)
    print(Hf)
    ax1.plot(res.y[0,:]*f*Hf**3/1e6, res.x*Hf)
    ax2.plot(-res.y[2,:]*f**2*Hf, res.x*Hf)
    plt.ylim((-3e3,0))
    ax1.set_xlim((-5,20))
    ax2.set_xlim((-0.01,0.04))
plt.show()
