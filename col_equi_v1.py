"""
One dimensionsional ocean column model at equilibrium.
"""

import math
import numpy as np
import scipy
from scipy import integrate

f = 1.2e-4 # Coriolis parameter (in the north of basin)
b_s = 0.025 # surface buoyancy in the basin
# kappa is the value used for constant diffusivity:
kappa = 6e-5 
#  the folowing three are parameters for 
#  vertically varying diffusivity
#  add or remove comments in main loop to use constant or varying diff.
kappa_back = 1e-5
kappa_s = 3e-5
kappa_4k = 3e-4
a = 6.37e6 # planetary radius (used in basin area calculation)
# Area over which diapycnal upwelling occurs:
A = 2*np.pi*a**2*59/360*(np.sin(math.radians(69)) - np.sin(math.radians(-48)))
#  Diffusive buoyancy loss to abyss (vector to compute multiple cases):
B_int = [3e3, 1.2e4, 3e3, 1.2e4]
# Depth of SO upwelling (vector to compute multiple cases):
H_max_so = [2000, 2000, 1500, 1500]
# Maximum value of SO upwelling (in SV):
psi_max_so = 4

N = len(B_int)
z = np.zeros(())
zz = np.linspace(-1,0,100) # vertical grid

for i in range(0, N):
    zi = np.linspace(0, H_max_so[i], 500) # vertical z grid on which to solve the ode
    kap = lambda z, H: kappa/(H**2/f)
    dkap_dz = lambda z, H: 0
    # For vertically varying kappa, use the next two lines
    # kap = lambda z, H: 1/H**2/f*(kappa_back + kappa_s*np.exp(z*H/100) + kappa_4k*np.exp(-z*H/1000 - 4)) # nondim kapp
    # dkap_dz = lambda z, H: 1/H/f*(kappa_s/100*np.exp(z*H/100) - kappa_4k/1000*np.exp(-z*H/1000 - 4)) # nondimensional d(kappa)/dz  
    a = lambda z, H: 1 / (kap(z, H) * A) * (H**2) # 1/(kappa A) nondimensionalized
    # nondim. surface buoyancy: (still needs to be div. by H, which 
    # will be done in BC itself, so it doesn't need to be function)
    b = -b_s/(f**2)
    # (Nondim.) stratification at bottom of cell:
    bz = lambda H: B_int[i]/(A*kap(-1,H))/(f**3)/(H**2)
    # Southern ocean upwelling
    psi_so = lambda z, H: psi_max_so * 1e6 / f / H[0]**3 * \
        np.sin(list(map(lambda x: -np.pi * max(x * H[0], -H_max_so[i]) / H_max_so[i], z)))**2
    # Initial guess for solver:
    # if the solver doesn't converge, it often helps to play a little with the initial guess
    sol_init = np.ones((4,100)) * np.expand_dims([1,0,-0.1,-bz(H_max_so[i] + 200)], 1)
    # This is the actual differential equation we are solving
    # (see Jansen ????)
    # (y(1)=psi,y(2)=psi_z,y(3)=psi_zz,y(4)=psi_zzz) 
    ode = lambda z, y, H: np.vstack((y[1], y[2] ,y[3], a(z,H)*(y[0]- psi_so(z,H) - A/H**2*dkap_dz(z,H))*y[3]))
    # Boundary conditions for differential equation:
    # (Notice that the eq. is 4th order, but we are also solving for H
    # hence 4 BCs are needed)
    bc = lambda ya, yb, H: [ya[0], ya[1], ya[3] + bz(H), yb[0], yb[2] - b/H] # fixed bottom b_z,psi, and psi_z
    # This is where the equation actually gets solved:
    # sol = bvp4c(ode,bc,solinit);
    res = integrate.solve_bvp(ode, bc, zz, sol_init, [H_max_so[i] + 200])
    # evaluate solution at points in z for plotting
    y = res.sol(zz)
    print(y)
    # H = sol.parameters % H is the depth of the overturning cell
