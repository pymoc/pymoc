"""
One dimensionsional ocean column model at equilibrium.
"""

import math
import numpy as np
import scipy
from scipy import integrate
from matplotlib import pyplot as plt

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
psi_max_so = 4.

N = len(B_int)
z = np.zeros(())
zz = np.linspace(-1,0,100) # vertical grid

def bound_cond(ya, yb, p, bs, f, B, A, kappa): 
    H = p[0]
    return np.array([ya[0], yb[0], ya[1], ya[2] + b_s/f, yb[3] + B/(A*f*kappa)])
def fun(z, y, H, a, A, dkap_dz, psi_so):
    return np.vstack((y[1], y[2], y[3], y[3]*(y[0] - psi_so(z, H).T - A*dkap_dz(z, H))))
def psi(z, H, H_max_so, f, psi_max_so):
    d = list(map(lambda x: (np.pi * max(x*H, -H_max_so) / H_max_so), z))
    return psi_max_so * np.sin(d)**2

for i in range(0, N):
    # zi = np.linspace(0, H_max_so[i], zz.size) # vertical z grid on which to solve the ode
    zi = H_max_so[i] * zz
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
    bz = lambda H: B_int[i]/(A*kappa/H)
    # Southern ocean upwelling
    psi_so = lambda z, H: psi(z, H, H_max_so[i], f, psi_max_so)
    # psi_so = lambda z, H: psi_max_so * 1e6 / f / H[0]**3 * \
    #     np.sin(list(map(lambda x: -np.pi * max(x * H[0], -H_max_so[i]) / H_max_so[i], z)))**2
    # Initial guess for solver:
    # if the solver doesn't converge, it often helps to play a little with the initial guess
    sol_init = np.zeros((4, zi.size))
    sol_init[0,:] = np.ones((zi.size)) * f * H_max_so[i]**2
    sol_init[2,:] = -0.1*np.ones((zi.size)) * f
    sol_init[3,:] = B_int[i] / (A*H_max_so[i] * f**2) * np.ones((zi.size))
    # sol_init = np.ones((4,100)) * np.expand_dims([1,0,-0.1,-bz(H_max_so[i] + 200)], 1)
    # This is the actual differential equation we are solving
    # (see Jansen ????)
    # (y(1)=psi,y(2)=psi_z,y(3)=psi_zz,y(4)=psi_zzz) 
    ode = lambda z, y, H: fun(z,y,H,a,A,dkap_dz,psi_so)
    # Boundary conditions for differential equation:
    # (Notice that the eq. is 4th order, but we are also solving for H
    # hence 4 BCs are needed)
    bc = lambda ya, yb, p: bound_cond(ya, yb, p, b_s, f, B_int[i], A, kappa)
    # This is where the equation actually gets solved:
    # sol = bvp4c(ode,bc,solinit);
    res = integrate.solve_bvp(ode, bc, zi, sol_init, p=[-H_max_so[i] - 200], max_nodes=1e20, verbose=0, tol=1e4)
    # evaluate solution at points in z for plotting
    # y = res.sol(zi)
    print(res.status)
    H = res.p[0]
    print(H)

    # if i == 0:
    # plt.plot(-res.y[0,:]/1E6, res.x)
    print(np.shape(res.y))
    plt.plot(res.y[0,:], res.x)
    # plt.plot(-(f**2)*H*res.y[2,:], -res.x*H)
    # plt.xlim((-5,20))
    plt.ylim((-3e3,0))
plt.show()
        # line(y(1,:)*H^3*f/1e6,zz*H,'linestyle','--','color',[Bint(i)./max(Bint),0.5*Bint(i)./max(Bint),1-Bint(i)./max(Bint)]); hold on;
        # ax1=gca;
        # xlabel('$\Psi_{N} \; [SV]$','interpreter','Latex','fontsize',20)
        # ylabel('depth $[m]$','interpreter','Latex','fontsize',20)
        # ax1_pos = ax1.Position; % position of first axes
        # set(ax1,'xlim',[-5 20],'ylim',[-3000 0],'fontsize',15)
        # ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
        # 'YAxisLocation','right','XTick',(-0.01:0.01:0.04),...
        # 'Color','none','xlim',[-0.01 0.04],'ylim',[-3000 0],'fontsize',15);
        # xlabel(ax2,'$b \; [ms^{-2}]$','interpreter','Latex','fontsize',20)
        # plt.plot() 
    # else:

    # H = sol.parameters % H is the depth of the overturning cell
