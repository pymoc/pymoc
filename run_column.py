from model import Model
import math
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

N = 4 # number of cases to be considered
H_max_so = [2000, 2000, 1500, 1500] # depth extent of upwelling in SO
B_int=[3e3,1.2e4, 3e3, 1.2e4] # abyssal buoyancy loss

psi_so_max=4.
            
# Area of the basin (this is for idealized single-basin model)
a=6.37e6; A = 2*np.pi*a**2*59/360*(np.sin(math.radians(69)) - np.sin(math.radians(-48)))

# The next few lines are a an example for a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar, and don't provide dkappa_dz - it won't be used if you do)    
# (if kappa is defined but dkappa_dz isn't, dkappa_dz will be computet numerically by the model)
kappa_back=1e-5
kappa_s=3e-5
kappa_4k=3e-4
kappa = lambda z: (kappa_back + kappa_s*np.exp(z/100)+kappa_4k*np.exp(-z/1000 - 4))  
dkappa_dz = lambda z: (kappa_s/100*np.exp(z/100)-kappa_4k/1000*np.exp(-z/1000-4))   

fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

# kappa = lambda z: 6e-5 
# kappa_back = 1e-5
# kappa_s = 3e-5
# kappa_4k = 3e-4
# kappa = lambda z: kappa_back + kappa_s*np.exp(z/100.) + kappa_4k*np.exp(-z/1000. - 4.)

# Main loop which solves column model subject to different parameters:
for i in range(0, N):
    psi_so= lambda z: (psi_so_max * 1e6 * np.sin([-np.pi * max(x, -H_max_so[i]) / H_max_so[i] for x in z] )**2 )
    m = Model(B_int=B_int[i], A=A, kappa=kappa, dkappa_dz=dkappa_dz, psi_so=psi_so)
    res = m.solve()
    H = res['H']
    zg = res['z']
    psi = res['psi']
    print(H)
    ax1.plot(psi[0,:], zg)
    ax2.plot(psi[2,:], zg)
    plt.ylim((-3e3,0))
    ax1.set_xlim((-5,20))
    ax2.set_xlim((-0.01,0.04))

plt.show()
