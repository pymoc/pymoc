from model import Model
import math
import numpy as np
from matplotlib import pyplot as plt

N = 4 # number of cases to be considered
H_max_so = [2000, 2000, 1500, 1500] # depth extent of upwelling in SO
B_int=[3e3,1.2e4, 3e3, 1.2e4] # abyssal buoyancy loss

# Area of the basin (this is for idealized single-basin model)
a=6.37e6; A = 2*np.pi*a**2*59/360*(np.sin(math.radians(69)) - np.sin(math.radians(-48)))

# The next few lines are a an example for a reasonable vertically varying kappa profile:
# (to use const. kappa, you can instead just provide kappa_const to the Model)    
kappa_back=1e-5
kappa_s=3e-5
kappa_4k=3e-4
kappa = lambda z, H: (kappa_back + kappa_s*np.exp(z*H/100)+kappa_4k*np.exp(-z*H/1000 - 4))  
dkappa_dz = lambda z, H: (kappa_s/100*np.exp(z*H/100)-kappa_4k/1000*np.exp(-z*H/1000-4))   


fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()  

# Main loop which solves column model subject to different parameters:
for i in range(0, N):
    m = Model(H_max_so=H_max_so[i], B_int=B_int[i], A=A, kappa=kappa, dkappa_dz=dkappa_dz)
    res = m.solve()
    H = res['H']
    z = res['z']
    psi = res['psi']
    print(H)
    ax1.plot(psi[0,:], z)
    ax2.plot(psi[2,:], z)
    plt.ylim((-3e3,0))
    ax1.set_xlim((-5,20))
    ax2.set_xlim((-0.01,0.04))

plt.show()
