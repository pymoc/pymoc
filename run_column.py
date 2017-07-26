from model import Model
import math
import numpy as np
from matplotlib import pyplot as plt

N = 4
H_max_so = [2000, 2000, 1500, 1500]
B_int=[3e3,1.2e4, 3e3, 1.2e4]

a=6.37e6
A = 2*np.pi*a**2*59/360*(np.sin(math.radians(69)) - np.sin(math.radians(-48)))


fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

for i in range(0, N):
    m = Model(H_max_so=H_max_so[i], B_int=B_int[i], A=A, diff_type='variable')
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