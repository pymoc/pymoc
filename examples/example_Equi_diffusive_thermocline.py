'''
This script shows an example of how to use the equilibrium column model
It solves for the diffusive overturning circulation with prescribed bottom
depth and urface and bottom buoyancies.     
'''

import sys
from pymoc.modules import Equi_Column
from matplotlib import pyplot as plt

#Define some key parameters:
A = 1.0e14    # Area of the basin
kappa = 5.0e-5    # vertical diffusivity
f = 1e-4    # Coriolis parameter
b_s = 0.02    # surface buoyancy
b_bot = 0.0    # bottom buoyancy
H = 4000.0    # depth of convection

# create column model instance
m = Equi_Column(A=A, b_bot=b_bot, b_s=b_s, f=f, kappa=kappa, H=H, nz=200)
# solve the model:
m.solve()

# plot overturning and buoyancy profile:
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(m.psi, m.z, '-r', label='$\Psi$')
ax2.plot(m.b, m.z, '-b', label='$b$')
plt.ylim((-4e3, 0))
ax1.set_xlim((-5, 20))
ax2.set_xlim((-0.01, 0.04))
ax1.set_xlabel('$\Psi$')
ax2.set_xlabel('$b$')
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, frameon=False)
plt.show()
