'''
This script shows an example of how to use the psi_thermwind model class 
to solve for the overturning circulation, given the buoyancy profile
in the basin and in the northern deep water formation region  
'''
import sys
from pymoc.modules import Psi_Thermwind
import numpy as np
from matplotlib import pyplot as plt


# buoyancy profile in the basin:
def b_basin(z):
  return 0.03 * np.exp(z / 300.) - 0.0004


# We will here assume b_N=0 (the default)

z = np.asarray(np.linspace(-4000, 0, 100))

# the next line would turn b_basin from a function to an array,
# which is also a valid input to Model_PsiNA
#b_basin=b_basin(z)

# create column model instance:
m = Psi_Thermwind(z=z, b1=b_basin)
# solve the model:
m.solve()

# Plot results:
fig = plt.figure(figsize=(6, 10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax2.plot(b_basin(m.z), m.z, color='b')
ax1.plot(m.Psi, m.z, color='r')
plt.ylim((-4e3, 0))
ax1.set_xlim((-5, 20))
ax2.set_xlim((-0.01, 0.04))
ax1.set_xlabel('$\Psi$', fontsize=14)
ax2.set_xlabel('b', fontsize=14)
plt.show()
