#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script shows an example for how to make 2D plots of the
column model solution, using the interpolation tools
"""

import sys
sys.path.append('../Modules')
from pymoc.modules import Psi_Thermwind, Psi_SO
from pymoc.plotting import Interpolate_channel, Interpolate_twocol
import numpy as np
from matplotlib import pyplot as plt

plt.close('all')

Diag=np.load('diags.npz')

# Notice that some of these variables may need to be adjusted to match the
# simulations to be plotted:
L=2e7
l=2e6;
nb=500
b_basin=1.0*Diag['arr_2'][:,-1];
b_north=1.0*Diag['arr_3'][:,-1];
bs_SO=1.0*Diag['arr_4'][:,-1];
z=1.0*Diag['arr_5'];
y=1.0*Diag['arr_7'];
tau=1.0*Diag['arr_9'];
kapGM=1.0*Diag['arr_10'];

AMOC = Psi_Thermwind(z=z,b1=b_basin,b2=b_north,f=1.2e-4)
AMOC.solve()
PsiSO=Psi_SO(z=z,y=y,b=b_basin,bs=bs_SO,tau=tau,f=1.2e-4,L=L,KGM=kapGM)
PsiSO.solve()

blevs=np.arange(-0.01,0.03,0.001) 
plevs=np.arange(-28.,28.,2.0)

bs=1.*bs_SO;bn=1.*b_basin;

if bs[0]>bs[1]:
   # due to the way the time-stepping works bs[0] can be infinitesimally larger 
   # than bs[0] here, which messe up interpolation 
   bs[0]=bs[1]
if bs[0]<bn[0]:
   # Notice that bn[0] can at most be infinitesimally larger than bs[0] 
   # (since bottom water formation from the channel should be happening in this case)
   # but for the interpolation to work, we need it infinitesimally smaller than bs[0]  
   bn[0]=bs[0]; 

# first interpolate buoyancy in channel along constant-slope isopycnals: 
bint=Interpolate_channel(y=y,z=z,bs=bs,bn=bn)
bsouth=bint.gridit()
# buoyancy in the basin is all the same:
lbasin=12000.
lnorth=1000.
lchannel=l/1e3
ybasin=np.linspace(lbasin/60.,lbasin,60)+lchannel;
bbasin=np.tile(b_basin,(len(ybasin),1))
# finally, interpolate buoyancy in the north:
ynorth=np.linspace(lnorth/10.,lnorth,10)+lchannel+lbasin;
bn=b_north.copy();
bn[0]=b_basin[0];# Notice that the interpolation procedure assumes that the bottom
#buoyancies in both colums match - which may not be exactly the case depending
# on when in teh time-step data is saved 
bint=Interpolate_twocol(y=ynorth*1000.-ynorth[0]*1000.,z=z,bs=b_basin,bn=bn)
bnorth=bint.gridit()
# now stick it all together:
ynew=np.concatenate((y/1e3,ybasin,ynorth))
bnew=np.concatenate((bsouth,bbasin,bnorth))

# Compute z-coordinate, b-coordinate and residual overturning streamfunction at all latitudes:
psiarray_b=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_res=np.zeros((len(ynew),len(z))) # "residual" overturning - i.e. isopycnal overturning mapped into z space
psiarray_z=np.zeros((len(ynew),len(z)))  # z-space, "eulerian" overturning
for iy in range(1,len(y)):
    # in the channel, interpolate PsiSO.Psi onto local isopycnal depth:
    psiarray_res[iy,:]=np.interp(bnew[iy,:],b_basin,PsiSO.Psi)
    psiarray_z[iy,:]=psiarray_res[iy,:]
    psiarray_b[iy,b_basin<bs_SO[iy]]=PsiSO.Psi[b_basin<bs_SO[iy]]
for iy in range(len(y),len(y)+len(ybasin)):
    # in the basin, linearly interpolate between Psi_SO and Psi_AMOC:
    psiarray_res[iy,:]=((ynew[iy]-lchannel)*AMOC.Psibz(nb=nb)[0]+(lchannel+lbasin-ynew[iy])*PsiSO.Psi)/lbasin   
    psiarray_z[iy,:]=((ynew[iy]-lchannel)*AMOC.Psi+(lchannel+lbasin-ynew[iy])*PsiSO.Psi)/lbasin  
    psiarray_b[iy,:]=((ynew[iy]-lchannel)*AMOC.Psibz(nb=nb)[0]+(lchannel+lbasin-ynew[iy])*PsiSO.Psi)/lbasin  
for iy in range(len(y)+len(ybasin),len(ynew)):
    # in the north, interpolate AMOC.psib to local isopycnal depth:
    psiarray_res[iy,:]=np.interp(bnew[iy,:],AMOC.bgrid,AMOC.Psib(nb=nb))
    psiarray_z[iy,:]=((lchannel+lbasin+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_b[iy,b_basin<bnew[iy,-1]]=AMOC.Psibz(nb=nb)[0][b_basin<bnew[iy,-1]]      
psiarray_res[-1,:]=0.;


# plot z-coord. overturning:
fig = plt.figure(figsize=(7,4))
ax1 = fig.add_subplot(111)
CS=ax1.contour(ynew,z,bnew.transpose(),levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10)
ax1.contour(ynew,z,psiarray_z.transpose(),levels=plevs,colors='0.5',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_z.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-max(plevs), vmax=max(plevs))
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Depth-averaged Overturning',fontsize=12)
fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
fig.tight_layout()   
#fig.savefig('psi_b_2D_depth.png', format='png', dpi=600)

# plot b-coord. overturning:
fig = plt.figure(figsize=(7,4))
ax1 = fig.add_subplot(111)
CS=ax1.contourf(ynew,z,psiarray_b.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-max(plevs), vmax=max(plevs))
ax1.contour(ynew,z,psiarray_b.transpose(),levels=plevs,colors='0.5',linewidths=0.5)
ax1.set_xlim([0,ynew[-1]])
#ax1.plot(ynew,np.interp(bnew[:,-1],b_basin,z),'k',linewidth=1,colors='0.5')
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005, 0., -0.001 , -0.002, -0.003],b_basin,z))
ax1.set_yticklabels([0.02, 0.005, 0., -0.001 , -0.002, -0.003])
ax1.set_title('Isopycnal Overturning',fontsize=12)
fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
fig.tight_layout()  
#fig.savefig('psi_b_2D_iso.png', format='png', dpi=600)
       

# Plot profiles:
fig = plt.figure(figsize=(4,4))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4.5e3,0))
ax1.set_ylabel('Depth [m]', fontsize=13)
ax1.set_xlim((-28,28))
ax2.set_xlim((-0.028,0.028))
ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
ax1.plot(AMOC.Psi, AMOC.z,linewidth=2,color='m',linestyle='--',label='$\Psi_N$')
ax1.plot(PsiSO.Psi, PsiSO.z,linewidth=2,color='c',linestyle='--',label='$\Psi_{SO}$')
ax2.plot(b_north, z, linewidth=2,color='r',label='$b_N$')
ax2.plot(b_basin, z, linewidth=2,color='b',label='$b_B$')
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(0.*z, z,linewidth=0.5,color='k',linestyle=':')
fig.tight_layout()
#fig.savefig('profiles.png', format='png', dpi=600)

plt.show()
