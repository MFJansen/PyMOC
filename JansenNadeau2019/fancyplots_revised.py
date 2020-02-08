#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This script makes fancy looking 2D plots of the column model solution.
It's mostly an exercise in interpolation, but useful to visualize solutions
"""

import sys
sys.path.append('../Modules')
from psi_thermwind import Psi_Thermwind
from psi_SO import Psi_SO
from SO_ML import SO_ML
from interp_channel import Interpolate_channel
from interp_twocol import Interpolate_twocol
import numpy as np
from matplotlib import pyplot as plt

plt.close('all')

Diag=np.load('diags_ref.npz')
adiabatic=False

if adiabatic:
   L=3.*4e6 # zonal length of SO
else:
   L=4e6
kapGM=800. #GM diffusivity in SO (m^2/s)
tau=0.12
l=2.e6;
nb=500

b_basin=1.0*Diag['arr_2'][:,-1];
b_north=1.0*Diag['arr_3'][:,-1];
bs_SO=1.0*Diag['arr_4'][:,-1];
z=1.0*Diag['arr_5'];
y=1.0*Diag['arr_7'];

AMOC = Psi_Thermwind(z=z,b1=b_basin,b2=b_north,f=1.2e-4)
AMOC.solve()
PsiSO=Psi_SO(z=z,y=y,b=b_basin,bs=bs_SO,tau=tau,f=1.2e-4,L=L,KGM=kapGM)
PsiSO.solve()
channel=SO_ML(y=y,L=L,bs=bs_SO)
channel.timestep(b_basin=b_basin,Psi_b=PsiSO.Psi,dt=0.0) #this module doesn't
  # currently have a separate method to just compute Psi_s (although that should
  # probably be implemented) so I'm just calling the time-stepping (with zero time-step))
  # during which Psi_s will be computed


#blevs=np.arange(-0.01,0.03,0.001) 
blevs=np.concatenate((np.arange(-0.01,0.004,0.001),np.arange(0.004,0.04,0.004))) 
plevs=np.arange(-15.,16,1.0)

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
ltrans=1000.
lnorth=400.
lchannel=l/1e3
ybasin=np.linspace(lbasin/60.,lbasin,60)+lchannel;
bbasin=np.tile(b_basin,(len(ybasin),1))
# interpolate buoyancy in northern ransition region:
ytrans=np.linspace(ltrans/20.,ltrans,20)+lchannel+lbasin;
bn=b_north.copy();
bn[0]=b_basin[0];# Notice that the interpolation procedure assumes that the bottom
#buoyancies in both colums match - which may not be exactly the case depending
# on when in teh time-step data is saved 
bint=Interpolate_twocol(y=ytrans*1000.-ytrans[0]*1000.,z=z,bs=b_basin,bn=bn)
btrans=bint.gridit()
# finally set buyancy in northern deep water formation region:
ynorth=np.linspace(lnorth/20.,lnorth,20)+lchannel+lbasin+ltrans;
bnorth=np.tile(bn,(len(ynorth),1))
# now stick it all together:
ynew=np.concatenate((y/1e3,ybasin,ytrans,ynorth))
bnew=np.concatenate((bsouth,bbasin,btrans,bnorth))

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
for iy in range(len(y)+len(ybasin),len(y)+len(ybasin)+len(ytrans)):
    # in the northern transition region, interpolate AMOC.psib to local isopycnal depth
    # for psi_res, while keeping psi_z constant and psi_z constant on non-outcropped isopycnals:
    psiarray_res[iy,:]=np.interp(bnew[iy,:],AMOC.bgrid,AMOC.Psib(nb=nb))
    psiarray_z[iy,:]=AMOC.Psi      
    #psiarray_z[iy,:]=((lchannel+lbasin+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_b[iy,b_basin<bnew[iy,-1]]=AMOC.Psibz(nb=nb)[0][b_basin<bnew[iy,-1]]      
for iy in range(len(y)+len(ybasin)+len(ytrans),len(ynew)):
    # in the northern sinking region, all psi decrease linearly to zero:
    psiarray_res[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psibz(nb=nb)[1])/lnorth
    psiarray_z[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_b[iy,b_basin<bnew[iy,-1]]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psibz(nb=nb)[0][b_basin<bnew[iy,-1]])/lnorth      
psiarray_res[-1,:]=0.;


# plot z-coord. overturning:
fig = plt.figure(figsize=(6.5,3))
ax1 = fig.add_subplot(111)
ax1.plot(np.array([l/1e3,l/1e3]),np.array([z[0],z[-1]]),color='0.5',linewidth=0.7,linestyle='dashed')
CS=ax1.contour(ynew,z,bnew.transpose(),levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10)
CS=ax1.contourf(ynew,z,psiarray_z.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-17, vmax=17)
ax1.contour(ynew,z,psiarray_z.transpose(),levels=plevs,colors='0.5',linewidths=0.5)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Depth-averaged Overturning',fontsize=12)
fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
fig.tight_layout()   
fig.savefig('psi_b_2D_depth_adiabatic_revised.png', format='png', dpi=600)

# plot b-coord. overturning:
fig = plt.figure(figsize=(6.5,3))
ax1 = fig.add_subplot(111)
ax1.plot(np.array([l/1e3,l/1e3]),np.array([z[0],z[-1]]),color='0.5',linewidth=0.7,linestyle='dashed')
CS=ax1.contourf(ynew,z,psiarray_b.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-17, vmax=17)
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
fig.savefig('psi_b_2D_iso_adiabatic_revised.png', format='png', dpi=600)

## plot residual overturning:
#fig = plt.figure(figsize=(6,3))
#ax1 = fig.add_subplot(111)
#CS=ax1.contour(ynew,z,bnew.transpose(),levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#ax1.contour(ynew,z,psiarray_res.transpose(),levels=plevs,colors='k',linewidths=0.5)
#CS=ax1.contourf(ynew,z,psiarray_res.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-17, vmax=17)
#ax1.set_xlim([0,ynew[-1]])
#ax1.set_xlabel('y [km]',fontsize=12)
#ax1.set_ylabel('Depth [m]',fontsize=12)
#ax1.set_title('Residual Overturning',fontsize=12)
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
#fig.tight_layout()   
#       

# Plot profiles:
fig = plt.figure(figsize=(3.5,3.5))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_ylabel('Depth [m]', fontsize=13)
ax1.set_xlim((-8,14))
ax2.set_xlim((-0.016,0.028))
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


# Plot SO ML results:
fig = plt.figure(figsize=(4,3.1))
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
plt.xlim((0,2000))
ax1.set_xlabel('y [km]', fontsize=13)
if adiabatic:
   ax1.set_ylim((-8.,10.))
   ax2.set_ylim((-0.02,0.025))
else:
   ax1.set_ylim((-3.,5.))
   ax2.set_ylim((-0.015,0.025))
ax1.set_ylabel('$\Psi_{SO}$ [SV]', fontsize=13)
ax2.set_ylabel('$b_{SO}$ [m s$^{-2}$]', fontsize=13)
ax1.plot(y/1000.,channel.Psi_s,linewidth=2,color='c',linestyle='--',label='$\Psi_{SO}$')
ax2.plot(y/1000.,channel.bs, linewidth=2,color='r',label='$b_{SO}$')
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(y/1000.,0.*y,linewidth=0.5,color='k',linestyle=':')
fig.tight_layout()   
#fig.savefig('bs_Psi_SO.png', format='png', dpi=600)


