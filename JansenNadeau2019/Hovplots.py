#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 10:06:41 2018

@author: malte
"""
import numpy as np
from matplotlib import pyplot as plt

plt.close('all')


dt=86400.*30.                                   # time-step for vert. adv. diff. calc.
total_iters=int(np.ceil(12000.*360*86400./dt))  # total number of timesteps
MOC_up_iters=int(np.floor(1.*360.*86400./dt))   # multiplier for MOC time-step 
Diag_iters=10*MOC_up_iters # multiplier for Diags - needs to be multiple of MOC_up_iters
time=np.arange(total_iters/Diag_iters)*Diag_iters*dt/86400./360.



case=['diags_p006.npz','+3 K','cyan','-']

# get reference solution:
Diag_ref=np.load('diags_ref.npz')
b_basin_ref=1.0*Diag_ref['arr_2'];
z=1.0*Diag_ref['arr_5'];
Psi_ref=1.0*Diag_ref['arr_0'];
#Psi_ref=1.0*Diag_ref['arr_8'];

# get reference solution:
Diag=np.load(case[0])
b_basin=1.0*Diag['arr_2'];
Psi=1.0*Diag['arr_0'];
#Psi=1.0*Diag['arr_8'];


blevs=np.arange(-0.01,0.03,0.001) 
plevs=np.arange(-15.,16,1.0)

time=np.insert(time,0,-50)


psiplot=np.concatenate((Psi_ref[:,:2], Psi[:,1:]), axis=1)
# Plot AMOC Hovmoeller plot
fig = plt.figure(figsize=(10,2.3))
ax1 = plt.subplot2grid((1, 8), (0, 0), colspan=3)
ax2 = plt.subplot2grid((1, 8), (0, 3), colspan=5)
cmap = plt.cm.get_cmap('bwr', 30) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[14] = (1.,1.,1.,1.)
cmaplist[15] = (1.,1.,1.,1.)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
ax1.contour(time,z,psiplot,levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(time,z,psiplot,levels=plevs,cmap=cmap, vmin=-15, vmax=15)
ax1.set_xlim([-20,500])
ax1.set_ylabel('Depth [m]')
ax1.set_xlabel('time [years]')
cmap = plt.cm.get_cmap('bwr', 30) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[14] = (1.,1.,1.,1.)
cmaplist[15] = (1.,1.,1.,1.)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
ax2.contour(time,z,psiplot,levels=plevs,colors='k',linewidths=0.5)
CS=ax2.contourf(time,z,psiplot,levels=plevs,cmap=cmap, vmin=-15, vmax=15)
#plt.colorbar(CS)
fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
ax2.set_xlabel('time [years]')
#ax2.set_title('North Atlantic Overturning')
ax2.set_xlim([500,10000])
ax2.get_yaxis().set_ticklabels([])
fig.subplots_adjust(left=0.1, bottom=0.2, right=0.98, top=0.95,
                wspace=0.1, hspace=None)
#fig.savefig('hovmoeller_p006.png', format='png', dpi=1200)


plevs=np.arange(-15.,16,1.0)#/5.
plevs_cont=np.concatenate((np.arange(-15.,0.,1.0),np.arange(1.,100.,1.0)))#/5.


psiplot=np.concatenate((0.*Psi_ref[:,0][:,np.newaxis], Psi[:,1:]-Psi_ref[:,1:]), axis=1)
# Plot AMOC Hovmoeller plot
fig = plt.figure(figsize=(10,2.3))
ax1 = plt.subplot2grid((1, 8), (0, 0), colspan=3)
ax2 = plt.subplot2grid((1, 8), (0, 3), colspan=5)
cmap = plt.cm.get_cmap('bwr', 30) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[14] = (1.,1.,1.,1.)
cmaplist[15] = (1.,1.,1.,1.)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
ax1.contour(time[1:],z,psiplot,levels=plevs_cont,colors='k',linewidths=0.5)
CS=ax1.contourf(time[1:],z,np.clip(psiplot, min(plevs), max(plevs)),levels=plevs,cmap=cmap, vmin=min(plevs), vmax=max(plevs))
ax1.set_xlim([0,500])
ax1.set_ylabel('Depth [m]')
ax1.set_xlabel('time [years]')
cmap = plt.cm.get_cmap('bwr', 30) 
cmaplist = [cmap(i) for i in range(cmap.N)]
cmaplist[14] = (1.,1.,1.,1.)
cmaplist[15] = (1.,1.,1.,1.)
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
ax2.contour(time[1:],z,psiplot,levels=plevs_cont,colors='k',linewidths=0.5)
CS=ax2.contourf(time[1:],z,np.clip(psiplot, min(plevs), max(plevs)),levels=plevs,cmap=cmap, vmin=min(plevs), vmax=max(plevs))
fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
ax2.set_xlabel('time [years]')
#ax2.set_title('North Atlantic Overturning')
ax2.set_xlim([500,10000])
ax2.get_yaxis().set_ticklabels([])
fig.subplots_adjust(left=0.1, bottom=0.2, right=0.98, top=0.95,
                wspace=0.1, hspace=None)
#fig.tight_layout()   
#fig.savefig('hovmoeller_psi_ano_p006.png', format='png', dpi=600)




