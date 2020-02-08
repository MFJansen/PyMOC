#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 15:32:58 2019

@author: malte
"""
from matplotlib import pyplot as plt
import numpy as np

# choose whether to plot dependence on kappa or b_AABW:
do='b_AABW'

if do=='Kappa':
   cases=[('kapx0p1.npz','$0.1\kappa_0$',0.1),
          ('kapx0p25.npz','$0.25\kappa_0$',0.25),
          ('kapx0p5.npz','$0.5\kappa_0$',0.5),
          ('kapx0p75.npz','$0.75\kappa_0$',0.75),
          ('ref.npz','$1\kappa_0$',1.0),
          ('kapx1p25.npz','$1.25\kappa_0$',1.25),
          ('kapx1p5.npz','$1.5\kappa_0$',1.5),
          ('kapx1p75.npz','$1.75\kappa_0$',1.75),
          ('kapx2.npz','$2\kappa_0$',2.0),
          ('kapx2p25.npz','$2.25\kappa_0$',2.25),
          ('kapx2p5.npz','$2.5\kappa_0$',2.5),
          ('kapx2p75.npz','$2.75\kappa_0$',2.75),
          ('kapx3.npz','$3\kappa_0$',3.0)]
   refcase=3       
elif do=='b_AABW':
   cases=[('b_AABW_0.npz','0',0.0),
          ('b_AABW_0003.npz','-0.3',-0.0003),
          ('b_AABW_0006.npz','-0.6',-0.0006),
          ('b_AABW_0009.npz','-0.3',-0.0009),
          ('b_AABW_0012.npz','-1.2',-0.0012),
          ('b_AABW_0015.npz','-1.2',-0.0015),
          ('b_AABW_0018.npz','-1.8',-0.0018),
          ('b_AABW_0024.npz','-2.4',-0.0024),
          ('b_AABW_0030.npz','-3.0',-0.0030),
          ('b_AABW_0036.npz','-3.6',-0.0036),
          ('b_AABW_0042.npz','-4.2',-0.0042),
          ('b_AABW_0048.npz','-4.8',-0.0048),
          ('b_AABW_0060.npz','-6.0',-0.0060)]
   refcase=6
else:
    raise ValueError('Can\'t do {}'.format(do))

# path where to save figure
path='/Users/malte/Dropbox/Column_model_paper/2Basin/Figs-pymoc/kappa_tanh_new/'
path=path+do

plt.close('all')

z=np.asarray(np.linspace(-4000, 0, 80))

A_Atl=7e13 
A_north=5.5e12
A_Pac=1.7e14

fig = plt.figure(figsize=(10,7.5))
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

b_Atl=1.0*np.load(cases[refcase][0])['b_Atl']
b_Pac=1.0*np.load(cases[refcase][0])['b_Pac']
bgrid_ZOC=1.0*np.load(cases[refcase][0])['bgrid_ZOC']
b_basin_grid=(A_Atl*b_Atl+A_Pac*b_Pac)/(A_Atl+A_Pac);

b_NADW=[]
b_SA=[]
b_SO=[]
z_NADW=[]
z_SA=[]
z_SO=[]
Psimax_N=[]
Psimax_SO=[]
Psi_NADW_bSO=[]
Psi_SA_bSO=[]
Psi_dia_max=[]

for case in cases:
   Psi_SO_Atl=1.0*np.load(case[0])['Psi_SO_Atl']
   Psi_SO_Pac=1.0*np.load(case[0])['Psi_SO_Pac']
   Psi_ZOC=1.0*np.load(case[0])['Psi_ZOC']
   Psib_ZOC=1.0*np.load(case[0])['Psib_ZOC']
   Psi_AMOC=1.0*np.load(case[0])['Psi_AMOC']
   Psib_AMOC=1.0*np.load(case[0])['Psib_AMOC']
   b_Atl=1.0*np.load(case[0])['b_Atl']
   b_Pac=1.0*np.load(case[0])['b_Pac']
   bgrid_ZOC=1.0*np.load(case[0])['bgrid_ZOC']
   bgrid_AMOC=1.0*np.load(case[0])['bgrid_AMOC']
   
   
   #b_basin=(A_Atl*b_Atl+A_Pac*b_Pac)/(A_Atl+A_Pac);
   
   Psi_SA_b=( np.interp(b_basin_grid,b_Atl,Psi_SO_Atl,left=np.NaN,right=np.NaN)
                -np.interp(b_basin_grid,bgrid_ZOC,Psib_ZOC,left=np.NaN,right=np.NaN) )
   #Psi_SA_res=( Psi_SO_Atl -np.interp(b_Atl,bgrid_ZOC,Psib_ZOC) )
   Psi_SO_b=( np.interp(b_basin_grid,b_Atl,Psi_SO_Atl,left=np.NaN,right=np.NaN)
             +np.interp(b_basin_grid,b_Pac,Psi_SO_Pac,left=np.NaN,right=np.NaN) )
   Psi_AMOC_b=np.interp(b_basin_grid,bgrid_AMOC,Psib_AMOC,left=np.NaN,right=np.NaN)
   
   #b_NADW.append(np.min(bgrid_AMOC[Psib_AMOC>1.0]))
   b_NADW.append(np.interp(1.0,Psib_AMOC[bgrid_AMOC<0.005],bgrid_AMOC[bgrid_AMOC<0.005]))
   #b_SA.append(np.min(b_basin_grid[Psi_SA_b>0.0]))
   where=np.logical_and(z<-100.,np.isfinite(Psi_SA_b))
   b_SA.append(np.interp(0.0,Psi_SA_b[where],b_basin_grid[where]))
   #b_SO.append(np.min(b_basin_grid[Psi_SO_b>0.0]))
   where=np.logical_and(z<-100.,np.isfinite(Psi_SO_b))
   b_SO.append(np.interp(0.0,Psi_SO_b[where],b_basin_grid[where]))
   #z_NADW.append(np.min(z[Psi_AMOC>1.0]))
   z_NADW.append(np.interp(1.0,Psi_AMOC[z<-200.],z[z<-200.]))
   #z_SA.append(np.min(z[Psi_SO_Atl-Psi_ZOC>0.0]))
   z_SA.append(np.interp(0.0,Psi_SO_Atl[z<-100.]-Psi_ZOC[z<-100.],z[z<-100.]))
   #z_SO.append(np.min(z[Psi_SO_Atl+Psi_SO_Pac>0.0]))
   z_SO.append(np.interp(0.0,Psi_SO_Atl[z<-100.]+Psi_SO_Pac[z<-100.],z[z<-100.]))
   
   Psimax_N.append(np.nanmax(Psib_AMOC))
   Psimax_SO.append(np.nanmax(Psi_SO_b))
   Psi_dia_max.append(np.nanmax( Psi_AMOC_b-np.minimum(np.maximum(Psi_SO_b,0),Psi_AMOC_b) ))
   Psi_NADW_bSO.append(np.interp(b_SO[-1],bgrid_AMOC,Psib_AMOC))
   Psi_SA_bSO.append(np.interp(b_SO[-1],b_basin_grid,Psi_SA_b))
   
   ax1.plot(Psi_SO_Atl+Psi_SO_Pac,z)
   ax2.plot(Psi_SO_Atl-Psi_ZOC,z)
   ax3.plot(Psi_AMOC,z)
   ax4.plot(Psi_SO_b,z)
   ax5.plot(Psi_SA_b,z)
   ax6.plot(Psi_AMOC_b,z)
   if do=='Kappa':
      ax4.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005],b_basin_grid,z))
      ax4.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005])
      ax5.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005],b_basin_grid,z))
      ax5.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005])
      ax6.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005],b_basin_grid,z))
      ax6.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005])
   else:
      ax4.set_yticks(np.interp([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005],b_basin_grid,z))
      ax4.set_yticklabels([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005])
      ax5.set_yticks(np.interp([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005],b_basin_grid,z))
      ax5.set_yticklabels([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005])
      ax6.set_yticks(np.interp([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005],b_basin_grid,z))
      ax6.set_yticklabels([0.02, 0.005,0.002, 0., -0.001, -0.002, -0.003, -0.004, -0.005])
   plt.pause(0.01)
  
ax1.legend([case[1] for case in cases],fontsize=10,frameon=False) 
ax1.plot(0.*z, z,color='black',linewidth=0.5)
ax1.set_xlabel('$\Psi$ [Sv]')
ax1.set_ylabel('Depth [m]')
ax1.set_title('Southern Ocean, z')
ax2.plot(0.*z, z,color='black',linewidth=0.5)
ax2.set_xlabel('$\Psi$ [Sv]')
ax2.set_ylabel('Depth [m]')
ax2.set_title('South Atlantic, z')
ax3.plot(0.*z, z,color='black',linewidth=0.5)
ax3.set_xlabel('$\Psi$ [Sv]')
ax3.set_ylabel('Depth [m]')
ax3.set_title('North Atlantic, z')
ax4.plot(0.*z, z,color='black',linewidth=0.5)
ax4.set_xlabel('$\Psi$ [Sv]')
ax4.set_ylabel('b [m s$^{-2}$]')
ax4.set_title('Southern Ocean, b')
ax5.plot(0.*z, z,color='black',linewidth=0.5)
ax5.set_xlabel('$\Psi$ [Sv]')
ax5.set_ylabel('b [m s$^{-2}$]')
ax5.set_title('South Atlantic, b')
ax6.plot(0.*z, z,color='black',linewidth=0.5)
ax6.set_xlabel('$\Psi$ [Sv]')
ax6.set_ylabel('b [m s$^{-2}$]')
ax6.set_title('North Atlantic, b')
ax1.set_ylim((-4e3,0))
ax2.set_ylim((-4e3,0))
ax3.set_ylim((-4e3,0))
ax4.set_ylim((-4e3,0))
ax5.set_ylim((-4e3,0))
ax6.set_ylim((-4e3,0))
fig.tight_layout()

fig = plt.figure(figsize=(10,6.5))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
x=[case[2] for case in cases];
labels=[case[1] for case in cases]
ax1.plot(x,np.asarray(b_NADW)*1e3)
ax1.plot(x,np.asarray(b_SA)*1e3)
ax1.plot(x,np.asarray(b_SO)*1e3)
ax1.set_xticks(x[2::2])
ax1.set_xticklabels(labels[2::2])
ax1.legend(['$b_{NADW}$','$b_{0,SA}$','$b_{0,SO}$'],fontsize=10,frameon=False)   
ax1.set_ylabel('b [mm s$^{-2}$]')   
ax1.set_title('Buoyancy overlap')   
ax2.plot(x,z_NADW)
ax2.plot(x,z_SA)
ax2.plot(x,z_SO)
ax2.set_xticks(x[2::2])
ax2.set_xticklabels(labels[2::2])
ax2.legend(['$z_{NADW}$','$z_{0,SA}$','$z_{0,SO}$'],fontsize=10,frameon=False)      
ax2.set_ylabel('Depth [m]')   
ax2.set_title('Depth overlap')   
ax3.plot(x,Psi_dia_max)
ax3.plot(x,Psimax_N)
ax3.set_xticks(x[2::2])
ax3.set_xticklabels(labels[2::2])
ax3.legend(['max($\Psi_{dia}$)','max($\Psi_N$)'],fontsize=10,frameon=False)      
ax3.set_ylabel('$\Psi$ [Sv]')   
ax3.set_title('Diabatic Overturning')   
ax4.plot(x,Psi_NADW_bSO)
ax4.plot(x,Psi_SA_bSO)
ax4.set_xticks(x[2::2])
ax4.set_xticklabels(labels[2::2])
ax4.legend(['$\Psi_N(b_{0,SO})$','$\Psi_{SA}(b_{0,SO})$'],fontsize=10,frameon=False)      
ax4.set_ylabel('$\Psi$ [Sv]')   
ax4.set_title('Overturning overlap')
if do=='b_AABW':
  ax1.set_xlabel('b_AABW [mm s$^{-2}$]') 
  ax2.set_xlabel('b_AABW [mm s$^{-2}$]') 
  ax3.set_xlabel('b_AABW [mm s$^{-2}$]') 
  ax4.set_xlabel('b_AABW [mm s$^{-2}$]') 
fig.tight_layout()

ticks=(0,4,8,-1)

fig = plt.figure(figsize=(10,3))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
x=[case[2] for case in cases];
labels=[case[1] for case in cases]
ax1.plot(x,np.asarray(b_NADW)*1e3)
#ax1.plot(x,np.asarray(b_SA)*1e3)
ax1.plot(x,np.asarray(b_SO)*1e3)
ax1.set_xticks([x[i] for i in ticks])
ax1.set_xticklabels([labels[i] for i in ticks])
ax1.legend(['$b_{NADW}$','$b_{0}$'],fontsize=14,frameon=False)   
ax1.set_ylabel('b [mm s$^{-2}$]',fontsize=13)   
ax1.set_title('Bottom Buoyancy')   
ax2.plot(x,z_NADW)
#ax2.plot(x,z_SA)
ax2.plot(x,z_SO)
ax2.set_xticks([x[i] for i in ticks])
ax2.set_xticklabels([labels[i] for i in ticks])
ax2.legend(['$z_{NADW}$','$z_{0}$'],fontsize=14,frameon=False)      
ax2.set_ylabel('Depth [m]',fontsize=12)   
ax2.set_title('Depth')   
ax3.plot(x,Psi_NADW_bSO)
ax3.plot(x,Psi_SA_bSO)
ax3.set_xticks([x[i] for i in ticks])
ax3.set_xticklabels([labels[i] for i in ticks])
ax3.legend(['$\Psi_N(b_{0})$','$\Psi_{SA}(b_{0})$'],fontsize=13,frameon=False)      
ax3.set_ylabel('$\Psi$ [Sv]',fontsize=12)   
ax3.set_title('Overturning overlap')
if do=='b_AABW':
  ax1.set_xlabel('$b_{AABW}$ [mm s$^{-2}$]',fontsize=13) 
  ax2.set_xlabel('$b_{AABW}$ [mm s$^{-2}$]',fontsize=13) 
  ax3.set_xlabel('$b_{AABW}$ [mm s$^{-2}$]',fontsize=13) 
else:
  ax1.set_xlabel('$\kappa$',fontsize=13) 
  ax2.set_xlabel('$\kappa$',fontsize=13) 
  ax3.set_xlabel('$\kappa$',fontsize=13) 
fig.tight_layout()
#plt.savefig(path+'/overlap.png', format='png', dpi=400)


