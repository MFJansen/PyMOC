#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:01:22 2019

@author: Malte
"""


import sys
sys.path.append('../Modules')
import numpy as np
from matplotlib import pyplot as plt

kapxcases=[0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
bAABWcases=[0.00036,-0.0,-0.0006,-0.0012,-0.0018,-0.0024,-0.0030,
            -0.0036,-0.0042,-0.0048,-0.0060,-0.0072,-0.0084]

z=np.asarray(np.linspace(-4000, 0, 80))

A_Atl=7e13 
A_north=5.5e12
A_Pac=1.7e14

Nkap=len(kapxcases);NbAABW=len(bAABWcases);

b_NADW=np.zeros([Nkap,NbAABW])
b_SA=np.zeros([Nkap,NbAABW])
b_SO=np.zeros([Nkap,NbAABW])
z_NADW=np.zeros([Nkap,NbAABW])
z_SA=np.zeros([Nkap,NbAABW])
z_SO=np.zeros([Nkap,NbAABW])
Psimax_N=np.zeros([Nkap,NbAABW])
Psimax_SO=np.zeros([Nkap,NbAABW])
Psi_NADW_bSO=np.zeros([Nkap,NbAABW])
Psi_SA_bSO=np.zeros([Nkap,NbAABW])
Psi_dia_max=np.zeros([Nkap,NbAABW])

kapi=0
for kapx in kapxcases:
    bAABWi=0
    for bAABW in bAABWcases:
       #plt.close('all')
         
       diag_file=('./diags/'
                  +'kapx'+"{:.0f}".format(np.floor(kapx))+'p'+"{:.0f}".format(kapx*100%100)
                  +'_bAABW00'+"{:02.0f}".format(-bAABW*1e4)+'.npz') 
    
       Psi_SO_Atl=1.0*np.load(diag_file)['Psi_SO_Atl']
       Psi_SO_Pac=1.0*np.load(diag_file)['Psi_SO_Pac']
       Psi_ZOC=1.0*np.load(diag_file)['Psi_ZOC']
       Psib_ZOC=1.0*np.load(diag_file)['Psib_ZOC']
       Psi_AMOC=1.0*np.load(diag_file)['Psi_AMOC']
       Psib_AMOC=1.0*np.load(diag_file)['Psib_AMOC']
       b_Atl=1.0*np.load(diag_file)['b_Atl']
       b_Pac=1.0*np.load(diag_file)['b_Pac']
       bgrid_ZOC=1.0*np.load(diag_file)['bgrid_ZOC']
       bgrid_AMOC=1.0*np.load(diag_file)['bgrid_AMOC']
             
       b_basin=(A_Atl*b_Atl+A_Pac*b_Pac)/(A_Atl+A_Pac);
       
       Psi_SA_b=( np.interp(b_basin,b_Atl,Psi_SO_Atl,left=np.NaN,right=np.NaN)
                 -np.interp(b_basin,bgrid_ZOC,Psib_ZOC,left=np.NaN,right=np.NaN) )
       Psi_SO_b=( np.interp(b_basin,b_Atl,Psi_SO_Atl,left=np.NaN,right=np.NaN)
                 +np.interp(b_basin,b_Pac,Psi_SO_Pac,left=np.NaN,right=np.NaN) )
       Psi_AMOC_b=np.interp(b_basin,bgrid_AMOC,Psib_AMOC,left=np.NaN,right=np.NaN)
#       
       b_NADW[kapi,bAABWi]=(np.interp(1.0,Psib_AMOC[bgrid_AMOC<0.005],bgrid_AMOC[bgrid_AMOC<0.005]))
       where=np.logical_and(z<-100.,np.isfinite(Psi_SA_b))
       b_SA[kapi,bAABWi]=(np.interp(0.0,Psi_SA_b[where],b_basin[where]))
       where=np.logical_and(z<-100.,np.isfinite(Psi_SO_b))
       b_SO[kapi,bAABWi]=(np.interp(0.0,Psi_SO_b[where],b_basin[where]))
       z_NADW[kapi,bAABWi]=(np.interp(1.0,Psi_AMOC[z<-200.],z[z<-200.]))
       z_SA[kapi,bAABWi]=(np.interp(0.0,Psi_SO_Atl[z<-100.]-Psi_ZOC[z<-100.],z[z<-100.]))
       z_SO[kapi,bAABWi]=(np.interp(0.0,Psi_SO_Atl[z<-100.]+Psi_SO_Pac[z<-100.],z[z<-100.]))
       Psimax_N[kapi,bAABWi]=(np.nanmax(Psib_AMOC))
       Psimax_SO[kapi,bAABWi]=(np.nanmax(Psi_SO_b))
       Psi_dia_max[kapi,bAABWi]=(np.nanmax( Psi_AMOC_b-np.minimum(np.maximum(Psi_SO_b,0),Psi_AMOC_b) ))
       Psi_NADW_bSO[kapi,bAABWi]=(np.interp(b_SO[kapi,bAABWi],bgrid_AMOC,Psib_AMOC))
       Psi_SA_bSO[kapi,bAABWi]=(np.interp(b_SO[kapi,bAABWi],b_basin,Psi_SA_b))
       #print kapi
       #print bAABWi
       bAABWi+=1
    kapi+=1 
      
#******************************************************************************    
# GCM data for comparison
#******************************************************************************
kappa =   [1,1,1,1,1,1,0.25, 0.5, 2, 3]
bAABW =   np.array([-0.7772E-2, -0.6021E-2, -0.4296E-2, -0.2482E-2,-0.1086E-2,0.0391E-2,
                    -0.1588E-2,-0.1405E-2,-0.0085E-2, 0.0841E-2 ])
bNADW =   [-0.0376E-2,-0.0139E-2, 0.0219E-2, 0.0612E-2, 0.0586E-2, 0.0561E-2,
           0.0240E-2, 0.0378E-2, 0.0619E-2, 0.0395E-2]
b0    =   np.array([-0.1335E-2,-0.0544E-2, 0.0219E-2, 0.0872E-2, 0.1145E-2, 0.1606E-2,
                    0.0059E-2, 0.0378E-2, 0.3154E-2, 0.4434E-2])
bSA   =   np.array([-0.1987E-2,-0.1292E-2,-0.0414E-2, 0.0248E-2, 0.0378E-2, 0.0863E-2,
                    -0.0212E-2, 0.0010E-2, 0.1783E-2, 0.2465E-2 ])
zNADW =   np.array([-2047., -2290., -2812., -3532., -3685., -3685.,
                    -2677.,  -3091.,  -3841.,  -3841.])
z0    =   np.array([-1816., -1930., -1930., -2047., -1705., -1816,
                    -1816., -1816., -2047., -2047.])
zSA   =   np.array([-2167., -2290., -2416., -2545., -2416., -3235.,
                    -2047., -2167., -3532., -3532.])
psiN  =   [-0.7189, -0.2000, 1.1471, 7.0356, 12.6424, 14.8468,
           -0.1225, 2.0679, 17.3247, 17.9337]
psiSA =   [2.8234, 3.5982, 5.1028, 7.5824, 9.0466, 11.5372,
           2.3794, 5.0167, 14.2445, 18.1607] 
#*******************************************************************************    
    
levs=np.arange(-28,29,2)   
fig = plt.figure(figsize=(10,3.8))
ax1 = fig.add_subplot(121)
CS=ax1.contourf(bAABWcases,kapxcases,Psi_NADW_bSO,cmap=plt.cm.bwr,levels=levs, vmin=-27, vmax=27)
#ax1.scatter(bAABWcases,kapxcases,Psi_NADW_bSO,cmap=plt.cm.bwr, vmin=-25, vmax=25)
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$\Psi_{N}(b_0)$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1.scatter(bAABW, kappa,c=psiN,cmap=plt.cm.bwr, vmin=-27, vmax=27, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1 = fig.add_subplot(122)
CS=ax1.contourf(bAABWcases,kapxcases,Psi_SA_bSO,cmap=plt.cm.bwr,levels=levs, vmin=-27, vmax=27)
ax1.scatter(bAABW, kappa,c=psiSA,cmap=plt.cm.bwr, vmin=-27, vmax=27, s=100,edgecolors='k')
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$\Psi_{SA}(b_0)$',fontsize=12)
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
fig.colorbar(CS, orientation='vertical')
fig.tight_layout()


levs=np.arange(-0.0044,0.0044,0.0004)   

fig = plt.figure(figsize=(13,3.8))
ax1 = fig.add_subplot(131)
CS=ax1.contourf(bAABWcases,kapxcases,b_NADW,levels=levs, vmin=-0.004, vmax=0.003)
ax1.scatter(bAABW, kappa,c=bNADW, vmin=-0.004, vmax=0.003, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$b_{NADW}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(132)
CS=ax1.contourf(bAABWcases,kapxcases,b_SA,levels=levs, vmin=-0.004, vmax=0.003)
ax1.scatter(bAABW, kappa,c=bSA, vmin=-0.004, vmax=0.003, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$b_{SA}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(133)
CS=ax1.contourf(bAABWcases,kapxcases,b_SO,levels=levs, vmin=-0.004, vmax=0.003)
ax1.scatter(bAABW, kappa,c=b0, vmin=-0.004, vmax=0.003, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$b_{0}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
fig.tight_layout()
#fig.subplots_adjust(bottom=0.2)

levs=np.arange(-0.0048,0.0049,0.0002)   

fig = plt.figure(figsize=(10,3.8))
ax1 = fig.add_subplot(121)
CS=ax1.contourf(bAABWcases,kapxcases,b_SO-b_NADW,cmap=plt.cm.bwr,levels=levs, vmin=-0.0039, vmax=0.0039)
ax1.scatter(bAABW, kappa,c=b0-bNADW,cmap=plt.cm.bwr,vmin=-0.0039, vmax=0.0039, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$b_0-b_{NADW}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(122)
CS=ax1.contourf(bAABWcases,kapxcases,b_SO-b_SA,cmap=plt.cm.bwr,levels=levs,vmin=-0.0039, vmax=0.0039)
ax1.scatter(bAABW, kappa,c=b0-bSA,cmap=plt.cm.bwr,vmin=-0.0039, vmax=0.0039, s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$b_0-b_{SA}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
fig.tight_layout()

levs=np.arange(-4000,-1001,200)   

fig = plt.figure(figsize=(13,3.8))
ax1 = fig.add_subplot(131)
CS=ax1.contourf(bAABWcases,kapxcases,z_NADW,levels=levs, vmin=min(levs), vmax=max(levs))
ax1.scatter(bAABW, kappa,c=zNADW, vmin=min(levs), vmax=max(levs), s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$z_{NADW}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(132)
CS=ax1.contourf(bAABWcases,kapxcases,z_SA,levels=levs, vmin=min(levs), vmax=max(levs))
ax1.scatter(bAABW, kappa,c=zSA, vmin=min(levs), vmax=max(levs), s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$z_{SA}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(133)
CS=ax1.contourf(bAABWcases,kapxcases,z_SO,levels=levs, vmin=min(levs), vmax=max(levs))
ax1.scatter(bAABW, kappa,c=z0, vmin=min(levs), vmax=max(levs), s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$z_{0}$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
fig.tight_layout()
#fig.subplots_adjust(bottom=0.2)

levs=np.arange(-2900,2901,200)   

fig = plt.figure(figsize=(10,3.8))
ax1 = fig.add_subplot(121)
CS=ax1.contourf(bAABWcases,kapxcases,z_SO-z_NADW,cmap=plt.cm.bwr,levels=levs, vmin=min(levs), vmax=max(levs))
ax1.scatter(bAABW, kappa,c=z0-zNADW,cmap=plt.cm.bwr, vmin=min(levs), vmax=max(levs), s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$z_{NADW}-z_0$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
ax1 = fig.add_subplot(122)
CS=ax1.contourf(bAABWcases,kapxcases,z_SO-z_SA,cmap=plt.cm.bwr,levels=levs, vmin=min(levs), vmax=max(levs))
ax1.scatter(bAABW, kappa,c=z0-zSA,cmap=plt.cm.bwr, vmin=min(levs), vmax=max(levs), s=100,edgecolors='k')
ax1.set_xlim(np.min(bAABWcases),np.max(bAABWcases))
ax1.set_ylim(np.min(kapxcases),np.max(kapxcases))
ax1.set_ylabel('$\kappa/\kappa_0$',fontsize=12)
ax1.set_xlabel('$b_{AABW}$',fontsize=12)
ax1.set_title('$z_{SA}-z_0$',fontsize=12)
fig.colorbar(CS, orientation='vertical')
fig.tight_layout()

