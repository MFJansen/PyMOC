#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 10:06:41 2018

@author: malte
"""
import numpy as np
from matplotlib import pyplot as plt

plt.close('all')


dt=86400.*30.                                  # time-step for vert. adv. diff. calc.
total_iters=int(np.ceil(12000.*360*86400./dt))  # total number of timesteps
MOC_up_iters=int(np.floor(1.*360.*86400./dt))  # multiplier for MOC time-step 
Diag_iters=10*MOC_up_iters # multiplier for Diags - needs to be multiple of MOC_up_iters
time=np.arange(total_iters/Diag_iters)*Diag_iters*dt/86400./360.

db=0.0006

cases=[#('diags_m006.npz','-3 K','cyan','-'),
       #('diags_m0006.npz','-0.3 K','blue','-')],
        ('diags_p006.npz','+3 K','blue','-'),
        ('diags_p0006.npz','+0.3 K','red','-')]

fig = plt.figure(figsize=(8,3.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_xlabel('time [yr]', fontsize=12)
ax1.set_ylabel('b\'/|$\Delta$b|', fontsize=12)
ax2.set_xlabel('time [yr]', fontsize=12)
#ax2.set_ylabel('b/$\Delta$b', fontsize=12)
ax1.set_xlim((0,500))
ax2.set_xlim((0,10000))
#ax1.set_ylim((-1.2,1.2))
#ax2.set_ylim((-1.2,1.2))
ax1.set_ylim((0.,1.3))
ax2.set_ylim((0.,1.3))
ax1.set_title('Thermocline', fontsize=12)
ax2.set_title('Abyss', fontsize=12)


# get reference solution:
Diag=np.load('diags_ref.npz')
b_basin_ref=1.0*Diag['arr_2'];
z=1.0*Diag['arr_5'];

nabyss=40; 
ntherm=70;

for case in cases:
    Diag=np.load(case[0])
    b_basin=1.0*Diag['arr_2'];
    b_abyss=np.mean(b_basin[:nabyss,:]-b_basin_ref[:nabyss,:],axis=0)
    b_therm=np.mean(b_basin[ntherm:,:]-b_basin_ref[ntherm:,:],axis=0)
    print b_abyss[-1]
    print b_therm[-1]
    b_abyss=b_abyss/np.abs(b_abyss[-1])
    b_therm=b_therm/np.abs(b_therm[-1])
    b_abyss[0]=0; b_therm[0]=0; # set anomaly at t=0 to 0 - in the diags i't not exactly,
                                # because surface warming has already been applied and conv adjustment has happened 
    ax1.plot(time, b_therm,color=case[2],linewidth=1,linestyle=case[3])
    ax2.plot(time, b_abyss,color=case[2],linewidth=1,linestyle=case[3])
    plt.pause(0.01)    
  
ax1.legend([case[1] for case in cases],fontsize=12,frameon=False)   
ax1.plot(time, 0.*time+1.0,color='black',linewidth=0.5)
ax1.plot(time, 0.*time-1.0,color='black',linewidth=0.5)
ax2.plot(time, 0.*time+1.0,color='black',linewidth=0.5)
ax2.plot(time, 0.*time-1.0,color='black',linewidth=0.5) 
fig.tight_layout()   
fig.savefig('db_vs_t_amplitudes.png', format='png', dpi=600)
    


    

fig = plt.figure(figsize=(8,3.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_xlabel('time [yr]', fontsize=12)
ax1.set_ylabel('b\'/|$\Delta$b|', fontsize=12)
ax2.set_xlabel('time [yr]', fontsize=12)
ax1.set_xlim((0,500))
ax2.set_xlim((0,10000))
ax1.set_title('Thermocline', fontsize=12)
ax2.set_title('Abyss', fontsize=12)
ax1.set_ylim((-0.1,1.3))
ax2.set_ylim((-0.1,1.3))

cases=[('diags_p0006.npz','+0.3 K','red','-','diags_ref.npz'),
#       ('diags.npz','+0.3 K test','y','--','diags_ref.npz'),
       ('diags_p0006_fixPsiSO.npz','+0.3 K fixed $\Psi_{SO}$','g','--','diags_ref.npz'),
#       ('diags_p0006_fix_bSO_psiSO.npz','+0.3 K fixed $\Psi_{SO}$ and  $b_{SO}$','blue','--','diags_ref.npz'),      
       ('diags_p0006_fixPsiN.npz','+0.3 K fixed $\Psi_{N}$','b',':','diags_ref.npz'),
#       ('diags_m0006_fixpsiN.npz','-0.3 K fixed $\Psi_{N}$','blue',':','diags_ref.npz'),      
       ('diags_p0006_fixbSO.npz','+0.3 K fixed $b_{SO}$','magenta','-.','diags_ref.npz')]
 #      ('diags_m0006_fixbSO.npz','-0.3 K fixed $b_{s,SO}$','blue','-.','diags_ref.npz')]      
          
for case in cases:
    b_basin=1.0*np.load(case[0])['arr_2']
    b_basin_ref=1.0*np.load(case[4])['arr_2']
    b_abyss=np.mean(b_basin[1:nabyss,:]-b_basin_ref[1:nabyss,:],axis=0)
    b_therm=np.mean(b_basin[ntherm:,:]-b_basin_ref[ntherm:,:],axis=0)
    #b_abyss=np.mean(b_basin[1:nabyss,:],axis=0)
    #b_therm=np.mean(b_basin[ntherm:,:],axis=0)
    #b_abyss_ref=np.mean(b_basin_ref[1:nabyss,:],axis=0)
    #b_therm_ref=np.mean(b_basin_ref[ntherm:,:],axis=0)
    #print np.abs(b_abyss[-1])
    #print np.abs(b_therm[-1])
    #b_abyss=b_abyss/np.abs(b_abyss[-1])
    #b_therm=b_therm/np.abs(b_therm[-1])
    b_abyss=b_abyss/db
    b_therm=b_therm/db
    #b_abyss=(b_abyss-b_abyss_ref[0])/0.0006
    #b_therm=(b_therm-b_therm_ref[0])/0.0006
    #b_abyss_ref=(b_abyss_ref-b_abyss_ref[0])/0.0006
    #b_therm_ref=(b_therm_ref-b_therm_ref[0])/0.0006
    b_abyss[0]=0; b_therm[0]=0; # set anomaly at t=0 to 0 - in the diags it's not exactly,
                                # because surface warming has already been applied and conv adjustment has happened 
    ax1.plot(time, b_therm,color=case[2],linewidth=1,linestyle=case[3])
    ax2.plot(time, b_abyss,color=case[2],linewidth=1,linestyle=case[3])
    #ax1.plot(time, b_therm_ref,color='k',linewidth=1,linestyle=':')
    #ax2.plot(time, b_abyss_ref,color='k',linewidth=1,linestyle=':')
    plt.pause(0.01)
    
ax1.legend([case[1] for case in cases[:]], fontsize=11,frameon=False)   
ax1.plot(time, 0.*time+1.,color='black',linewidth=0.5)
ax1.plot(time, 0.*time-1.,color='black',linewidth=0.5)
ax2.plot(time, 0.*time+1.,color='black',linewidth=0.5)
ax2.plot(time, 0.*time-1.,color='black',linewidth=0.5)
fig.tight_layout()   
fig.savefig('db_vs_t_sens.png', format='png', dpi=600)


fig = plt.figure(figsize=(8,3.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_xlabel('time [yr]', fontsize=12)
ax1.set_ylabel('b\'/|$\Delta$b|', fontsize=12)
ax2.set_xlabel('time [yr]', fontsize=12)
ax1.set_xlim((0,500))
ax2.set_xlim((0,10000))
ax1.set_title('Thermocline', fontsize=12)
ax2.set_title('Abyss', fontsize=12)
ax1.set_ylim((-0.1,1.5))
ax2.set_ylim((-0.1,1.5))

cases=[('diags_adiabatic_p0006.npz','+0.3 K ','red','-','diags_adiabatic_ref.npz'),
 #      ('diags.npz','+0.3 K test','y','--','diags_adiabatic_ref.npz'),
       ('diags_adiabatic_p0006_fixPsiSO.npz','+0.3 K fixed $\Psi_{SO}$','green','--','diags_adiabatic_ref.npz'),
 #      ('diags_adiabatic_m0006_fixpsiSO.npz','-0.03 K fixed $\Psi_{SO}$','blue','--','diags_adiabatic_ref.npz'),      
       ('diags_adiabatic_p0006_fixPsiN.npz','+0.3 K fixed $\Psi_{N}$','blue',':','diags_adiabatic_ref.npz'),
 #      ('diags_adiabatic_m0006_fixpsiN.npz','-0.03 K fixed $\Psi_{N}$','blue',':','diags_adiabatic_ref.npz'),      
       ('diags_adiabatic_p0006_fixbSO.npz','+0.3 K fixed $b_{SO}$','magenta','-.','diags_adiabatic_ref.npz')]
 #      ('diags_adiabatic_m0006_fixbSO.npz','-0.03 K fixed $b_{s,SO}$','blue','-.','diags_adiabatic_ref.npz')]      

          
for case in cases:
    b_basin=1.0*np.load(case[0])['arr_2']
    b_basin_ref=1.0*np.load(case[4])['arr_2']
    b_abyss=np.mean(b_basin[1:nabyss,:]-b_basin_ref[1:nabyss,:],axis=0)
    b_therm=np.mean(b_basin[ntherm:,:]-b_basin_ref[ntherm:,:],axis=0)
    #print np.abs(b_abyss[-1])
    #print np.abs(b_therm[-1])
    #b_abyss=b_abyss/np.abs(b_abyss[-1])
    #b_therm=b_therm/np.abs(b_therm[-1])
    b_abyss=b_abyss/db
    b_therm=b_therm/db
    b_abyss[0]=0; b_therm[0]=0; # set anomaly at t=0 to 0 - in the diags it's not exactly,
                                # because surface warming has already been applied and conv adjustment has happened 
    ax1.plot(time, b_therm,color=case[2],linewidth=1,linestyle=case[3])
    ax2.plot(time, b_abyss,color=case[2],linewidth=1,linestyle=case[3])
    plt.pause(0.01)

ax1.legend([case[1] for case in cases[:]], fontsize=10,frameon=False)   
ax1.plot(time, 0.*time+1.,color='black',linewidth=0.5)
ax1.plot(time, 0.*time-1.,color='black',linewidth=0.5)
ax2.plot(time, 0.*time+1.,color='black',linewidth=0.5)
ax2.plot(time, 0.*time-1.,color='black',linewidth=0.5)
fig.tight_layout()   
fig.savefig('db_vs_t_sens_adiabatic.png', format='png', dpi=300)

