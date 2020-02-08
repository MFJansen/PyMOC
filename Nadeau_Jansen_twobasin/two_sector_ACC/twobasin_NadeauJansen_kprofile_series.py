'''
This script shows an example of a three column model for the 
overturning circulation in two basins connected to a channel in the south.
Two columns represent the two basins, while the third column represents
the northern sinking region, connected to one of the basins. The meridional
overtunning circulation is computed at the northern end of the first basin
(i.e at the interface to the northern sinking region), and at the 
the southern end of each basin (i.e. at the interface to the channel). 
Moreover a zonal overturning circulation between the two basins in the channel
is computed (using the thermal wind relation - analog to the computation of the 
overturning between the basin and the northern sinking region)
Adiabatic mapping is used for the overturning between the differn columns
'''

import sys
sys.path.append('../../Modules')
from psi_thermwind import Psi_Thermwind
from psi_SO import Psi_SO
from column import Column
import numpy as np
#from matplotlib import pyplot as plt

kapxcases=[0.1,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0]
#kapxcases=[0.1]
bAABWcases=[0.00036,-0.0,-0.0006,-0.0012,-0.0018,-0.0024,-0.0030,
           -0.0036,-0.0042,-0.0048,-0.0060,-0.0072,-0.0084]
#bAABWcases=[0.00036,-0.0]



for kapx in kapxcases:
    for bAABW in bAABWcases:
#        plt.close('all')
         
        diag_file=('./diags/'
                  +'kapx'+"{:.0f}".format(np.floor(kapx))+'p'+"{:.0f}".format(kapx*100%100)
                  +'_bAABW00'+"{:02.0f}".format(-bAABW*1e4)) 
        
        # boundary conditions:
        bs=0.02;#0.015;
        bs_north=0.00036;
        bbot=min(bAABW,bs_north)
        
        
        # S.O. surface boundary conditions and grid:
        y=np.asarray(np.linspace(0,3.e6, 51))
        tau=0.16
        #tau=0.305*np.sin(np.pi*(y-1.11e5)/4.7e6)**2-0.068 
        offset=0.0345*(1-np.cos(np.pi*(5.55e5-1.e5)/8e6));
        bs_SO=(0.0345*(1-np.cos(np.pi*(y-1.e5)/8e6))*(y>5.55e5)
               +(bAABW-offset)/5.55e5*np.maximum(0,5.55e5-y)+offset*(y<5.55e5));

        
        #bs_SO=bbot+(bs-bbot)*y/(y[-1]-y[0])
        # time-stepping parameters:
        dt=86400.*20.                               # time-step for vert. adv. diff. calc.
        MOC_up_iters=int(np.floor(1*360*86400/dt)) # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
        plot_iters= int(np.ceil(2500*360*86400/dt))  # plotting frequency (in iterations)
        total_iters=int(np.ceil(10000*360*86400/dt))# total number of timesteps
        
        # Effective diffusivity profile
        def kappaeff(z): # effective diffusivity profile with tapering in BBL
                #return kapx*( 1e-4*(1.1-np.tanh(np.maximum(z+2000.,0)/700. + np.minimum(z+2000.,0)/1000.))
                #        *(1.-np.maximum(-4000.-z+800.,0.)/800.)**2 ) 
                return kapx*( 1e-4*(1.1-np.tanh(np.maximum(z+2000.,0)/1000. + np.minimum(z+2000.,0)/1300.))
                        *(1.-np.maximum(-4000.-z+600.,0.)/600.)**2 ) 
        
                
        #1.5 =10Sv overlap        
              
        A_Atl=7e13#6e13  
        A_north=5.5e12
        A_Pac=1.7e14#1.44e14
        
        Lx = 1.3e+07  #(length of the channel)
        Latl=6./21.*Lx;
        Lpac=15./21.*Lx;
        K = 1800. # 1500 for variable wind stress % 1700 for const. wind stress
        N2min=2e-7
 
        # create vertical grid:
        z=np.asarray(np.linspace(-4000, 0, 80))
        
#        fig = plt.figure(figsize=(2.5,4))
#        plt.semilogx(kappaeff(z),z)
#        plt.xlim(1e-5,5e-4)
#        plt.xlabel('$\kappa_{eff}$ [m$^2$s$^{-1}$]',fontsize=12)
#        plt.ylabel('$z$ [m]',fontsize=13)
#        fig.tight_layout() 
        
        
        # create initial guess for buoyancy profile in the Atl
        def b_Atl(z): return bs*np.exp(z/300.)+z/z[0]*bbot
        #def b_Atl(z): return 0.3*bs+z/z[0]*(bbot-0.3*bs)
        def b_Pac(z): return bs*np.exp(z/300.)+z/z[0]*bbot
        
        # create N.A. overturning model instance
        AMOC = Psi_Thermwind(z=z,b1=b_Atl,b2=0.01*b_Atl(z))
        # and solve for initial overturning streamfunction:
        AMOC.solve()
        # map North Atlantic overturning to isopycnal space:
        [Psi_iso_Atl,Psi_iso_N]=AMOC.Psibz()
        
        # create interbasin zonal overturning model instance
        ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=1e-4)
        #ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=2e-6)
        # and solve for initial overturning streamfunction:
        ZOC.solve()
        # map inter-basin overturning to isopycnal space:
        [Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
        
        
        # create S.O. overturning model instance for Atlantic sector
        #SO_Atl=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=4.6e6,KGM=K.,Hsill=500.,HEk=100.,Htapertop=100.,Htaperbot=500.)
        SO_Atl=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Latl,KGM=K)
        SO_Atl.solve()
        
        # create S.O. overturning model instance for Pacific sector
        #SO_Pac=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Lx*2./3.,KGM=K.,Hsill=500.,HEk=100.,Htapertop=100.,Htaperbot=500.)
        SO_Pac=Psi_SO(z=z,y=y,b=b_Pac(z),bs=bs_SO,tau=tau,L=Lpac,KGM=K)
        SO_Pac.solve()
        
        
        # create adv-diff column model instance for Atl
        Atl= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs,bbot=bbot,Area=A_Atl)
        # create adv-diff column model instance for northern sinking region
        north= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs_north,bbot=bbot,Area=A_north,N2min=N2min)
        # create adv-diff column model instance for Pac
        Pac= Column(z=z,kappa=kappaeff,b=b_Pac,bs=bs,bbot=bbot,Area=A_Pac)
        
        
#        # Create figure:
#        fig = plt.figure(figsize=(6,9))
#        ax1 = fig.add_subplot(111)
#        ax2 = ax1.twiny()
#        plt.ylim((-4e3,0))
#        ax1.set_xlim((-20,30))
#        ax2.set_xlim((-0.02,0.030))
        
        # loop to iteratively find equilibrium solution
        for ii in range(0, total_iters):    
           # update buoyancy profile
           # using isopycnal overturning:
           wA_Atl=(Psi_iso_Atl+Psi_zonal_Atl-SO_Atl.Psi)*1e6
           wAN=-Psi_iso_N*1e6
           wA_Pac=(-Psi_zonal_Pac-SO_Pac.Psi)*1e6
           Atl.timestep(wA=wA_Atl,dt=dt)
           north.timestep(wA=wAN,dt=dt,do_conv=True)
           Pac.timestep(wA=wA_Pac,dt=dt)
           
           if ii%MOC_up_iters==0:
              # update overturning streamfunction (can be done less frequently)
              AMOC.update(b1=Atl.b,b2=north.b)
              AMOC.solve()
              [Psi_iso_Atl, Psi_iso_N]=AMOC.Psibz()
              ZOC.update(b1=Atl.b,b2=Pac.b)
              ZOC.solve()
              [Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
              SO_Atl.update(b=Atl.b)
              SO_Atl.solve()
              SO_Pac.update(b=Pac.b)
              SO_Pac.solve()
             
 #          if ii%plot_iters==0:
 #             # Plot current state:
 #             ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=0.5)
 #             #ax1.plot(SO_Atl.Psi, SO_Atl.z, ':b', linewidth=0.5)
 #             #ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=0.5)
 #             ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=0.5)
 #             ax1.plot(SO_Atl.Psi-ZOC.Psi, SO_Atl.z, '--b', linewidth=0.5)
 #             ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=0.5)
 #             ax1.plot(SO_Atl.Psi+SO_Pac.Psi, z, '--c', linewidth=0.5)
 #             ax2.plot(Atl.b, Atl.z, '-b', linewidth=0.5)
 #             ax2.plot(Pac.b, Pac.z, '-g', linewidth=0.5)
 #             ax2.plot(north.b, north.z, '-r', linewidth=0.5)
 #             plt.pause(0.01)
         
        
#        ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
#        #ax1.plot(SO_Atl.Psi, SO_Atl.z, ':b', linewidth=1.5)
#        #ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
#        ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
#        ax1.plot(SO_Atl.Psi-ZOC.Psi, SO_Atl.z, '--b', linewidth=1.5)
#        ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=1.5)
#        ax1.plot(SO_Atl.Psi+SO_Pac.Psi, z, '--c', linewidth=1.5)
#        ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
#        ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
#        ax2.plot(north.b, north.z, '-r', linewidth=1.5)
        
        
        
        if diag_file is not None:
           np.savez(diag_file, b_Atl=Atl.b, b_Pac=Pac.b, b_north=north.b,
                    Psi_SO_Atl=SO_Atl.Psi, Psi_SO_Pac=SO_Pac.Psi,
                    Psi_ZOC=ZOC.Psi, Psib_ZOC=ZOC.Psib(), bgrid_ZOC=ZOC.bgrid,
                    Psi_AMOC=AMOC.Psi, Psib_AMOC=AMOC.Psib(), bgrid_AMOC=AMOC.bgrid)
        print(diag_file+' done')
