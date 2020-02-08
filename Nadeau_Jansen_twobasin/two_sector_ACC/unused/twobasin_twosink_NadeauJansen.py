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
The SO surface b-profile has two linear parts
one linear dependency from the basin buoyancy to the buoyancy of norther sinking
and a second (applied over the southernmost 5 gridpints matching to a prescribed
bottom buoyancy) - varying the second shows how the abyssal cell (and spec. abyssal strat)
affects the uper cell
'''
import sys
sys.path.append('../Modules')
from psi_thermwind import Psi_Thermwind
from psi_SO import Psi_SO
from column import Column
import numpy as np
from matplotlib import pyplot as plt

# boundary conditions:
bs=0.015; bs_northAtl=0.0; bs_northPac=0.008; bbot= -0.0015 

# S.O. surface boundary conditions and grid:
y=np.asarray(np.linspace(0,2.2e6, 40))
tau=0.12 #*np.sin(np.pi*y/2.e6)
#bs_SO=9*y; bs_SO[2:]= bs*(y[2:]-y[2])/(y[-1]-y[2])
#bs_SO[:2]=bbot*(y[2]-y[:2])/(y[2]-y[0])
bs_SO=bbot+(bs-bbot)*y/(y[-1]-y[0])


# time-stepping parameters:
dt=86400.*30.                                 # time-step for vert. adv. diff. calc.
MOC_up_iters=int(np.floor(2.*360*86400/dt)) # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters= int(np.ceil(500*360*86400/dt))  # plotting frequency (in iterations)
total_iters=int(np.ceil(8000*360*86400/dt)) # total number of timesteps

# The next few lines define a reasonable vertically varying kappa profile:
# (to use const. kappa, simply define kappa as scalar)
def kappaeff(z): # effective diffusivity profile with tapering in BBL
        return ( 5e-5*(1.-np.maximum(-4000.-z+1000.,0.)/1000.)**2 ) 
   

A_Atl=6.1e13  
A_northAtl=5.1e12
A_Pac=1.5e14
A_northPac=1e13

Lx =  1.13e+07  #(length of the channel)
K = 1000
# create vertical grid:
z=np.asarray(np.linspace(-4000, 0, 80))


# create initial guess for buoyancy profile in the Atl
def b_Atl(z): return bs*np.exp(z/300.)+z/z[0]*bbot
#def b_Atl(z): return 0.3*bs+z/z[0]*(bbot-0.3*bs)
def b_Pac(z): return bs*np.exp(z/300.)+z/z[0]*bbot

# create N.A. overturning model instance
AMOC = Psi_Thermwind(z=z,b1=b_Atl,b2=0.01*b_Atl(z))
# and solve for initial overturning streamfunction:
AMOC.solve()
# map North Atlantic overturning to isopycnal space:
[Psi_iso_Atl,Psi_iso_NAtl]=AMOC.Psibz()

# create N.P. overturning model instance
PMOC = Psi_Thermwind(z=z,b1=b_Pac,b2=0.01*b_Pac(z))
# and solve for initial overturning streamfunction:
PMOC.solve()
# map North Pacific overturning to isopycnal space:
[Psi_iso_Pac,Psi_iso_NPac]=PMOC.Psibz()


# create interbasin zonal overturning model instance
ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=8e-5)
#ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=2e-6)
# and solve for initial overturning streamfunction:
ZOC.solve()
# map inter-basin overturning to isopycnal space:
[Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()


# create S.O. overturning model instance for Atlantic sector
#SO_Atl=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=4.6e6,KGM=K.,Hsill=500.,HEk=100.,Htapertop=100.,Htaperbot=500.)
SO_Atl=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Lx/3.,KGM=K)
SO_Atl.solve()

# create S.O. overturning model instance for Pacific sector
#SO_Pac=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Lx*2./3.,KGM=K.,Hsill=500.,HEk=100.,Htapertop=100.,Htaperbot=500.)
SO_Pac=Psi_SO(z=z,y=y,b=b_Pac(z),bs=bs_SO,tau=tau,L=Lx*2./3.,KGM=K)
SO_Pac.solve()


# create adv-diff column model instance for Atl
Atl= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs,bbot=bbot,Area=A_Atl)
# create adv-diff column model instance for northern sinking region
northAtl= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs_northAtl,bbot=bbot,Area=A_northAtl)
northPac= Column(z=z,kappa=kappaeff,b=b_Pac,bs=bs_northPac,bbot=bbot,Area=A_northPac)
# create adv-diff column model instance for Pac
Pac= Column(z=z,kappa=kappaeff,b=b_Pac,bs=bs,bbot=bbot,Area=A_Pac)


# Create figure:
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))

# loop to iteratively find equilibrium solution
for ii in range(0, total_iters):    
   # update buoyancy profile
    # using isopycnal overturning:
   wA_Atl=(Psi_iso_Atl+Psi_zonal_Atl-SO_Atl.Psi)*1e6
   wANAtl=-Psi_iso_NAtl*1e6
   wA_Pac=(Psi_iso_Pac-Psi_zonal_Pac-SO_Pac.Psi)*1e6
   wANPac=-Psi_iso_NPac*1e6
   Atl.timestep(wA=wA_Atl,dt=dt)
   northAtl.timestep(wA=wANAtl,dt=dt,do_conv=True)
   Pac.timestep(wA=wA_Pac,dt=dt)
   northPac.timestep(wA=wANPac,dt=dt,do_conv=True)
   # with hor. advection: (not to be used with isopycnal overturning)
   #vdx=0*AMOC.Psi
   #vdx[1:-1]=-1e6*(AMOC.Psi[0:-2]-AMOC.Psi[2:])/(AMOC.z[0:-2]-AMOC.z[2:])
   #Atl.timestep(wA=wAb,dt=dt,vdx_in=-vdx,b_in=north.b)
   #north.timestep(wA=wAN,dt=dt,do_conv=True,vdx_in=vdx,b_in=Atl.b)
   
   if ii%MOC_up_iters==0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(b1=Atl.b,b2=northAtl.b)
      AMOC.solve()
      [Psi_iso_Atl, Psi_iso_NAtl]=AMOC.Psibz()
      PMOC.update(b1=Pac.b,b2=northPac.b)
      PMOC.solve()
      [Psi_iso_Pac, Psi_iso_NPac]=PMOC.Psibz()
      ZOC.update(b1=Atl.b,b2=Pac.b)
      ZOC.solve()
      [Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
      SO_Atl.update(b=Atl.b)
      SO_Atl.solve()
      SO_Pac.update(b=Pac.b)
      SO_Pac.solve()
     
   if ii%plot_iters==0:
      # Plot current state:
      ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=0.5)
      ax1.plot(PMOC.Psi, PMOC.z,'--m', linewidth=0.5)
      #ax1.plot(SO_Atl.Psi, SO_Atl.z, ':b', linewidth=0.5)
      #ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=0.5)
      ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=0.5)
      ax1.plot(SO_Atl.Psi-ZOC.Psi, SO_Atl.z, '--b', linewidth=0.5)
      ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=0.5)
      ax1.plot(SO_Atl.Psi+SO_Pac.Psi, z, '--c', linewidth=0.5)
      ax2.plot(Atl.b, Atl.z, '-b', linewidth=0.5)
      ax2.plot(Pac.b, Pac.z, '-g', linewidth=0.5)
      ax2.plot(northAtl.b, northAtl.z, '-r', linewidth=0.5)
      ax2.plot(northPac.b, northPac.z, '-m', linewidth=0.5)
      plt.pause(0.01)
 

ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
ax1.plot(PMOC.Psi, PMOC.z,'--m', linewidth=1.5)
#ax1.plot(SO_Atl.Psi, SO_Atl.z, ':b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
ax1.plot(SO_Atl.Psi-ZOC.Psi, SO_Atl.z, '--b', linewidth=1.5)
ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=1.5)
ax1.plot(SO_Atl.Psi+SO_Pac.Psi, z, '--c', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(northAtl.b, northAtl.z, '-r', linewidth=1.5)
ax2.plot(northPac.b, northPac.z, '-m', linewidth=1.5)

fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))
ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
ax1.plot(PMOC.Psi, PMOC.z,'--m', linewidth=1.5)
#ax1.plot(SO_Atl.Psi, SO_Atl.z, ':b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
ax1.plot(SO_Atl.Psi-ZOC.Psi, SO_Atl.z, '--b', linewidth=1.5)
ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=1.5)
ax1.plot(SO_Atl.Psi+SO_Pac.Psi, z, '--c', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(northAtl.b, northAtl.z, '-r', linewidth=1.5)
ax2.plot(northPac.b, northPac.z, '-m', linewidth=1.5)
#     
#
## interpolate results for fancy plots:
#barray=np.transpose(np.array([Atl.b,Atl.b,north.b]))
#barray_p=np.transpose(np.array([Pac.b,Pac.b]))
#barray_glob=np.transpose(np.array([(A_Pac*Pac.b+A_Atl*Atl.b)/(A_Pac+A_Atl),
#                                   (A_Pac*Pac.b+A_Atl*Atl.b)/(A_Pac+A_Atl),north.b]))
#barray_iso=np.transpose(np.array([Atl.b,Atl.b,north.b]))
#psiarray=np.transpose(np.array([SO_Atl.Psi-ZOC.Psi,AMOC.Psi,0.*AMOC.Psi]))
#psiarray_p=np.transpose(np.array([SO_Pac.Psi+ZOC.Psi,0*ZOC.Psi]))
#psiarray_glob=np.transpose(np.array([SO_Atl.Psi+SO_Pac.Psi,AMOC.Psi,0*AMOC.Psi]))
#psiarray_iso=np.transpose(np.array([SO_Atl.Psi-Psi_zonal_Atl,Psi_iso_Atl,Psi_iso_N]))
#yplot=np.array([0,10000,10700])
#yplot_p=np.array([0,10000])
#yplot_iso=np.array([0,9000,11000])
#blevs=np.arange(-2.0,3,0.1)
#plevs=np.arange(-21,23,2.0)
#
## plot results:
#fig = plt.figure(figsize=(10,4))
#ax1 = fig.add_subplot(111)
#CS=ax1.contour(yplot,SO_Atl.z,barray*100,levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#ax1.contour(yplot,SO_Atl.z,psiarray,levels=plevs,colors='k',linewidths=0.5)
##ax1.clabel(CS,  fontsize=10)
#CS=ax1.contourf(yplot,SO_Atl.z,psiarray,levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
#fig.colorbar(CS, ticks=plevs, orientation='vertical')
#ax1.set_xlabel('y [km]')
#ax1.set_ylabel('Depth [m]')
#
#fig = plt.figure(figsize=(10,4))
#ax1 = fig.add_subplot(111)
#CS=ax1.contour(yplot_p,SO_Pac.z,barray_p*100,levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#ax1.contour(yplot_p,SO_Pac.z,psiarray_p,levels=plevs,colors='k',linewidths=0.5)
##ax1.clabel(CS,  fontsize=10)
#CS=ax1.contourf(yplot_p,SO_Pac.z,psiarray_p,levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
#fig.colorbar(CS, ticks=plevs, orientation='vertical')
#ax1.set_xlabel('y [km]')
#ax1.set_ylabel('Depth [m]')
#
#fig = plt.figure(figsize=(10,4))
#ax1 = fig.add_subplot(111)
#CS=ax1.contour(yplot,SO_Atl.z,barray_glob*100,levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#ax1.contour(yplot,SO_Atl.z,psiarray_glob,levels=plevs,colors='k',linewidths=0.5)
##ax1.clabel(CS,  fontsize=10)
#CS=ax1.contourf(yplot,SO_Atl.z,psiarray_glob,levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
#fig.colorbar(CS, ticks=plevs, orientation='vertical')
#ax1.set_xlabel('y [km]')
#ax1.set_ylabel('Depth [m]')

       
#fig = plt.figure(figsize=(10,5))
#ax1 = fig.add_subplot(111)
#CS=ax1.contour(yplot_iso,SO_Atl.z,barray_iso*100,levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#ax1.contour(yplot_iso,SO_Atl.z,psiarray_iso,levels=plevs,colors='k',linewidths=0.5)
##ax1.clabel(CS,  fontsize=10)
#CS=ax1.contourf(yplot_iso,SO_Atl.z,psiarray_iso,levels=plevs,cmap=plt.cm.bwr, vmin=-16, vmax=16)
#fig.colorbar(CS, ticks=plevs, orientation='vertical')
#ax1.set_xlabel('y [km]')
#ax1.set_ylabel('Depth [m]')
       



#fig = plt.figure(figsize=(5,8))
#ax1 = fig.add_subplot(111)
#ax2 = ax1.twiny()
#ax1.plot(AMOC.Psi, AMOC.z, linewidth=0.5, color='r')
#ax1.plot(Psi_iso_atl, AMOC.z, linewidth=0.5, color='g')
#ax1.plot(Psi_iso_n, AMOC.z, linewidth=0.5, color='y')
#ax1.plot(SO.Psi, SO.z, linewidth=0.5, color='m')
#ax2.plot(Atl.b, Atl.z, linewidth=0.5,color='b')
#ax2.plot(north.b, north.z, linewidth=0.5,color='c')




