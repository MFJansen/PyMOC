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
sys.path.append('../../../Modules')
from psi_thermwind import Psi_Thermwind
from psi_SO import Psi_SO
from column import Column
import numpy as np
from matplotlib import pyplot as plt
from interp_channel import Interpolate_channel
from interp_twocol import Interpolate_twocol

diag_file=None

# boundary conditions:
bs=0.02;#0.015;
bs_north=0.0006; bbot= -0.0012

# S.O. surface boundary conditions and grid:
y=np.asarray(np.linspace(0,3.e6, 51))
tau=0.17
#tau=0.305*np.sin(np.pi*(y-1.11e5)/4.7e6)**2-0.068 
offset=0.0345*(1-np.cos(np.pi*(-1.e5)/8e6))+0.*1e-10;# Notice the 1e-10 is to
# make sure that bs_SO[0] is ever so slightly smaller than bbot so that bottom water comes in at teh bottom
bs_SO=(0.0345*(1-np.cos(np.pi*(y-1.e5)/8e6))
       +(bbot-offset)/5.55e5*np.maximum(0,5.55e5-y));

if bbot==0:
    bs_SO=np.maximum(bs_SO,0.)

GCM_b=np.array([-0.0015,-0.0004,0.0007,0.0013,0.0022,0.0032,0.0045,0.0060,
                0.0078,0.0099,0.0123,0.0144,0.0174,0.0212,0.0243,0.0266])

GCM_tau=np.array(
        [-0.0680,-0.0614,-0.0420,-0.0117,0.0271,0.0709,0.1159,0.1582,
          0.1942,0.2207,0.2354,0.2371,0.2255,0.2018,0.1679,0.1269])


GCM_y=1.11e5*np.array([1.,3.,5.,7.,9.,11.,13.,15.
                      ,17.,19.,21.,23.,25.,27.,29.,31.])

fig = plt.figure(figsize=(6,6))
plt.plot(GCM_y,GCM_b)
plt.plot(y,bs_SO)    
#
#fig = plt.figure(figsize=(6,6))
#plt.plot(GCM_y,GCM_tau)
#plt.plot(y,tau)    


#bs_SO=bbot+(bs-bbot)*y/(y[-1]-y[0])
# time-stepping parameters:
dt=86400.*30.                               # time-step for vert. adv. diff. calc.
MOC_up_iters=int(np.floor(1.*360*86400/dt)) # multiplier for MOC time-step (MOC is updated every MOC_up_iters time steps)
plot_iters= int(np.ceil(500*360*86400/dt))  # plotting frequency (in iterations)
total_iters=int(np.ceil(6000*360*86400/dt))# total number of timesteps

# Effective diffusivity profile
def kappaeff(z): # effective diffusivity profile with tapering in BBL
        #return (1e-5+3e-4*(np.exp((-4000-z)/1000.)))*(1.-np.maximum(-4000.-z+500.,0.)/500.)**2  
        #return (2e-5+1e-4*(z<-1500.))*(1.-np.maximum(-4000.-z+500.,0.)/500.)**2  
        return 1.0*( 1e-4*(1.1-np.tanh(np.maximum(z+2000.,0)/700. + np.minimum(z+2000.,0)/1000.))
                *(1.-np.maximum(-4000.-z+800.,0.)/800.)**2 ) 

        
#1.5 =10Sv overlap        
      
A_Atl=7e13#6e13  
A_north=5.5e12
A_Pac=1.7e14#1.44e14

Lx = 1.3e+07  #(length of the channel)
Latl=6./21.*Lx;
Lpac=15./21.*Lx;
K = 1700. # 1500 for variable wind stress % 1700 for const. wind stress

# create vertical grid:
z=np.asarray(np.linspace(-4000, 0, 80))

fig = plt.figure(figsize=(2.5,4))
plt.semilogx(kappaeff(z),z)
plt.xlim(1e-5,5e-4)
plt.xlabel('$\kappa_{eff}$ [m$^2$s$^{-1}$]',fontsize=12)
plt.ylabel('$z$ [m]',fontsize=13)
fig.tight_layout() 


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
ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=8.34e-5) # 8.34e-5 corresponds to 35S
#ZOC = Psi_Thermwind(z=z,b1=b_Atl,b2=b_Pac, f=2e-6)
# and solve for initial overturning streamfunction:
ZOC.solve()
# map inter-basin overturning to isopycnal space:
[Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()


# create S.O. overturning model instance for Atlantic sector
SO=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Latl+Lpac,KGM=K,Hsill=200.,HEk=200.,Htapertop=200.,Htaperbot=200.)
#SO=Psi_SO(z=z,y=y,b=b_Pac(z),bs=bs_SO,tau=tau,L=Latl+Lpac,KGM=K)
SO.solve()

# create S.O. overturning model instance for Pacific sector
#SO_Pac=Psi_SO(z=z,y=y,b=b_Atl(z),bs=bs_SO,tau=tau,L=Lx*2./3.,KGM=K.,Hsill=500.,HEk=100.,Htapertop=100.,Htaperbot=500.)
#SO_Pac=Psi_SO(z=z,y=y,b=b_Pac(z),bs=bs_SO,tau=tau,L=Lpac,KGM=K)
#SO_Pac.solve()


# create adv-diff column model instance for Atl
Atl= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs,bbot=bbot,Area=A_Atl)
# create adv-diff column model instance for northern sinking region
north= Column(z=z,kappa=kappaeff,b=b_Atl,bs=bs_north,bbot=bbot,Area=A_north,N2min=1e-7 )
# create adv-diff column model instance for Pac
Pac= Column(z=z,kappa=kappaeff,b=b_Pac,bs=bs,bbot=bbot,Area=A_Pac)


# Create figure:
fig = plt.figure(figsize=(6,9))
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))

# loop to iteratively find equilibrium solution
for ii in range(0, total_iters):    
   # update buoyancy profile
   #wAb=(AMOC.Psi-SO.Psi)*1e6
   #wAN=-AMOC.Psi*1e6
   # using isopycnal overturning: (Notice PsiSO_iso_atl is SO overturning interpolated to Atlantic buoyancy surfaces)
   PsiSO_iso_Atl=(np.interp(Atl.b,Pac.b,SO.Psi));
   # Notice that net Psi at southern end of Atl is PsiSO_batl-Psi_zonal_Atl
   wA_Atl=(Psi_iso_Atl-PsiSO_iso_Atl+Psi_zonal_Atl)*1e6
   wAN=-Psi_iso_N*1e6
   wA_Pac=(-Psi_zonal_Pac)*1e6
   Atl.timestep(wA=wA_Atl,dt=dt)
   north.timestep(wA=wAN,dt=dt,do_conv=True)
   Pac.timestep(wA=wA_Pac,dt=dt)
   # with hor. advection: (not to be used with isopycnal overturning)
   #vdx=0*AMOC.Psi
   #vdx[1:-1]=-1e6*(AMOC.Psi[0:-2]-AMOC.Psi[2:])/(AMOC.z[0:-2]-AMOC.z[2:])
   #Atl.timestep(wA=wAb,dt=dt,vdx_in=-vdx,b_in=north.b)
   #north.timestep(wA=wAN,dt=dt,do_conv=True,vdx_in=vdx,b_in=Atl.b)
   
   if ii%MOC_up_iters==0:
      # update overturning streamfunction (can be done less frequently)
      AMOC.update(b1=Atl.b,b2=north.b)
      AMOC.solve()
      [Psi_iso_Atl, Psi_iso_N]=AMOC.Psibz()
      ZOC.update(b1=Atl.b,b2=Pac.b)
      ZOC.solve()
      [Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
      SO.update(b=Pac.b)
      SO.solve()
      
   if ii%plot_iters==0:
      # Plot current state:
      ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=0.5)
      #ax1.plot(SO.Psi, SO.z, ':b', linewidth=0.5)
      #ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=0.5)
      ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=0.5)
      ax1.plot(SO.Psi-ZOC.Psi, SO.z, '--b', linewidth=0.5)
      ax1.plot(SO.Psi, z, '--c', linewidth=0.5)
      ax2.plot(Atl.b, Atl.z, '-b', linewidth=0.5)
      ax2.plot(Pac.b, Pac.z, '-g', linewidth=0.5)
      ax2.plot(north.b, north.z, '-r', linewidth=0.5)
      plt.pause(0.01)
 

ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
#ax1.plot(SO.Psi, SO.z, ':b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
ax1.plot(SO.Psi-ZOC.Psi, SO.z, '--b', linewidth=1.5)
ax1.plot(SO.Psi, z, '--c', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(north.b, north.z, '-r', linewidth=1.5)



if diag_file is not None:
   np.savez(diag_file, b_Atl=Atl.b, b_Pac=Pac.b, b_north=north.b,
            Psi_SO=SO.Psi, 
            Psi_ZOC=ZOC.Psi, Psib_ZOC=ZOC.Psib(), bgrid_ZOC=ZOC.bgrid,
            Psi_AMOC=AMOC.Psi, Psib_AMOC=AMOC.Psib(), bgrid_AMOC=AMOC.bgrid)


#*****************************************************************************
# Below is all just for fancy plots
#*****************************************************************************

#blevs=np.arange(-0.01,0.03,0.001) 
blevs=np.array([-0.004,-0.002,0.0,0.004,0.01,0.018])
#np.concatenate((np.arange(-0.01,0.006,0.002),np.arange(0.006,0.04,0.004))) 
#plevs=np.arange(-20.,22,2.0)
plevs=np.arange(-20.,22.,2.0)
nb=500;

b_basin=(A_Atl*Atl.b+A_Pac*Pac.b)/(A_Atl+A_Pac);
b_north=north.b;
# Compute total SO residual overturning by first interpolating circulation 
# in both sectors to mean buoyancy profile values and then summing (at constant b).
PsiSO=(np.interp(b_basin,Pac.b,SO.Psi));
AMOC_bbasin=np.interp(b_basin,AMOC.bgrid,AMOC.Psib(nb=nb));
# this is to Plot adiabatic overturning only:
PsiSO_bgrid=(np.interp(AMOC.bgrid,Pac.b,SO.Psi));
Psi_ad=np.maximum(np.minimum(AMOC.Psib(),PsiSO_bgrid),0.);
PsiSO_ad=np.interp(b_basin,AMOC.bgrid,Psi_ad);
#Psi_atl_ad=np.interp(Atl.b,AMOC.bgrid,Psi_ad);
Psi_north_ad=np.interp(north.b,AMOC.bgrid,Psi_ad);
AMOC_bbasin_ad=np.interp(b_basin,AMOC.bgrid,Psi_ad);

l=y[-1];

bs=1.*bs_SO;bn=1.*Pac.b;
bs[-1]=1.*bn[-1]; # if bs is tiny bit warmer there is a zero slope point at the surface which makes the interpolation fail...

# first interpolate buoyancy in channel along constant-slope isopycnals: 
bint=Interpolate_channel(y=y,z=z,bs=bs,bn=bn)
bsouth=bint.gridit()
# buoyancy in the basin is all the same:
lbasin=11000.
ltrans=1500.
lnorth=400.
lchannel=l/1e3
ybasin=np.linspace(0,lbasin,60)+lchannel;
bbasin=np.tile(b_basin,(len(ybasin),1))
bAtl=np.tile(Atl.b,(len(ybasin),1))
bPac=np.tile(Pac.b,(len(ybasin),1))
# interpolate buoyancy in northern ransition region:
ytrans=np.linspace(ltrans/20.,ltrans,20)+lchannel+lbasin;
bn=b_north.copy();
bn[0]=b_basin[0];# Notice that the interpolation procedure assumes that the bottom
#buoyancies in both colums match - which may not be exactly the case depending
# on when in teh time-step data is saved 
bint=Interpolate_twocol(y=ytrans*1000.-ytrans[0]*1000.,z=z,bs=Atl.b,bn=bn)
btrans=bint.gridit()
# finally set buyancy in northern deep water formation region:
ynorth=np.linspace(lnorth/20.,lnorth,20)+lchannel+lbasin+ltrans;
bnorth=np.tile(bn,(len(ynorth),1))
# now stick it all together:
ynew=np.concatenate((y/1e3,ybasin,ytrans,ynorth))
bnew=np.concatenate((bsouth,bbasin,btrans,bnorth))
bnew_Atl=np.concatenate((np.nan*bsouth,bAtl,btrans,bnorth))
bnew_Pac=np.concatenate((np.nan*bsouth,bPac,np.nan*btrans,np.nan*bnorth))

# Compute z-coordinate and residual overturning streamfunction at all latitudes:
psiarray_b=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_b_Atl=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_b_Pac=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_b_ad=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_b_Atl_ad=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_b_Pac_ad=np.zeros((len(ynew),len(z))) # overturning in b-coordinates
psiarray_res=np.zeros((len(ynew),len(z))) # "residual" overturning - i.e. isopycnal overturning mapped into z space
psiarray_Atl=np.zeros((len(ynew),len(z))) # "residual" overturning in ATl.
psiarray_Pac=np.zeros((len(ynew),len(z))) # "residual" overturning in Pac.
psiarray_z=np.zeros((len(ynew),len(z)))  # z-space, "eulerian" overturning
psiarray_z_Atl=np.zeros((len(ynew),len(z)))  # z-space, "eulerian" overturning
psiarray_z_Pac=np.zeros((len(ynew),len(z)))  # z-space, "eulerian" overturning
for iy in range(1,len(y)):
    # Psi in SO channel is just total SO overturning (same field is used to fill Atl and Pac fields)
    psiarray_res[iy,:]=np.interp(bnew[iy,:],Pac.b,SO.Psi)
    psiarray_z[iy,:]=np.interp(bnew[iy,:],Pac.b,SO.Psi)
    psiarray_z_Atl[iy,:]=psiarray_z[iy,:]
    psiarray_z_Pac[iy,:]=psiarray_z[iy,:]
    psiarray_b[iy,b_basin<bs_SO[iy]]=PsiSO[b_basin<bs_SO[iy]]
    psiarray_b_ad[iy,b_basin<bs_SO[iy]]=PsiSO_ad[b_basin<bs_SO[iy]]
    psiarray_Atl[iy,:]=psiarray_res[iy,:];    
    psiarray_b_Atl[iy,:]=psiarray_b[iy,:];    
    psiarray_b_Atl_ad[iy,:]=psiarray_b_ad[iy,:];    
    psiarray_Pac[iy,:]=psiarray_res[iy,:];    
    psiarray_b_Pac[iy,:]=psiarray_b[iy,:];    
    psiarray_b_Pac_ad[iy,:]=psiarray_b_ad[iy,:];    
for iy in range(len(y),len(y)+len(ybasin)):
    # in the basin, linearly interpolate between Psi at southern and northern end of respective basin (or global ocean):
    psiarray_res[iy,:]=((ynew[iy]-lchannel)*np.interp(b_basin,AMOC.bgrid,AMOC.Psib(nb=nb))+(lchannel+lbasin-ynew[iy])*PsiSO)/lbasin   
    psiarray_z[iy,:]=((ynew[iy]-lchannel)*AMOC.Psi+(lchannel+lbasin-ynew[iy])*(SO.Psi))/lbasin  
    psiarray_z_Atl[iy,:]=((ynew[iy]-lchannel)*AMOC.Psi+(lchannel+lbasin-ynew[iy])*(SO.Psi-ZOC.Psi))/lbasin  
    psiarray_z_Pac[iy,:]=((lchannel+lbasin-ynew[iy])*(ZOC.Psi))/lbasin  
    psiarray_b[iy,:]=((ynew[iy]-lchannel)*AMOC_bbasin
                      +(lchannel+lbasin-ynew[iy])*PsiSO)/lbasin  
    psiarray_b_ad[iy,:]=((ynew[iy]-lchannel)*AMOC_bbasin_ad
                      +(lchannel+lbasin-ynew[iy])*PsiSO_ad)/lbasin  
    psiarray_Atl[iy,:]=((ynew[iy]-lchannel)*AMOC.Psibz(nb=nb)[0]+(lchannel+lbasin-ynew[iy])*(SO.Psi-ZOC.Psibz()[0]))/lbasin   
    psiarray_b_Atl[iy,:]=((ynew[iy]-lchannel)*AMOC_bbasin+
                    (lchannel+lbasin-ynew[iy])*(np.interp(b_basin,Pac.b,SO.Psi)
                      -np.interp(b_basin,ZOC.bgrid,ZOC.Psib())))/lbasin   
    psiarray_b_Atl_ad[iy,:]=psiarray_b_ad[iy,:]
    psiarray_Pac[iy,:]=((lchannel+lbasin-ynew[iy])*(ZOC.Psibz()[1]))/lbasin   
    psiarray_b_Pac[iy,:]=((lchannel+lbasin-ynew[iy])*(
                        np.interp(b_basin,ZOC.bgrid,ZOC.Psib())))/lbasin   
    psiarray_b_Pac_ad[iy,:]=0.*psiarray_b_ad[iy,:]
for iy in range(len(y)+len(ybasin),len(y)+len(ybasin)+len(ytrans)):
    # in the northern transition region, interpolate AMOC.psib to local isopycnal depth
    # for psi_res, while keeping psi_z constant and psi_z constant on non-outcropped isopycnals:
    psiarray_res[iy,:]=np.interp(bnew[iy,:],AMOC.bgrid,AMOC.Psib(nb=nb))
    psiarray_z[iy,:]=AMOC.Psi      
    psiarray_z_Atl[iy,:]=AMOC.Psi      
    psiarray_z_Pac[iy,:]=np.NaN      
    #psiarray_z[iy,:]=((lchannel+lbasin+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_b[iy,b_basin<bnew[iy,-1]]=AMOC_bbasin[b_basin<bnew[iy,-1]]      
    psiarray_b_ad[iy,b_basin<bnew[iy,-1]]=AMOC_bbasin_ad[b_basin<bnew[iy,-1]]      
    psiarray_Atl[iy,:]=np.interp(bnew[iy,:],AMOC.bgrid,AMOC.Psib(nb=nb))
    psiarray_b_Atl[iy,b_basin<bnew[iy,-1]]=psiarray_b[iy,b_basin<bnew[iy,-1]]
    psiarray_b_Atl_ad[iy,:]=psiarray_b_ad[iy,:]      
    psiarray_Pac[iy,:]=np.NaN
    psiarray_b_Pac[iy,:]=np.NaN
    psiarray_b_Pac_ad[iy,:]=np.NaN
for iy in range(len(y)+len(ybasin)+len(ytrans),len(ynew)):
    # in the northern sinking region, all psi decrease linearly to zero:
    psiarray_res[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psibz(nb=nb)[1])/lnorth
    psiarray_z[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_z_Atl[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psi)/lnorth      
    psiarray_z_Pac[iy,:]=np.NaN      
    psiarray_b[iy,b_basin<bnew[iy,-1]]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])
        *AMOC_bbasin[b_basin<bnew[iy,-1]])/lnorth      
    psiarray_b_ad[iy,b_basin<bnew[iy,-1]]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])
        *AMOC_bbasin_ad[b_basin<bnew[iy,-1]])/lnorth      
    psiarray_Atl[iy,:]=((lchannel+lbasin+ltrans+lnorth-ynew[iy])*AMOC.Psibz(nb=nb)[1])/lnorth
    psiarray_b_Atl[iy,b_basin<bnew[iy,-1]]=psiarray_b[iy,b_basin<bnew[iy,-1]]
    psiarray_b_Atl_ad[iy,:]=psiarray_b_ad[iy,:]
    psiarray_Pac[iy,:]=np.NaN
    psiarray_b_Pac[iy,:]=np.NaN
#psiarray_res[-1,:]=0.;

# plot z-coordinate overturning and buoyancy structure:
fig = plt.figure(figsize=(10.8,6.8))
ax1 =fig.add_axes([0.1,.57,.33,.36])
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_ylabel('Depth [m]', fontsize=13)
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))
ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
ax1.plot(SO.Psi-ZOC.Psi, SO.z, '--b', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, '--g', linewidth=1.5)
ax1.plot(SO.Psi, z, '--k', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(0.*z, z,linewidth=0.5,color='k',linestyle=':')

ax1 = fig.add_subplot(222)
CS=ax1.contour(ynew,z,bnew.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_z.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_z.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Global',fontsize=12)
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
       
ax1 = fig.add_subplot(223)
CS=ax1.contour(ynew,z,bnew_Pac.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_z_Pac.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_z_Pac.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
plt.gca().invert_xaxis()
#ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Pacific',fontsize=12)

ax1 = fig.add_subplot(224)
CS=ax1.contour(ynew,z,bnew_Atl.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_z_Atl.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_z_Atl.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Atlantic',fontsize=12)
cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
fig.tight_layout()   
fig.subplots_adjust(right=0.92)
fig.colorbar(CS, cax=cbar_ax,ticks=plevs[0::5], orientation='vertical')
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')




# plot global residual overturning and global mean buoyancy structure:
        
#fig = plt.figure(figsize=(10,7.5))
fig = plt.figure(figsize=(10.8,6.8))


# Notice that this isn't all really "residual", but it's also not clear
# which buoyancy to interpolate residual overturning between two columns to
ax1 =fig.add_axes([0.1,.57,.33,.36])
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_ylabel('Depth [m]', fontsize=13)
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))
ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
ax1.plot(AMOC.Psi, AMOC.z,'--r', linewidth=1.5)
#ax1.plot(SO.Psi, SO.z, ':b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
ax1.plot(SO.Psi-ZOC.Psibz()[0], SO.z, '--b', linewidth=1.5)
ax1.plot(ZOC.Psibz()[1], ZOC.z, '--g', linewidth=1.5)
ax1.plot(PsiSO, z, '--k', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(0.*z, z,linewidth=0.5,color='k',linestyle=':')

ax1 = fig.add_subplot(222)
CS=ax1.contour(ynew,z,bnew.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_res.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_res.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Global',fontsize=12)
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
       
ax1 = fig.add_subplot(223)
CS=ax1.contour(ynew,z,bnew_Pac.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_Pac.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_Pac.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Pacific',fontsize=12)
plt.gca().invert_xaxis()

ax1 = fig.add_subplot(224)
CS=ax1.contour(ynew,z,bnew_Atl.transpose()/2e-3,levels=blevs/2e-3,colors='k',linewidths=1.0,linestyles='solid')
ax1.clabel(CS,fontsize=10,fmt='%1.0f')
ax1.contour(ynew,z,psiarray_Atl.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_Atl.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('Depth [m]',fontsize=12)
ax1.set_title('Atlantic',fontsize=12)
cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
fig.tight_layout()   
fig.subplots_adjust(right=0.92)
fig.colorbar(CS, cax=cbar_ax,ticks=plevs[0::5], orientation='vertical')
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')


# Plot Isopycnal overturning

fig = plt.figure(figsize=(10.8,6.8))
#ax1 = fig.add_subplot(421, colspan = 3)
#ax1 =fig.add_axes([0.1,.57,.31,.36])
ax1 =fig.add_axes([0.1,.57,.33,.36])
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_ylabel('b [m s$^{-2}$]', fontsize=13)
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))
ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
ax1.plot(AMOC_bbasin, z,'--r', linewidth=1.5)
ax1.plot(np.interp(b_basin,ZOC.bgrid,ZOC.Psib()), z, ':c', linewidth=1.5)
ax1.plot(np.interp(b_basin,Pac.b,SO.Psi)
        -np.interp(b_basin,ZOC.bgrid,ZOC.Psib()), z, '--b', linewidth=1.5)
ax1.plot(np.interp(b_basin,ZOC.bgrid,ZOC.Psib()), z, '--g', linewidth=1.5)
ax1.plot(PsiSO, z, '--k', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(0.*b_basin, b_basin,linewidth=0.5,color='k',linestyle=':')
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])

ax1 = fig.add_subplot(222)
ax1.contour(ynew,z,psiarray_b.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Global',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
       
ax1 = fig.add_subplot(223)
ax1.contour(ynew,z,psiarray_b_Pac.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Pac.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Pacific',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
plt.gca().invert_xaxis()
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')

ax1 = fig.add_subplot(224)
ax1.contour(ynew,z,psiarray_b_Atl.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Atl.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Atlantic',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
fig.tight_layout()   
fig.subplots_adjust(right=0.92)
fig.colorbar(CS, cax=cbar_ax,ticks=plevs[0::5], orientation='vertical')






# Plot adiabatic Isopycnal overturning

fig = plt.figure(figsize=(10.8,6.8))
##ax1 = fig.add_subplot(421, colspan = 3)
##ax1 =fig.add_axes([0.1,.57,.31,.36])
#ax1 =fig.add_axes([0.1,.57,.33,.36])
#ax2 = ax1.twiny()
#plt.ylim((-4e3,0))
#ax1.set_ylabel('b [m s$^{-2}$]', fontsize=13)
#ax1.set_xlim((-20,30))
#ax2.set_xlim((-0.02,0.030))
#ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
#ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
#ax1.plot(PsiSO_ad, z,'--r', linewidth=1.5)
##ax1.plot(np.interp(b_basin,ZOC.bgrid,0.*ZOC.Psib()), z, ':c', linewidth=1.5)
##ax1.plot(np.interp(b_basin,Atl.b,SO.Psi), z, '--b', linewidth=1.5)
##ax1.plot(np.interp(b_basin,Pac.b,SO_Pac.Psi)
##          +np.interp(b_basin,ZOC.bgrid,ZOC.Psib()), z, '--g', linewidth=1.5)
#ax1.plot(PsiSO_ad, z, '--k', linewidth=1.5)
#ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
#ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
#ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
#h1, l1 = ax1.get_legend_handles_labels()
#h2, l2 = ax2.get_legend_handles_labels()
#ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
#ax1.plot(0.*b_basin, b_basin,linewidth=0.5,color='k',linestyle=':')
#ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
#ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])

ax1 = fig.add_subplot(222)
ax1.contour(ynew,z,psiarray_b.transpose()-psiarray_b_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b.transpose()-psiarray_b_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Global',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
       
ax1 = fig.add_subplot(223)
ax1.contour(ynew,z,psiarray_b_Pac.transpose()-psiarray_b_Pac_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Pac.transpose()-psiarray_b_Pac_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Pacific',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])

#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')

ax1 = fig.add_subplot(224)
ax1.contour(ynew,z,psiarray_b_Atl.transpose()-psiarray_b_Atl_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Atl.transpose()-psiarray_b_Atl_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Atlantic',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
plt.gca().invert_xaxis()
cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
fig.tight_layout()   
fig.subplots_adjust(right=0.92)
fig.colorbar(CS, cax=cbar_ax,ticks=plevs[0::5], orientation='vertical')


# Plot adiabatic Diapycnal overturning

fig = plt.figure(figsize=(10.8,6.8))
#ax1 = fig.add_subplot(421, colspan = 3)
#ax1 =fig.add_axes([0.1,.57,.31,.36])
ax1 =fig.add_axes([0.1,.57,.33,.36])
ax2 = ax1.twiny()
plt.ylim((-4e3,0))
ax1.set_ylabel('b [m s$^{-2}$]', fontsize=13)
ax1.set_xlim((-20,30))
ax2.set_xlim((-0.02,0.030))
ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
ax1.plot(PsiSO_ad, z,'--r', linewidth=1.5)
#ax1.plot(np.interp(b_basin,ZOC.bgrid,0.*ZOC.Psib()), z, ':c', linewidth=1.5)
#ax1.plot(np.interp(b_basin,Atl.b,SO.Psi), z, '--b', linewidth=1.5)
#ax1.plot(np.interp(b_basin,Pac.b,SO_Pac.Psi)
#          +np.interp(b_basin,ZOC.bgrid,ZOC.Psib()), z, '--g', linewidth=1.5)
ax1.plot(PsiSO_ad, z, '--k', linewidth=1.5)
ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
ax1.plot(0.*b_basin, b_basin,linewidth=0.5,color='k',linestyle=':')
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])

ax1 = fig.add_subplot(222)
ax1.contour(ynew,z,psiarray_b_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Global',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
       
ax1 = fig.add_subplot(223)
ax1.contour(ynew,z,psiarray_b_Pac_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Pac_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Pacific',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
plt.gca().invert_xaxis()
cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
fig.tight_layout()   
fig.subplots_adjust(right=0.92)
fig.colorbar(CS, cax=cbar_ax,ticks=plevs[0::5], orientation='vertical')

ax1 = fig.add_subplot(224)
ax1.contour(ynew,z,psiarray_b_Atl_ad.transpose(),levels=plevs,colors='k',linewidths=0.5)
CS=ax1.contourf(ynew,z,psiarray_b_Atl_ad.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-20, vmax=20)
ax1.plot(np.array([y[-1]/1000.,y[-1]/1000.]),np.array([-4000.,0.]),color='k')
ax1.set_xlim([0,ynew[-1]])
ax1.set_xlabel('y [km]',fontsize=12)
ax1.set_ylabel('b [m s$^{-2}$]',fontsize=12)
ax1.set_title('Atlantic',fontsize=12)
ax1.set_yticks(np.interp([0.02, 0.005,0.002,0.001,0.0005, 0.,-0.0005, -0.001 ],b_basin,z))
ax1.set_yticklabels([0.02, 0.005,0.002, 0.001,0.0005, 0.,-0.0005, -0.001 ])
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')


#fig = plt.figure(figsize=(9,6))
#ax1 =fig.add_subplot(121)
#plt.ylim((-0.005,0.015))
#ax1.set_ylabel('b [m s$^{-2}$]', fontsize=13)
#ax1.set_xlim((-30,15))
#ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
#ax1.plot(-AMOC.Psib(), AMOC.bgrid,'-b', linewidth=1.5)
#ax1.plot(SO.Psi, SO.b(SO.z),'--b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.b(SO_Pac.z),'--g', linewidth=1.5)
#[Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
#ax1.plot((AMOC.Psibz()[0]+Psi_zonal_Atl-SO.Psi), Atl.b,':b', linewidth=1.5)
#ax1.plot((-Psi_zonal_Pac-SO_Pac.Psi), Atl.b,':g', linewidth=1.5)
#ax1.plot(0.*AMOC.bgrid, AMOC.bgrid,linewidth=0.5,color='k',linestyle='-')
#ax1 =fig.add_subplot(122)
#plt.ylim((-0.005,0.005))
#ax1.set_ylabel('b [m s$^{-2}$]', fontsize=13)
#ax1.set_xlim((-25,15))
#ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
#ax1.plot(Psi_zonal_Atl-SO.Psi,Atl.b,'-b', linewidth=1.5)
#ax1.plot(SO.Psi, SO.b(SO.z),'--b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.b(SO_Pac.z),'--g', linewidth=1.5)
#[Psi_zonal_Atl,Psi_zonal_Pac]=ZOC.Psibz()
#ax1.plot((-Psi_zonal_Pac-SO_Pac.Psi), Atl.b,':g', linewidth=1.5)
#ax1.plot(0.*AMOC.bgrid, AMOC.bgrid,linewidth=0.5,color='k',linestyle='-')

#ax1.plot(SO.Psi, SO.z, ':b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi, SO_Pac.z, ':g', linewidth=1.5)
#ax1.plot(ZOC.Psi, ZOC.z, ':c', linewidth=1.5)
#ax1.plot(SO.Psi-ZOC.Psi, SO.z, '--b', linewidth=1.5)
#ax1.plot(SO_Pac.Psi+ZOC.Psi, SO_Pac.z, '--g', linewidth=1.5)
#ax1.plot(SO.Psi+SO_Pac.Psi, z, '--c', linewidth=1.5)
#ax2.plot(Atl.b, Atl.z, '-b', linewidth=1.5)
#ax2.plot(Pac.b, Pac.z, '-g', linewidth=1.5)
#ax2.plot(north.b, north.z, '-r', linewidth=1.5)     
#h1, l1 = ax1.get_legend_handles_labels()
#h2, l2 = ax2.get_legend_handles_labels()
#ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)





## plot z-coord. overturning:
#fig = plt.figure(figsize=(6.5,3))
#ax1 = fig.add_subplot(111)
#ax1.plot(np.array([l/1e3,l/1e3]),np.array([z[0],z[-1]]),color='0.5',linewidth=0.7,linestyle='dashed')
#CS=ax1.contour(ynew,z,bnew.transpose(),levels=blevs,colors='k',linewidths=1.0,linestyles='solid')
#ax1.clabel(CS,fontsize=10)
#CS=ax1.contourf(ynew,z,psiarray_z.transpose(),levels=plevs,cmap=plt.cm.bwr, vmin=-17, vmax=17)
#ax1.contour(ynew,z,psiarray_z.transpose(),levels=plevs,colors='0.5',linewidths=0.5)
#ax1.set_xlim([0,ynew[-1]])
#ax1.set_xlabel('y [km]',fontsize=12)
#ax1.set_ylabel('Depth [m]',fontsize=12)
#ax1.set_title('Depth-averaged Overturning',fontsize=12)
#fig.colorbar(CS, ticks=plevs[0::5], orientation='vertical')
#fig.tight_layout()   
##fig.savefig('psi_b_2D_depth_adiabatic_revised.png', format='png', dpi=600)
#

## Plot profiles:
#fig = plt.figure(figsize=(3.5,3.5))
#ax1 = fig.add_subplot(111)
#ax2 = ax1.twiny()
#plt.ylim((-4e3,0))
#ax1.set_ylabel('Depth [m]', fontsize=13)
#ax1.set_xlim((-8,14))
#ax2.set_xlim((-0.016,0.028))
#ax1.set_xlabel('$\Psi$ [SV]', fontsize=13)
#ax2.set_xlabel('$b_B$ [m s$^{-2}$]', fontsize=13)
#ax1.plot(AMOC.Psi, AMOC.z,linewidth=2,color='m',linestyle='--',label='$\Psi_N$')
#ax1.plot(PsiSO, PsiSO.z,linewidth=2,color='c',linestyle='--',label='$\Psi_{SO}$')
#ax2.plot(b_north, z, linewidth=2,color='r',label='$b_N$')
#ax2.plot(b_basin, z, linewidth=2,color='b',label='$b_B$')
#h1, l1 = ax1.get_legend_handles_labels()
#h2, l2 = ax2.get_legend_handles_labels()
#ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
#ax1.plot(0.*z, z,linewidth=0.5,color='k',linestyle=':')
#fig.tight_layout()
##fig.savefig('profiles.png', format='png', dpi=600)
#
#
## Plot SO ML results:
#fig = plt.figure(figsize=(4,3.1))
#ax1 = fig.add_subplot(111)
#ax2 = ax1.twinx()
#plt.xlim((0,2600))
#ax1.set_xlabel('y [km]', fontsize=13)
#ax1.set_ylim((-5.,10.))
#ax2.set_ylim((-0.01,0.02))
#ax1.set_ylabel('$\Psi_{SO}$ [SV]', fontsize=13)
#ax2.set_ylabel('$b_{SO}$ [m s$^{-2}$]', fontsize=13)
##ax1.plot(y/1000.,channel.Psi_s,linewidth=2,color='c',linestyle='--',label='$\Psi_{SO}$')
#ax2.plot(y/1000.,bs_SO, linewidth=2,color='r',label='$b_{SO}$')
#h1, l1 = ax1.get_legend_handles_labels()
#h2, l2 = ax2.get_legend_handles_labels()
#ax1.legend(h1+h2, l1+l2, loc=4, frameon=False)
#ax1.plot(y/1000.,0.*y,linewidth=0.5,color='k',linestyle=':')
#fig.tight_layout()   
##fig.savefig('bs_Psi_SO.png', format='png', dpi=600)



