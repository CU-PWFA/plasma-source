#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 10:06:04 2019

@author: keenan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:59:42 2019

@author: mike
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 16:45:17 2019

@author: mike
"""


"""
This code is used to analyze the butterfly pattern generated at
the spectrometer screen for various matching scenarios in
the PWFA. The beam starts at the beginning of the flat-top
with alpha=0 and beta determined by the desired B-mag value.
The exit ramp is shaped to produce a virtual vacuum waist with
the value beta_x0 for twice the initial beam energy. The code
assumes that the imaging spectrometer is automatically tuned
to image the centroid energy from the waist location obtained
by analyzing the full, projected beam Twiss parameters at QS0.
"""


# standard python code
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import patches
import scipy.constants as const
from scipy import optimize
import psutil
import timeit
from time import time
tic = timeit.default_timer() # time the process

# add WARGSim directory location to $PATH
sys.path.insert(0, "/home/keenan/WARGSim")

# import WARGSim packages
from beams import electronbeam
from beams import betatronbeam
from elements import pwfa as pwfa
from elements import drift as drift
from elements import quadrupole as quadrupole
from elements import dipole as dipole
from elements import beamline as beamline
#from optics import profile
from interactions import interactions

# Define output location. Don't forget final slash!
#path = '/home/mike/work/github_repos/CU-PWFA/plasma-source/litos/butterfly/'
path = '/media/keenan/Data_Storage/WARGSim/Dump';

# Get the number of available cores
Ncores = psutil.cpu_count()
numThread  = Ncores # int, number of threads to run on



Nslice = 1000


#%%
# define electron beam

N       = int(1e6) # number of particles in the beam
# eV, electron rest mass energy:
me = const.physical_constants['electron mass energy equivalent in MeV'][0]*(1e6)
beamE   = 10e9 # beam energy in eV
gb0     = beamE/me # reltivistic Lorentz factor (centroid value)
eps_n0   = 3.0e-6 # m-rad, normalized emittance
beta_x0  = 0.25   # m, x beta function at z=0
beta_y0  = 0.25   # m, y beta function at z=0
alpha_x0 = 0.00   # x alpha function at z=0
alpha_y0 = 0.00   # y alpha function at z=0
#gamma_x0 = (1+alpha_x0**2)/beta_x0 # 1/m, x gamma function at z=0
#gamma_y0 = (1+alpha_y0**2)/beta_y0 # 1/m, y gamma function at z=0
rms_z0  = 0.00   # m, rms bunch length
rms_gb0 = 0.01 # relative energy spread


# Useful functions for checking emittance growth

def match_beta(np0,gb0):
    np0 = np0*(1e17)
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = 1.0/kbeta
    return beta

def calc_Bmag(beta,beta_m):
    return (1/2)*(beta/beta_m + beta_m/beta)

def calc_beta_fact(Bmag):
    return (Bmag+np.sqrt((Bmag**2) - 1))




# force matched Twiss params:
Bmag_target = 5.00



n0 = 0.34 # make sure this agrees with value used to define plasma!!
beta_m = match_beta(n0,gb0)
beta_x00 = calc_beta_fact(Bmag_target)*beta_m
beta_y00 = calc_beta_fact(Bmag_target)*beta_m
Bmag = calc_Bmag(beta_x0,beta_m)
print(f'beta factor = {beta_x0/beta_m :.2f}')
print(f'Bmag = {Bmag :.2f}')


electronParams = {
        'name'      : 'ElectronBeam',
        'path'      : path,
        'load'      : False,
        'N'         : N,
        'shape'     : 'Gauss',
        'gb0'       : gb0,
        'rms_gb'    : rms_gb0,
        'rms_z'     : rms_z0,
        'eps_nx'    : eps_n0,
        'eps_ny'    : eps_n0,
        'beta_x'    : beta_x00,
        'beta_y'    : beta_y00,
        'alpha_x'   : alpha_x0,
        'alpha_y'   : alpha_y0,
        'seed'      : 3680995680
    }
my_ebeam = electronbeam.ElectronBeam(electronParams)
print("Beam random seed =", my_ebeam.seed)
eps_nx0, eps_ny0 = my_ebeam.get_emit_n()
print('init eps_nx = ',eps_nx0)
print('init eps_ny = ',eps_ny0)
init_ptcls = np.copy(my_ebeam.ptcls);

#%%
# introduce offset at vacuum waist

sig_x = np.sqrt(eps_n0*beta_x0/gb0)
sig_xp = np.sqrt(eps_n0/(beta_x0*gb0))
sig_y = np.sqrt(eps_n0*beta_y0/gb0)
sig_yp = np.sqrt(eps_n0/(beta_y0*gb0))

#my_ebeam.ptcls[:,0] += 1.00*sig_x
#my_ebeam.ptcls[:,1] += 1.00*sig_xp
#my_ebeam.ptcls[:,2] += 1.00*sig_y
#my_ebeam.ptcls[:,3] += -1.00*sig_yp


#%%
# define the betatron radiation beam

betatronParams = {
        'name'    : 'BetatronBeam',
        'Nxp'      : 100,
        'Nyp'      : 100,
        'NEc'      : 100,
        'xp_range' : 0.0020,
        'yp_range' : 0.0020,
        'Ec_range' : 10000
    }
my_bbeam = betatronbeam.BetatronBeam(betatronParams)


#%%
# define the plasma

def match_hw(np0,gb0,beta0):
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = beta0*kbeta
    hw = (1.96e-3)*beta**2 + 1.49*beta - 3.08
    return hw/kbeta

def match_zw(np0,gb0,beta0):
    kbeta = (1.33e-4)*np.sqrt(np0/gb0)
    beta = beta0*kbeta
    zw = (-1.47e-2)*beta**2 - 4.34*beta + 17.9
    return zw/kbeta

#n0    = 0.34  # 1/(10^17 cm^-3), flat-top plasma density
shape = 'Gauss' # ramp shape
L_ft  = 0.40 # m, length of flat-top
hw_up = 0 # m, half-width of up-ramp
L_up  = 8*hw_up # m, full length of up-ramp
hw_dn = 1.00*match_hw(n0*(1e17),2*gb0,beta_x0) # m, half-width of down-ramp
L_dn  = min(1.600,8*hw_dn) #L_up  # m, full length of down-ramp
L_p   = L_up + L_ft + L_dn # m, full length of plasma region

pwfaParams ={
    'name'  : 'PWFA0',
    'L'     : L_p,
    'n0'    : n0,
    'autoNz': True,
    'gb0'   : gb0,
    'shape' : shape,
    'L_ft'  : L_ft,
    'hw_up' : hw_up,
    'hw_dn' : hw_dn
}
pwfa0 = pwfa.PWFA(pwfaParams)
#pwfa0.dz = np.array([pwfa0.dz[0], pwfa0.dz[0]])

#%%
dumpPeriod = 1e9
 # int, step period for electron beam dump
# propagate the beam through the plasma

# propagate the beam from the waist back to the start of the simulation
#waist      = match_zw(n0*(1e17),gb0,beta_x0) # -0.3884 # m, waist distance upstream of plasma flat-top start
#Z_back     = -(L_up+waist)
#interactions.electron_vacuum(my_ebeam, np.linspace(0,Z_back,2), dumpPeriod, numThread)

# propagate the beam through the PWFA
#interactions.ebeam_bbeam_pwfa(my_ebeam, my_bbeam, pwfa0, dumpPeriod, numThread)
interactions.ebeam_pwfa(my_ebeam, pwfa0, dumpPeriod, numThread)
#interactions.ebeam_pwfa2(my_ebeam, pwfa0, dumpPeriod, numThread)
#interactions.ebeam_pwfa_slice(my_ebeam, pwfa0, Nslice, numThread)

#
"""
#%%
# plot total radiated power

int_rad = np.zeros([np.shape(my_bbeam.rad)[0],np.shape(my_bbeam.rad)[1]])
for i in range(np.shape(my_bbeam.rad)[2]):
    int_rad += my_bbeam.rad[:,:,i]*i

int_rad = int_rad[1:-1,1:-1]

fig1, ax1 = plt.subplots(figsize=(6, 6))

cmap0 = plt.get_cmap('PuBu_r')

plt.imshow(int_rad,cmap=cmap0)

plt.colorbar()
#plt.clim(0,5e7)

print("sum: %1.3e" % sum(sum(int_rad)))
print("peak: %1.3e" % np.max(int_rad))


"""
#%%
# get the virtual vacuum waist info

my_CS = my_ebeam.get_CS()

gamma_wx = my_CS['gamma_x']
gamma_wy = my_CS['gamma_y']

beta_wx = 1.0/gamma_wx
beta_wy = 1.0/gamma_wy

print('virt. waist beta_x (cm): ',beta_wx/0.01)
print('virt. waist beta_y (cm): ',beta_wy/0.01)

beta_x = my_CS['beta_x']
beta_y = my_CS['beta_y']

z_wx = np.sqrt((beta_x-beta_wx)*beta_wx)
z_wy = np.sqrt((beta_y-beta_wy)*beta_wy)

print('virt. waist z_x (m): ',z_wx)
print('virt. waist z_y (m): ',z_wy)

gb0 = my_ebeam.get_gb0()
beamE  = gb0*me

print('post-plasma beam energy (GeV): ',beamE/(1e9))

eps_nx, eps_ny = my_ebeam.get_emit_n()
print('post-plasma eps_nx = ',eps_nx)
print('post-plasma eps_ny = ',eps_ny)

print(f'emit. ratio in x = {eps_nx/eps_nx0 :.2f}')
print(f'emit. ratio in y = {eps_nx/eps_nx0 :.2f}')




#%%
# propagate backward to the virtual vacuum waist

interactions.electron_vacuum(my_ebeam, np.linspace(0,-z_wx,2), dumpPeriod, numThread)


#%%
# plot transverse phase space at plasma exit virtual waist

def calc_Twiss(x,xp,weights=None):
    cen_x  = np.average(x, weights=weights)
    dx     = x - cen_x
    cen_xp = np.average(xp, weights=weights)
    dxp    = xp - cen_xp
    # Calculate the RMS sizes and the correlation
    sigmax2  = np.average(dx**2,  weights=weights)
    sigmaxp2 = np.average(dxp**2, weights=weights)
    sigmaxxp = np.average(dx*dxp, weights=weights)
    # Calculate the emittance
    ex = np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
    # Assign the outputs
    eps   = ex
    beta  = sigmax2/ex
    alpha = sigmaxxp/ex
    gamma = sigmaxp2/ex
    return eps, beta, alpha, gamma


fig1, ax1 = plt.subplots(figsize=(6, 6))

x  = my_ebeam.ptcls[:,0]*1e6
xp = my_ebeam.ptcls[:,1]*1e6
gb = my_ebeam.ptcls[:,5]

gb0     = np.mean(gb)
dgb     = (gb-gb0)/gb0
std_dgb = np.std(dgb)

std_x  = np.std(x)
std_xp = np.std(xp)

x_pad  = 3.5*std_x
xp_pad = 3.5*std_xp

Cx  = np.mean(x)
Cxp = np.mean(xp)

x_range  = [Cx-x_pad,Cx+x_pad]
xp_range = [Cxp-xp_pad,Cxp+xp_pad]

cmap0 = plt.get_cmap('coolwarm')

scat0 = plt.scatter(x,xp,c=gb,s=0.75,cmap=cmap0)

ax1.set_facecolor('black')

plt.xlim(x_range)
plt.ylim(xp_range)

plt.xlabel(r'$x \,{\rm (\mu m)}$')
plt.ylabel(r'$x^\prime \,{\rm (\mu rad)}$')


nslice = 5

dgb_edges = np.linspace(-2*std_dgb,+2*std_dgb,nslice+1)

for islice in range(nslice):
    slice_min = dgb_edges[islice]
    slice_max = dgb_edges[islice+1]
    irange = (dgb>slice_min)&(dgb<slice_max)
    eps, beta, alpha, gamma = calc_Twiss(x[irange],xp[irange])
    
    H      = (beta+gamma)/2
    maj_ax = np.sqrt(eps/2)*(np.sqrt(H+1)+np.sqrt(H-1))
    min_ax = np.sqrt(eps/2)*(np.sqrt(H+1)-np.sqrt(H-1))
    phi    = (1/2)*np.arctan2(-2*alpha,beta-gamma)    
    
    el  = patches.Ellipse((Cx, Cxp), 3*min_ax, 3*maj_ax,
                          angle=-(phi+np.pi/2)*180/np.pi,
                          linewidth=3, fill=False,
                          edgecolor=cmap0(islice/(nslice-1)))

    ax1.add_patch(el)


    zw = -np.sign(alpha)*np.sqrt((beta-(1/gamma))/gamma)

    print(f'zw (cm) = {100*zw :.2f}')


#%%
# introduce offset into dumpline
#my_ebeam.ptcls[:,0] += -100e-6
#my_ebeam.ptcls[:,1] += -100e-6
#my_ebeam.ptcls[:,2] += -100e-6
#my_ebeam.ptcls[:,3] += -100e-6


#%%
# propagate forward to the end of the plasma section

interactions.electron_vacuum(my_ebeam, np.linspace(0,z_wx,2), dumpPeriod, numThread)
    

#%%
# plot transverse phase space at plasma exit virtual waist

def calc_Twiss(x,xp,weights=None):
    cen_x  = np.average(x, weights=weights)
    dx     = x - cen_x
    cen_xp = np.average(xp, weights=weights)
    dxp    = xp - cen_xp
    # Calculate the RMS sizes and the correlation
    sigmax2  = np.average(dx**2,  weights=weights)
    sigmaxp2 = np.average(dxp**2, weights=weights)
    sigmaxxp = np.average(dx*dxp, weights=weights)
    # Calculate the emittance
    ex = np.sqrt(sigmax2*sigmaxp2 - sigmaxxp**2)
    # Assign the outputs
    eps   = ex
    beta  = sigmax2/ex
    alpha = sigmaxxp/ex
    gamma = sigmaxp2/ex
    return eps, beta, alpha, gamma


fig1, ax1 = plt.subplots(figsize=(6, 6))

x  = my_ebeam.ptcls[:,0]*1e6
xp = my_ebeam.ptcls[:,1]*1e6
gb = my_ebeam.ptcls[:,5]

gb0     = np.mean(gb)
dgb     = (gb-gb0)/gb0
std_dgb = np.std(dgb)

std_x  = np.std(x)
std_xp = np.std(xp)

x_pad  = 3.5*std_x
xp_pad = 3.5*std_xp

Cx  = np.mean(x)
Cxp = np.mean(xp)

x_range  = [Cx-x_pad,Cx+x_pad]
xp_range = [Cxp-xp_pad,Cxp+xp_pad]

cmap0 = plt.get_cmap('coolwarm')

scat0 = plt.scatter(x,xp,c=gb,s=0.75,cmap=cmap0)

ax1.set_facecolor('black')

plt.xlim(x_range)
plt.ylim(xp_range)

plt.xlabel(r'$x \,{\rm (\mu m)}$')
plt.ylabel(r'$x^\prime \,{\rm (\mu rad)}$')


nslice = 5

dgb_edges = np.linspace(-2*std_dgb,+2*std_dgb,nslice+1)

for islice in range(nslice):
    slice_min = dgb_edges[islice]
    slice_max = dgb_edges[islice+1]
    irange = (dgb>slice_min)&(dgb<slice_max)
    eps, beta, alpha, gamma = calc_Twiss(x[irange],xp[irange])
    
    H      = (beta+gamma)/2
    maj_ax = np.sqrt(eps/2)*(np.sqrt(H+1)+np.sqrt(H-1))
    min_ax = np.sqrt(eps/2)*(np.sqrt(H+1)-np.sqrt(H-1))
    phi    = (1/2)*np.arctan2(-2*alpha,beta-gamma)    
    
    el  = patches.Ellipse((Cx, Cxp), 3*min_ax, 3*maj_ax,
                          angle=-(phi+np.pi/2)*180/np.pi,
                          linewidth=3, fill=False,
                          edgecolor=cmap0(islice/(nslice-1)))

    ax1.add_patch(el)


    zw = -np.sign(alpha)*np.sqrt((beta-(1/gamma))/gamma)

    print(f'zw (cm) = {100*zw :.2f}')




#fig2, ax2 = plt.subplots(figsize=(6, 6))
#
#x  = my_ebeam.ptcls[:,0]*1e6
#y  = my_ebeam.ptcls[:,2]*1e6
#gb = my_ebeam.ptcls[:,5]
#
#gb0     = np.mean(gb)
#dgb     = (gb-gb0)/gb0
#std_dgb = np.std(dgb)
#
#std_x  = np.std(x)
#x_pad  = 3.5*std_x
#Cx  = np.mean(x)
#
#std_y  = np.std(y)
#y_pad  = 3.5*std_y
#Cy  = np.mean(y)
#
#x_range = [Cx-x_pad,Cx+x_pad]
#y_range = [Cy-y_pad,Cy+y_pad]
#
#cmap0 = plt.get_cmap('coolwarm')
#
#scat0 = plt.scatter(x,y,c=gb,s=2.0,cmap=cmap0)
#
#ax2.set_facecolor('black')
#
#plt.xlim(x_range)
#plt.ylim(y_range)
#
#plt.xlabel(r'$x \,{\rm (\mu m)}$')
#plt.ylabel(r'$y \,{\rm (\mu m)}$')


#%%
# define the beamline

LD0 = 1.6000 -z_wx # m, 1.6 m from virtual waist
LQ0 = 1.0000 # m
LD1 = 1.2214 # m
LQ1 = 1.0000 # m
LD2 = 1.2213 # m
LQ2 = 1.0000 # m
LD3 = 3.5215 # m
LB5 = 0.9780 # m
LD4 = 8.8530 # m

fact  = beamE/(10e9)
G0    = fact*12.3551 #fact*0.9962*12.4022 # T/m
FxDx0 = 'Dx'
G1    = fact*19.3048 #fact*1.0074*19.3048 # T/m
FxDx1 = 'Fx'
G2    = G0
FxDx2 = FxDx0

B5    = 0.2047 # T


#%%

# drift D0
drift0_Params  = {'name' : 'D0',  'L'    : LD0}
drift0 = drift.Drift(drift0_Params)

# quad Q0
quad0_Params   = {'name' : 'Q0',  'L'    : LQ0, 'FxDx' : FxDx0, 'G'    : G0}
quad0 = quadrupole.Quadrupole(quad0_Params)

# drift D1
drift1_Params  = {'name' : 'D1',  'L'    : LD1}
drift1 = drift.Drift(drift1_Params)

# quad Q1
quad1_Params   = {'name' : 'Q1',  'L'    : LQ1, 'FxDx' : FxDx1, 'G'    : G1}
quad1 = quadrupole.Quadrupole(quad1_Params)

# drift D2
drift2_Params  = {'name' : 'D2',  'L'    : LD2}
drift2 = drift.Drift(drift2_Params)

# quad Q2
quad2_Params   = {'name' : 'Q2',  'L'    : LQ2, 'FxDx' : FxDx2, 'G'    : G2}
quad2 = quadrupole.Quadrupole(quad2_Params)

# drift D3
drift3_Params  = {'name' : 'D3',  'L'    : LD3}
drift3 = drift.Drift(drift3_Params)

# dipole B5
dipole5_Params = {'name' : 'B5',  'L'    : LB5, 'B'    : B5}
dipole5 = dipole.Dipole(dipole5_Params)

# drift D4
drift4_Params  = {'name' : 'D4',  'L'    : LD4}
drift4 = drift.Drift(drift4_Params)

# beamline
bline_Params = {'name' : 'ImgSpect'}
bline = beamline.Beamline(bline_Params)

bline.add_element(drift0)
bline.add_element(quad0)
bline.add_element(drift1)
bline.add_element(quad1)
bline.add_element(drift2)
bline.add_element(quad2)
bline.add_element(drift3)
bline.add_element(dipole5)
bline.add_element(drift4)

bline.initialize_matrix(gb0)
R = bline.get_R_matrix()
T6 = bline.get_T6_matrix()

#%%
# propagate the beam through beamline
interactions.ebeam_beamline(my_ebeam, bline, dumpPeriod, numThread)


#%%
# Beam profile

x = my_ebeam.ptcls[:,0] # convert to microns
y = my_ebeam.ptcls[:,2] # convert to microns

#yp = my_ebeam.ptcls[:,3]

Cx = np.mean(x)
Cy = np.mean(y)

std_x = np.std(x)
std_y = np.std(y)

fig1 = plt.figure()
ax1  = fig1.add_subplot(111)

x_pad = 3*std_x
y_pad = 0.005*bline.R[2,5] # 5*std_y

x_range = [Cx-x_pad,Cx+x_pad]
y_range = [Cy-1.0*y_pad,Cy+1.0*y_pad]

#x_range = [-0.000150,+0.000150]
#y_range = [-0.0005,0.0005]

px_size = 4.5e-6 # um

Nbinx = int((x_range[1]-x_range[0])/px_size)
Nbiny = int((y_range[1]-y_range[0])/px_size)

cmap1 = plt.get_cmap('PuBu_r')

hst1 = ax1.hist2d(x, y, bins=(Nbinx, Nbiny), range=[x_range,y_range], cmap=cmap1)

ax1.set_xlabel(r'$x\mathrm{\ [\mu m]}$')
ax1.set_ylabel(r'$y\mathrm{\ [\mu m]}$')

print('std x: ',std_x)
print('std y: ',std_y)


#%%
# Butterfly fit

# calculate sigma_x^2 as a function of y
Cx   = np.zeros(Nbiny)
stdx = np.zeros(Nbiny)
for iy in range(0,Nbiny):
    
    Cx[iy] = 0
    wtot   = 0
    ww     = np.zeros(Nbinx)
    xx     = np.zeros(Nbinx)
    for ix in range(0,Nbinx):
        ww[ix] = hst1[0][ix][iy]
        xx[ix] = hst1[1][ix]
        Cx[iy] = Cx[iy] + ww[ix]*xx[ix]
        wtot = wtot + ww[ix]
    if wtot>0:
        Cx[iy] = Cx[iy]/wtot

    stdx[iy] = 0
    for ix in range(0,Nbinx):
        stdx[iy] = stdx[iy] + ww[ix]*((xx[ix]-Cx[iy])**2)
    if wtot>0:
        stdx[iy] = np.sqrt(stdx[iy]/wtot)
 
ybin_w = np.abs(hst1[2][1]-hst1[2][0])
ybin_C = hst1[2][0:len(hst1[2])-1]+0.5*ybin_w  
R36  = bline.R[2,5]
dgb  = ybin_C/R36

Mx    = bline.R[0,0]
stdx0 = np.min(stdx)

stdx2 = stdx**2


# butterfly 1: first order fit
def butterfly(dgb,eps0,beta0,alpha0,gb0,R11,T116,R12,T126):
    gamma0 = (1+alpha0**2)/beta0
    return  (eps0/gb0)*((R11**2)*beta0 - 2*R11*R12*alpha0 + (R12**2)*gamma0
            +dgb*2*(R11*T116*beta0 - (R11*T126+R12*T116)*alpha0 + R12*T126*gamma0) 
            +(dgb**2)*((T116**2)*beta0 - 2*T116*T126*alpha0 + (T126**2)*gamma0))

def scaled_butterfly(dgb,eps0,beta0,alpha0,gb0,R11,T116,R12,T126):
    return butterfly(dgb,eps0*(1e-6),beta0*(1e-2),alpha0,gb0,R11,T116,R12,T126)
    
    


# introduce error in centroid energy
gb0_err = 1.000 #0.947




bline.initialize_matrix(gb0_err*gb0)

R11  = bline.R[0,0]
R12  = bline.R[0,1]
T116 = bline.T6[0,0]
T126 = bline.T6[0,1]

mod_err = 1.00
R11  *= mod_err
R12  *= mod_err
T116 *= mod_err
T126 *= mod_err


params, params_cov = optimize.curve_fit(lambda x, eps0, beta0, alpha0: 
    scaled_butterfly(x,eps0,beta0,alpha0,gb0_err*gb0,R11,T116,R12,T126), dgb, stdx2,
    p0=[eps_n0/(1e-6),beta_x0/(1e-2),0]) #,bounds=[[0,0,-1],[100,1,1]]) #,method='dogbox') #,bounds=[[-np.inf,0,-10],[np.inf,1,2]])
    
#alpha0 = 0
#params, params_cov = optimize.curve_fit(lambda x, eps0, beta0: 
#    butterfly(x,eps0,beta0,alpha0,gb0,R11,T116,R12,T126), dgb, stdx2,
#    p0=[(3e-6),0.05],bounds=[[1e-9,1e-3],[1e-4,1]]) 

print('\n')    
print(params)

print('emittance = ',params[0])
print('beta0 = ',params[1])
print('alpha0 = ',params[2])
print('waist = ',params[1]*params[2])

fig2 = plt.figure()
ax2  = fig2.add_subplot(111)
ax2.plot(dgb,stdx2)
fit_dgb = np.linspace(dgb[0],dgb[-1],100)
ax2.plot(fit_dgb,scaled_butterfly(fit_dgb,params[0],params[1],params[2],
                            gb0_err*gb0,R11,T116,R12,T126))
#ax2.plot(fit_dgb,butterfly(fit_dgb,params[0],params[1],alpha0,
#                            gb0,R11,T116,R12,T126))
ax2.set_xlabel(r'$\delta_p$')
ax2.set_ylabel(r'$\sigma_x^2 \mathrm{\ [m^2]}$')



"""
# butterfly 2: exact matrix fit
def butterfly2(dgb,eps0,beta0,alpha0,gb0,bline):
    # for array input
    if np.size(dgb)>1:
        R11 = np.zeros(np.size(dgb))
        R12 = np.zeros(np.size(dgb))
        gb  = np.zeros(np.size(dgb))
        for i in range(0,np.size(dgb)):
            gb[i] = gb0*(1 + dgb[i])
            bline.initialize_matrix(gb[i])
            R = bline.get_R_matrix()
            R11[i] = R[0,0]
            R12[i] = R[0,1]
    # for single value input
    else:  
        gb = gb0*(1 + dgb)
        bline.initialize_matrix(gb)
        R = bline.get_R_matrix()
        R11 = R[0,0]
        R12 = R[0,1]

    R11  *= mod_err
    R12  *= mod_err

    gamma0 = (1+alpha0**2)/beta0
    sigx2 = (eps0/gb)*((R11**2)*beta0 - 2*R11*R12*alpha0 + (R12**2)*gamma0)
    return sigx2

def scaled_butterfly2(dgb,eps0,beta0,alpha0,gb0,bline):
    return butterfly2(dgb,eps0*(1e-6),beta0*(1e-2),alpha0,gb0,bline)

params, params_cov = optimize.curve_fit(lambda x, eps0, beta0, alpha0: 
    scaled_butterfly2(x,eps0,beta0,alpha0,gb0_err*gb0,bline), dgb, stdx2,
    p0=[eps_n0/(1e-6),beta_x0/(1e-2),0]) #,bounds=[[1e-9,0,-10],[1e-4,5,10]])
    
#alpha0 = 0
#params, params_cov = optimize.curve_fit(lambda x, eps0, beta0: 
#    butterfly2(x,eps0,beta0,alpha0,gb0,bline), dgb, stdx2,
#    p0=[(3e-6),0.05])

print('\n')
print(params)

print('emittance = ',params[0])
print('beta0 = ',params[1])
print('alpha0 = ',params[2])
print('waist = ',params[1]*params[2])

fig3 = plt.figure()
ax3  = fig3.add_subplot(111)
ax3.plot(dgb,stdx2)
fit_dgb = np.linspace(dgb[0],dgb[-1],100)
ax3.plot(fit_dgb,scaled_butterfly2(fit_dgb,params[0],params[1],params[2],
                            gb0_err*gb0,bline))
#ax3.plot(fit_dgb,butterfly2(fit_dgb,params[0],params[1],alpha0,
#                            gb0,bline))
ax3.set_xlabel(r'$\delta_p$')
ax3.set_ylabel(r'$\sigma_x^2 \mathrm{\ [m^2]}$')

"""
#%%
# print the elapsed time

toc=timeit.default_timer()
print('\n')
print('Wall clock run time: ',toc - tic) # elapsed time in seconds
