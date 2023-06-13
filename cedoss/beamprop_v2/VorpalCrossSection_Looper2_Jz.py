#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:27:08 2022

Loop over the cross section to track for all longitudinal pos:
    
This is a copy, but this newer version compares with the Control sim case
Also starting the move to putting everything on the hard drive

Specifically, here I am looking at Jz

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
#import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
#from scipy import optimize
#import scipy.integrate as integrate
#import scipy.constants as const
import scipy.interpolate as interp
import scipy.constants as const

c = const.physical_constants['speed of light in vacuum'][0]
e = const.physical_constants['elementary charge'][0]

#superpath = /home/chris/Desktop/'
superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'

#August 2021 n=2e16 runs

superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
cont_path = superpath + 'NERSC_n2e16_g0/'
cont_ind = 5
path = superpath + 'JzDump/'
grad = 8e17
ind = 5

offset_arr = np.arange(-140,90.5,2)
tranExtent = 200        #tranextent for the sim
threshold = 150         #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
start = 90 #um         #for correcting where drive beam is
trunc = 20              #for trimming the back of the wake
left = 0; right = 200   #for plotting r+ and r-
wanchor = -9            #for finding where r=0 to anchor yoff fit

reset = True #Set to False to not reload all control data
eringcalc = True #Looking at the fields is a bit more
# Define model function to be used to fit to the data above:
def coslike(x, *p):
    rp, b, re = p
    return rp-re*np.cos(x+b)

#offset_arr = np.arange(-80,90.5,2)
rmajor_arr = np.zeros(len(offset_arr))
rminor_arr = np.zeros(len(offset_arr))
n0sh_arr = np.zeros(len(offset_arr))
dndysh_arr = np.zeros(len(offset_arr))
dndysh_jz = np.zeros(len(offset_arr))
dndysh_o2 = np.zeros(len(offset_arr))
if reset:
    rcontr_arr = np.zeros(len(offset_arr))
    rerror_arr = np.zeros(len(offset_arr))
    n0sh_arr_cont = np.zeros(len(offset_arr))
    dndysh_arr_cont = np.zeros(len(offset_arr))
if eringcalc:
    eyvx_0 = np.zeros(len(offset_arr))
    eyvy_0 = np.zeros(len(offset_arr))
    ey_v_y_theory_0 = np.zeros(len(offset_arr))
    ey_v_x_theory_0 = np.zeros(len(offset_arr))

for i in range(len(offset_arr)):
    vmax = 1e18*(npcase/3.7e17)
    vmin = 1e16*(npcase/3.7e17)
    central_off = offset_arr[i]
    
    
    params = {'vmin' : vmin,
              'vmax' : vmax,
              'plasma' : 'electrons',
              'dumpInd' : ind,
              'path' : path,
              'simName' : 'MatchedBeams',
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'threshold' : threshold,
              'plot' : False
              }
    theta,r = plot.wake_cross_section(params)
    #r, theta, rhoPol = plot.wake_cross_section(params)
            
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
            
    xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
    rmajor_arr[i] = xcoeff[0]
    rminor_arr[i] = xcoeff[2]
    
    if xcoeff[0] > 0: #Otherwise there is no wake
        radmax = xcoeff[0]*1.02
        yoff = xcoeff[2]
    
        params_sh = {'vmin' : vmin,
              'vmax' : vmax,
              'plasma' : 'electrons',
              'dumpInd' : ind,
              'path' : path,
              'simName' : 'MatchedBeams',
              'zoom' : 4.0,
              'alphaCutoff' : 0.05,
              'centralOff' : central_off,
              'tranExtent' : tranExtent,
              'plot' : False,
              'radmax' : radmax,
              'yoff' : yoff,
              'drive' : 'rhoBeam'
              }
        x,y,n = plot.wake_cross_section_sheath(params_sh)
        sfit = np.polyfit(y,n,1)
        n0sh_arr[i] = sfit[1]
        dndysh_arr[i] = sfit[0]*1e4
        pfit = np.polyfit(y,n,2)
        dndysh_o2[i] = pfit[0]
        
        #Here is where I load in Jz and Rho
        Jz_plasma = 'JPlasma'
        rho_plasma = 'rhoPlasma'
        vector = 0
        params = {'plasma' : 'electrons',
          'dumpInd' : ind,
          'path' : path,
          'simName' : 'MatchedBeams',
          'zoom' : 4.0,
          'alphaCutoff' : 0.05,
          'centralOff' : central_off,
          'tranExtent' : tranExtent,
          'plot' : False,
          'vector' : vector,
          'field' : Jz_plasma
        }
        xa, ya, evx, evy, JzYZ = plot.wake_cross_section_field(params)
        params['field'] = rho_plasma
        params['vector'] = 0
        xb, yb, bvx, bvy, rhoYZ = plot.wake_cross_section_field(params)
        
        sumDen = -(rhoYZ-JzYZ/c) / (e*1e6)#(e*npcase*100**3)
        
        den_v_x = np.flip(sumDen[int(len(evx)*1/4):int(len(evx)*3/4),int(len(evx)/2)],0)
        axis1 = np.where(den_v_x==0)[0]+int(len(evx)*1/4)
        nvals1 = np.zeros(len(axis1)*2)
        yvals1 = np.zeros(len(axis1)*2)
        for k in range(len(axis1)):
            den_v_y = np.flip(sumDen[:,axis1[k]],0)
            left_ind = np.argmax(den_v_y[:int(len(den_v_y)/2)])
            right_ind = np.argmax(den_v_y[int(len(den_v_y)/2):])+int(len(den_v_y)/2)
            nvals1[2*k] = den_v_y[left_ind]
            nvals1[2*k+1] = den_v_y[right_ind]
            yvals1[2*k] = ya[left_ind]
            yvals1[2*k+1] = ya[right_ind]
            
        den_v_x = np.flip(sumDen[int(len(evx)/2),int(len(evx)*1/4):int(len(evx)*3/4)],0)
        axis2 = np.where(den_v_x==0)[0]+int(len(evx)*1/4)
        nvals2 = np.zeros(len(axis2)*2)
        yvals2 = np.zeros(len(axis2)*2)
        for k in range(len(axis2)):
            den_v_y = np.flip(sumDen[axis2[k],:],0)
            left_ind = np.argmax(den_v_y[:int(len(den_v_y)/2)])
            right_ind = np.argmax(den_v_y[int(len(den_v_y)/2):])+int(len(den_v_y)/2)
            nvals2[2*k] = den_v_y[left_ind]
            nvals2[2*k+1] = den_v_y[right_ind]
            yvals2[2*k] = -ya[axis2[k]]
            yvals2[2*k+1] = -ya[axis2[k]]
        
        nvals = np.append(nvals1,nvals2)
        yvals = np.append(yvals1,yvals2)
        
        if len(nvals) == 0:
            dndysh_jz[i] = 0
        else:
            sfit_jz = np.polyfit(yvals,nvals,1)
            dndysh_jz[i] = sfit_jz[0]*1e4
    
    if reset:
        params = {'vmin' : vmin,
                  'vmax' : vmax,
                  'plasma' : 'electrons',
                  'dumpInd' : cont_ind,
                  'path' : cont_path,
                  'simName' : 'MatchedBeams',
                  'zoom' : 4.0,
                  'alphaCutoff' : 0.05,
                  'centralOff' : central_off,
                  'tranExtent' : tranExtent,
                  'threshold' : threshold,
                  'plot' : False
                  }
        theta,r = plot.wake_cross_section(params)
        #r, theta, rhoPol = plot.wake_cross_section(params)
                
        # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
        p0 = [1., 0., 1.]
                
        xcoeff, var_matrix = curve_fit(coslike, theta, r, p0=p0)
        rcontr_arr[i] = xcoeff[0]
        rerror_arr[i] = xcoeff[2]
        
        if xcoeff[0] > 0: #Otherwise there is no wake
            radmax = xcoeff[0]*1.05
            yoff = xcoeff[2]
        
            params_sh = {'vmin' : vmin,
                  'vmax' : vmax,
                  'plasma' : 'electrons',
                  'dumpInd' : cont_ind,
                  'path' : cont_path,
                  'simName' : 'MatchedBeams',
                  'zoom' : 4.0,
                  'alphaCutoff' : 0.05,
                  'centralOff' : central_off,
                  'tranExtent' : tranExtent,
                  'plot' : False,
                  'radmax' : radmax,
                  'yoff' : yoff,
                  'drive' : 'rhoBeam'
                  }
            x,y,n = plot.wake_cross_section_sheath(params_sh)
            sfit = np.polyfit(y,n,1)
            n0sh_arr_cont[i] = sfit[1]
            dndysh_arr_cont[i] = sfit[0]*1e4
            
    if eringcalc:
        #Get Ey_0
        params_e = {'plasma' : 'electrons',
                  'dumpInd' : ind,
                  'path' : path,
                  'simName' : 'MatchedBeams',
                  'zoom' : 4.0,
                  'alphaCutoff' : 0.05,
                  'centralOff' : central_off,
                  'tranExtent' : tranExtent,
                  'plot' : False,
                  'vector' : 1,
                  'field' : 'edgeE'#'ElecFieldPlasma'
                  }
                  
        x, y, evx, evy, eYZ = plot.wake_cross_section_field(params_e)
        params_e['field'] = 'faceB'
        params_e['vector'] = 2
        xb, yb, bvx, bvy, bYZ = plot.wake_cross_section_field(params_e)
        eYZ = eYZ - bYZ*const.speed_of_light
        Nz = len(x)
        
        eyvx_1 = np.array(eYZ[int((Nz+1)/2)-1,:])
        eyvy_1 = np.array(eYZ[:,int((Nz+1)/2)-1])
        eyvx_2 = np.array(eYZ[int((Nz+1)/2)-2,:])
        eyvy_2 = np.array(eYZ[:,int((Nz+1)/2)-2])
        eyvx = (eyvx_1 + eyvx_2)/2
        eyvy = (eyvy_1 + eyvy_2)/2
        
        eyvx_0[i] = eyvx[int((Nz/2))]
        eyvy_0[i] = eyvy[int((Nz/2))]
        
        #Calculate the expected using ions and the measured yoff
        radius = rcontr_arr[i]*1e-6
        n_cent = npcase*100**3
        slope = grad*100**4
        ybar = 1/4*radius**2*slope/n_cent
        yoff = -1*rminor_arr[i]*1e-6
        eps = 8.854e-12
        e = const.physical_constants['elementary charge'][0]
        pi = np.pi
        
        ey_v_y_theory_0[i] = 1/(2*eps)*e*n_cent*(-2*ybar-yoff) + 1/(8*eps)*e*slope*(3*(-yoff)**2)
        ey_v_x_theory_0[i] = 1/(8*eps)*e*slope*(3*yoff**2) + 1/eps*e*n_cent*(-ybar-1/2*yoff)

varA = 0.3239 #Empirical Constant
varB = 2.482
delta = 0.134
e = 1.602e-19
eps = 8.854e-12
"""
#trunc = 0
#Showing the cross section radius is the same
plt.plot(offset_arr*dx*-1+start,rmajor_arr,label="Sim")
plt.plot(offset_arr*dx*-1+start,rcontr_arr,c='k',ls='dotted',label="Control")
plt.title("Cross-sectional wake radius")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.legend(); plt.show();
"""
x = offset_arr[trunc:]*dx*-1
y = rminor_arr[trunc:]
#weights = None
weights = np.zeros(len(offset_arr))+1; weights[wanchor]= 10; weights = weights[trunc:]
h, residuals, rank, singular_values, rcond = np.polyfit(x,y,1,w = weights, full=True)

radmax = max(rcontr_arr)

#varA = 0.30
rpl = radmax*np.sqrt((npcase + varA*1e-4*radmax*grad)/npcase)
rmn = radmax*np.sqrt((npcase - varA*1e-4*radmax*grad)/npcase)
"""
#Showing the R+ and R- sheath trajectories
plt.plot(offset_arr*dx*-1+start,rmajor_arr+rminor_arr,c='b',label="Top")
plt.plot(offset_arr*dx*-1+start,rmajor_arr-rminor_arr,c='r',label="Bottom")
plt.plot(offset_arr*dx*-1+start,rcontr_arr,c='k',ls='dotted',label="Control")
plt.plot((offset_arr*dx*-1+start)*rpl/radmax,(rcontr_arr)*rpl/radmax,c='green',ls='dotted',label="Control+")
plt.plot((offset_arr*dx*-1+start)*rmn/radmax,(rcontr_arr)*rmn/radmax,c='purple',ls='dotted',label="Control-")
plt.plot([left,right],[rpl,rpl],c='b',ls='--',label=r'$R_+$')
plt.plot([left,right],[rmn,rmn],c='r',ls='--',label=r'$R_-$')
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2-rminor_arr2)
#plt.plot(offset_arr*dx*-1+97,rmajor_arr2+rminor_arr2)
plt.title("Wake radius in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.legend(); plt.show();
"""
""" ### Just using L_half now, but this was from the before times
#varB = 3.7#3.4#2.7?#*0.80 #Empirical Constant, relates length of wake to kp
varB = 2.7
#pi=const.pi
#re=2.8179e-13
#kp = np.sqrt(2*pi*re*npcase)
kp = 1.883e-6 * np.sqrt(npcase) *1e-4
yoff = 0.5*(rpl-rmn)/(varB*kp**-1)
"""
lhalf = offset_arr[np.argmax(rcontr_arr)]*dx*-1+start
yoff = 0.5*(rpl-rmn)/lhalf

"""
  ### Testing if subtracting splines could do it, it didn't
bound=1
spl_pl = interp.interp1d((offset_arr*dx*-1+start)*rpl/radmax,(rcontr_arr)*rpl/radmax,kind='cubic')
spl_mn = interp.interp1d((offset_arr*dx*-1+start)*rmn/radmax,(rcontr_arr)*rmn/radmax,kind='cubic')
pl = spl_pl((offset_arr[trunc:-bound]*dx*-1+start))
mn = spl_mn((offset_arr[trunc:-bound]*dx*-1+start))

plt.plot((offset_arr[trunc:-bound]*dx*-1+start),pl)
plt.plot((offset_arr[trunc:-bound]*dx*-1+start),mn)
plt.plot((offset_arr[trunc:-bound]*dx*-1+start),pl-mn,ls='--')
plt.show()
"""
"""
plt.plot(offset_arr[trunc:]*dx*-1+start,rminor_arr[trunc:], label = 'Sim')
plt.plot(offset_arr[trunc:]*dx*-1+start,rerror_arr[trunc:], c='k',ls='dotted',label = 'Control')
plt.plot(offset_arr[trunc:]*dx*-1+start,x*h[0]+h[1],ls='--', label = 'Fit')
plt.plot(offset_arr[trunc:]*dx*-1+start,x*yoff+start*yoff,ls='--', label = 'Emperical')
#plt.plot((offset_arr[trunc:-bound]*dx*-1+start),0.5*(pl-mn),ls='dashdot',label='Spline Subtration')
plt.ylim([(offset_arr[-1]*dx*-1+start)*yoff*1.1,(offset_arr[trunc]*dx*-1+start)*yoff*1.5])
#plt.ylim([-0.05,2.5])
plt.title("Wake vertical offset in first bucket")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake vertical offset "+r'$(\mu m)$')
plt.legend(); plt.grid(); plt.show();
"""
yoff_slope, yoff_intercept, yoff_r_value, yoff_p_value, yoff_std_err = stats.linregress(offset_arr[trunc:]*dx*-1,rminor_arr[trunc:])
"""
plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh_arr[trunc:],label="Simulation")
plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh_arr_cont[trunc:],c='k',ls='dotted',label="Control")
plt.title("Sheath's central density longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath's central density "+r'$(cm^{-3})$')
plt.grid(); plt.legend(); plt.show();
"""
centloc = np.argmax(rcontr_arr)
offdlt = offset_arr[1]-offset_arr[0] #Since we don't load every slice
leftend = centloc + int(0.5*(offset_arr[centloc]*dx*-1+start)/(offdlt*dx))
rightend = centloc - int(0.5*(offset_arr[centloc]*dx*-1+start)/(offdlt*dx))

dndysh_ave = np.average(dndysh_arr[rightend:leftend])
dndysh_std = np.std(dndysh_arr[rightend:leftend])
#var3 =-4.748e-37 #-3.0789e-38#-3.1006e-38
#var2 = 2.750e-18 # 6.1568e-19#6.1973e-19
#var1 =-0.0390    # 1.0266#1.0085
#var0 = 5.429e16  # 2.5089e16#3.5021e16
#dndysh_mod = grad**3*(var3) + grad**2*(var2) + grad*(var1) + var0
dndysh_mod = grad*varB

trunc_sub = 5 #for when trunc doesnt get to the end of the wake

plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,dndysh_arr[trunc-trunc_sub:],label="Simulation")
plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,-dndysh_jz[trunc-trunc_sub:],c='b',ls='--',label="Using Jz")
plt.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_ave,dndysh_ave],c='r',ls='-',label="Selected Average")
plt.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_mod,dndysh_mod],c='g',ls='--',label="Empirical Model")
plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,dndysh_arr_cont[trunc-trunc_sub:],c='k',ls='dotted',label="Control")
#plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,np.abs(dndysh_o2[trunc-trunc_sub:])*500000,ls='--',label="Quadratic Term")
plt.ylim([-4e17,5e18])
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
plt.title("Sheath density gradient longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath density gradient "+r'$(cm^{-4})$')
plt.grid(); plt.legend(); plt.show();

#ey0_arr = (eyvy_0+eyvx_0)/2
ey0_arr = eyvx_0
"""
plt.title("Longitudinal Variation of Electric Field at Origin")
plt.plot(offset_arr[trunc:]*dx*-1+start,ey_v_x_theory_0[trunc:],label="Ringless Theory")
plt.plot(offset_arr[trunc:]*dx*-1+start,eyvy_0[trunc:],label='Ey_0 (vs y)')
plt.plot(offset_arr[trunc:]*dx*-1+start,eyvx_0[trunc:],label='Ey_0 (vs x)')
plt.plot(offset_arr[trunc:]*dx*-1+start,ey0_arr[trunc:],label='Ey_0 (ave)')
plt.ylabel("Electric Field "+r'$(V/m)$')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.grid();plt.legend();plt.show()
"""
leftend = centloc + int(0.5*(offset_arr[centloc]*dx*-1+start)/(offdlt*dx))
rightend = centloc - int(0.3*(offset_arr[centloc]*dx*-1+start)/(offdlt*dx))

ering_ave = np.average(ey0_arr[rightend:leftend])
ering_std = np.std(ey0_arr[rightend:leftend])

func = (ey_v_x_theory_0-ey0_arr)/((rcontr_arr*1e-6)**(2))*-1
ering_ave = np.average(func[rightend:leftend])
ering_std = np.std(func[rightend:leftend])
"""
plt.plot(offset_arr[trunc:]*dx*-1+start,func[trunc:],label='Normalized Field')
plt.ylim([0,(func[int(len(func)/2)])*2])
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
plt.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[ering_ave,ering_ave],c='r',ls='-',label="Selected Average")
#plt.ylabel("Electric Field "+r'$(V/m)$')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.grid();plt.legend();plt.show()
"""
"""
func = (ey_v_x_theory_0[trunc:]-eyvx_0[trunc:])/(rcontr_arr[trunc:]**(2))/(-dndysh_arr[trunc:]*100**4-slope)
plt.plot(offset_arr[trunc:]*dx*-1+start,func/func[int(len(func)/2)],label='full')
func = (ey_v_x_theory_0[trunc:]-eyvx_0[trunc:])/(rcontr_arr[trunc:]**(2))
plt.plot(offset_arr[trunc:]*dx*-1+start,func/func[int(len(func)/2)],label='nosh')
plt.ylim([0,2])
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
#plt.ylabel("Electric Field "+r'$(V/m)$')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.grid();plt.legend();plt.show()

power_arr = np.linspace(1,3,5)
for i in range(len(power_arr)):
    power = power_arr[i]
    func = (ey_v_x_theory_0[trunc:]-eyvx_0[trunc:])/(rcontr_arr[trunc:]**(power))#/(-dndysh_arr[trunc:]*100**4-slope)
    plt.plot(offset_arr[trunc:]*dx*-1+start,func/func[int(len(func)/2)],label=power)
plt.ylim([0,2])
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
#plt.ylabel("Electric Field "+r'$(V/m)$')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.grid();plt.legend();plt.show()

for i in range(len(power_arr)):
    power = power_arr[i]
    func = (ey_v_y_theory_0[trunc:]-eyvy_0[trunc:])/(rcontr_arr[trunc:]**(power))#/(-dndysh_arr[trunc:]*100**4-slope)
    plt.plot(offset_arr[trunc:]*dx*-1+start,func/func[int(len(func)/2)],label=power)
plt.ylim([0,2])
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
#plt.ylabel("Electric Field "+r'$(V/m)$')
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.grid();plt.legend();plt.show()
"""
"""
fig, (ax0,ax1) = plt.subplots(nrows=2,ncols=1,sharex=True)
fig.set_size_inches(5,6)
"""
pltfac = 1e18
dndysh_arr_plt = dndysh_arr/pltfac
dndysh_ave_plt = dndysh_ave/pltfac
dndysh_mod_plt = dndysh_mod/pltfac
dndysh_std_plt = dndysh_std/pltfac
"""
pl0 = ax0.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,dndysh_arr_plt[trunc-trunc_sub:],label="Simulation")
ax0.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_ave_plt,dndysh_ave_plt],c='r',ls='-',label="Selected Average")
ax0.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_mod_plt,dndysh_mod_plt],c='g',ls='--',label="Empirical Model")
#ax0.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
ax0.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),dndysh_ave_plt-dndysh_std_plt,dndysh_ave_plt+dndysh_std_plt,colors='r')
ax0.set_ylim([-1e17/pltfac,5e18/pltfac])
#ax0.set_ylabel("Sheath density gradient "+r'$(cm^{-4})$')
ax0.set_ylabel(r'$(\partial n / \partial y )_{sh} \ \mathrm{(10^{18}\times cm^{-4})}$')
ax0.legend(title="(a)")
"""
E_si_to_cgs = (1/3)*10**(-4)
L_si_to_cgs = 100
factor = 1/1e17#E_si_to_cgs/L_si_to_cgs/1e9
func_scl = func*factor
ering_ave_scl = ering_ave*factor
ering_emp_scl = e/2/eps*delta*(((varB-1)*grad)*100**4)*factor
"""
pl0 = ax1.plot(offset_arr[trunc:]*dx*-1+start,func_scl[trunc:],label='Normalized Field')
ax1.set_ylim([0,(func_scl[int(len(func_scl)/2)])*2])
#ax1.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),20,60,colors='r')
ax1.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),ering_ave_scl-ering_std*factor,ering_ave_scl+ering_std*factor,colors='r')
ax1.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[ering_ave_scl,ering_ave_scl],c='r',ls='-',label="Selected Average")
ax1.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[ering_emp_scl,ering_emp_scl],c='g',ls='--',label="Empirical Model")
#plt.ylabel("Electric Field "+r'$(V/m)$')
ax1.set_xlabel("Distance behind drive beam "+r'$(\mu m)$')
#ax1.set_ylabel(r'$\langle E_{ring}/R_p^2\rangle \ \mathrm{(GstatV/cm^3)}$')
ax1.set_ylabel(r'$E_{sh}/R^2 \ \mathrm{(10^{17}\times V/m^3)}$')
ax1.legend(title="(b)")

fig.subplots_adjust(hspace=0)
"""
plt.show()

print("np_case = ",npcase, " cm-3")
print("dn/dy = ",grad, " cm-4")
print("R_max,control = ",max(rcontr_arr)," um")
print("L_half = ",offset_arr[np.argmax(rcontr_arr)]*dx*-1+start," um")
print("offset = ",offset_arr[np.argmax(rcontr_arr)])
print("R_major = ",rmajor_arr[np.argmax(rcontr_arr)])
print("y_off = ",rminor_arr[np.argmax(rcontr_arr)])
print("R_+ = ",max(rmajor_arr+rminor_arr)," um")
print("R_- = ",max(rmajor_arr-rminor_arr)," um")
print("yoff_fit = ",h[0])
print("y_off std.err = ",yoff_std_err)
print("dn/dy_sh average = ",dndysh_ave," cm-4")
print("dn/dy_sh std.dev. = ",dndysh_std," cm-4")
print("ering/rp2 average = ",ering_ave," V/m 1/m^2")
print("ering/rp2 std.dev. = ",ering_std," V/m 1/m^2")




