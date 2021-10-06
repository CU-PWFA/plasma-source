#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:43:08 2021

Loop over the cross section to track for all longitudinal pos:
    
This is a copy, but this newer version compares with the Control sim case
Also starting the move to putting everything on the hard drive

@author: chris
"""

import sys
sys.path.insert(0, "../../python")
#import os
import numpy as np
from vsim import plot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy import optimize
#import scipy.integrate as integrate
#import scipy.constants as const
import scipy.interpolate as interp

#superpath = /home/chris/Desktop/'
superpath = '/media/chris/New Volume/VSimRuns/NonuniformPlasma/'

"""
#SepGrads
cont_path = superpath + 'NERSC_Sep_Control2/'
cont_ind = 5
#path = superpath + 'NERSC_Sep_Grad/'
#ind = 9
#grad = 8e17
path = superpath + 'NERSC_Mar_Grad0.001/'
grad = 2.5e15
ind = 10
#path = superpath + 'NERSC_Mar_Grad0.01/'
#grad = 2.5e16
#ind = 10
#path = superpath + 'NERSC_Dec_Grad/'
#grad = 8e17
#ind = 9
#path = superpath + '/VELA_Oct_Grad/'
#grad = 2e17
#ind = 5

offset_arr = np.arange(-120,100.5,2)
tranExtent = 200        #tranextent for the sim
threshold = 150         #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 2e16           #central density for the sim
dx=1.2                  #dx for the sim
start = 100 #um         #for correcting where drive beam is
trunc = 15              #for trimming the back of the wake
left = 0; right = 200   #for plotting r+ and r-
wanchor = -9            #for finding where r=0 to anchor yoff fit
"""

"""
#Jan Grads
cont_path = superpath + 'NERSC_Jan_Grad_Control/'
cont_ind = 8
path = superpath + 'NERSC_Jan_Grad/'
grad = 4e17
ind = 8

offset_arr = np.arange(-80,220.5,5) #JanGrad Full
tranExtent = 400        #tranextent for the sim
threshold = 150         #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 1e16           #central density for the sim
dx=1.2#0.25             #dx for the sim
start = 250 #um         #for correcting where drive beam is
trunc = 10#15              #for trimming the back of the wake
left = 20; right = 250  #for plotting r+ and r-
wanchor = -5            #for finding where r=0 to anchor yoff fit
"""
"""
#JunGrads (3e17)
cont_path = superpath + 'NERSC_Jun_Grad_Control/'
cont_ind = 6
#path = superpath + 'NERSC_Jun_Grad/'
#grad = 6e18
#ind = 6
path = superpath + 'NERSC_Jun_Grad_High/'
grad = 1.2e19
ind = 6

offset_arr = np.arange(-70,140.5,2)
tranExtent = 75        #tranextent for the sim
threshold = 40          #for use in cross section algorithm, make larger than rp and smaller than plasma extent
npcase = 3e17           #central density for the sim
dx=0.5                  #dx for the sim
start = 57 #um         #for correcting where drive beam is
trunc = 13#15           #for trimming the back of the wake
left = 10; right = 80  #for plotting r+ and r-
wanchor = -14            #for finding where r=0 to anchor yoff fit
"""
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#

#August 2021 n=2e16 runs
"""
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
cont_path = superpath + 'NERSC_n2e16_g0/'
cont_ind = 5
path = superpath + 'NERSC_n2e16_g8e17/'
grad = 8e17
#path = superpath + 'NERSC_n2e16_g2e17/'
#grad = 2e17
#path = superpath + 'NERSC_n2e16_g2.5e16/'
#grad = 2.5e16
#path = superpath + 'NERSC_n2e16_g2.5e15/'
#grad = 2.5e15
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
"""
"""
#August Linear Gradients, n=1e16 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
cont_path = superpath + 'NERSC_n1e16_g0/'
cont_ind=5
path = superpath + 'NERSC_n1e16_g4e17/'
ind = 5
grad = 4e17

tranExtent = 200
dx = 1.233 #um
start = 107 #um
threshold = 150
npcase = 1e16
offset_arr = np.arange(-140,90.5,2)
trunc = 5              #for trimming the back of the wake
left = 0; right = 200   #for plotting r+ and r-
wanchor = -9            #for finding where r=0 to anchor yoff fit
"""

#August Linear Gradients, n=1e17 sims
superpath = '/media/chris/New Volume/VSimRuns/AugustLinearGradient/'
cont_path = superpath + 'NERSC_n1e17_g0/'
cont_ind=5
path = superpath + 'NERSC_n1e17_g4e18/'
ind = 5
grad = 4e18
#path = superpath + 'NERSC_n1e17_g2e18/'
#ind = 5
#grad = 2e18

tranExtent = 95
dx = 0.5 #um
start = 45 #um
threshold = 70
npcase = 1e17
offset_arr = np.arange(-210,140.5,3)
trunc = 23              #for trimming the back of the wake
left = 0; right = 120   #for plotting r+ and r-
wanchor = -14#-15            #for finding where r=0 to anchor yoff fit

reset = False #Set to False to not reload all control data
# Define model function to be used to fit to the data above:
def coslike(x, *p):
    rp, b, re = p
    return rp-re*np.cos(x+b)

#offset_arr = np.arange(-80,90.5,2)
rmajor_arr = np.zeros(len(offset_arr))
rminor_arr = np.zeros(len(offset_arr))
n0sh_arr = np.zeros(len(offset_arr))
dndysh_arr = np.zeros(len(offset_arr))
if reset:
    rcontr_arr = np.zeros(len(offset_arr))
    rerror_arr = np.zeros(len(offset_arr))
    n0sh_arr_cont = np.zeros(len(offset_arr))
    dndysh_arr_cont = np.zeros(len(offset_arr))

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
        radmax = xcoeff[0]*1.05
    
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
              'radmax' : threshold,
              'drive' : 'rhoBeam'
              }
        x,y,n = plot.wake_cross_section_sheath(params_sh)
        sfit = np.polyfit(y,n,1)
        n0sh_arr[i] = sfit[1]
        dndysh_arr[i] = sfit[0]*1e4
    
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
                  'radmax' : threshold,
                  'drive' : 'rhoBeam'
                  }
            x,y,n = plot.wake_cross_section_sheath(params_sh)
            sfit = np.polyfit(y,n,1)
            n0sh_arr_cont[i] = sfit[1]
            dndysh_arr_cont[i] = sfit[0]*1e4

#Showing the cross section radius is the same
plt.plot(offset_arr*dx*-1+start,rmajor_arr,label="Sim")
plt.plot(offset_arr*dx*-1+start,rcontr_arr,c='k',ls='dotted',label="Control")
plt.title("Cross-sectional wake radius")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Wake radius "+r'$(\mu m)$')
plt.grid(); plt.legend(); plt.show();

x = offset_arr[trunc:]*dx*-1
y = rminor_arr[trunc:]
#weights = None
weights = np.zeros(len(offset_arr))+1; weights[wanchor]= 10; weights = weights[trunc:]
h = np.polyfit(x,y,1,w = weights)

radmax = max(rcontr_arr)
varA = 0.324 #Empirical Constant
varA = 0.30
rpl = radmax*np.sqrt((npcase + varA*1e-4*radmax*grad)/npcase)
rmn = radmax*np.sqrt((npcase - varA*1e-4*radmax*grad)/npcase)

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

plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh_arr[trunc:],label="Simulation")
plt.plot(offset_arr[trunc:]*dx*-1+start,n0sh_arr_cont[trunc:],c='k',ls='dotted',label="Control")
plt.title("Sheath's central density longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath's central density "+r'$(cm^{-3})$')
plt.grid(); plt.legend(); plt.show();

centloc = np.argmax(rcontr_arr)
offdlt = offset_arr[1]-offset_arr[0] #Since we don't load every slice
leftend = centloc + int(0.5*(offset_arr[centloc]*dx*-1+start)/(offdlt*dx))
rightend = centloc - int((offset_arr[centloc]*dx*-1+start)/(offdlt*dx))

dndysh_ave = np.average(dndysh_arr[rightend:leftend])
var3 =-3.0789e-38#-3.1006e-38
var2 = 6.1568e-19#6.1973e-19
var1 = 1.0266#1.0085
var0 = 2.5089e16#3.5021e16
dndysh_mod = grad**3*(var3) + grad**2*(var2) + grad*(var1) + var0
trunc_sub = 5 #for when trunc doesnt get to the end of the wake

plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,dndysh_arr[trunc-trunc_sub:],label="Simulation")
plt.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_ave,dndysh_ave],c='r',ls='-',label="Selected Average")
plt.plot([offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start],[dndysh_mod,dndysh_mod],c='g',ls='--',label="Empirical Model")
plt.plot(offset_arr[trunc-trunc_sub:]*dx*-1+start,dndysh_arr_cont[trunc-trunc_sub:],c='k',ls='dotted',label="Control")
plt.vlines((offset_arr[leftend]*dx*-1+start,offset_arr[rightend]*dx*-1+start),min(dndysh_arr[trunc:]),max(dndysh_arr[trunc:]),colors='r')
plt.title("Sheath density gradient longitudinally")
plt.xlabel("Distance behind drive beam "+r'$(\mu m)$')
plt.ylabel("Sheath density gradient "+r'$(cm^{-4})$')
plt.grid(); plt.legend(); plt.show();

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
print("dn/dy_sh average = ",dndysh_ave," cm-4")





