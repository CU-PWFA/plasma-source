#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  19 13:02:41 2017

Loops over flattop density for tolerance measurements

@author: chris
"""

import numpy as np
import PlasmaPropagation as PProp
import matplotlib.pyplot as plt

debug = 1
gamma_set = 2e4

n_center = 2.3e12*gamma_set

#npl0_arr = np.linspace(1e16, 2e17, num = 30)
npl0_arr = np.linspace(0.2*n_center, 4*n_center, num = 300)
#npl0_arr = np.linspace(1/2*n_center, 3/2*n_center, num = 300)
bmag_image = np.zeros((len(npl0_arr)))
nparts = 0

for i in range(len(npl0_arr)):
    dens = npl0_arr[i]
    
    #Make beam and bulk plasma just as in single_pass
    params = PProp.ReturnDefaultParams(npl0_change = dens, gbC_change = gamma_set)
    params['npart'] = nparts #Want 10,000 to 100,000 for a final figure
    params['L_ft'] = 0; params['L_dn'] = 0
    
    twiss = PProp.CallMakeTwiss(params)
    parts = PProp.CallMakeParts(twiss, params)
    ebeam0 = PProp.CallMakeBeam(twiss, parts, params)
    ebeam0 = PProp.PropagateBackwards(ebeam0, params)
    
    plasma0 = PProp.MakeBulkPlasma(params)
    ebeam0 = PProp.PropagatePlasma(ebeam0, plasma0)
    
    Bmag = PProp.CalcBmag(ebeam0, plasma0)
    if debug == 1: print('npl0: ',dens,'Bmag: ',Bmag)
    bmag_image[i] = Bmag
    
minloc = np.argmin(bmag_image)
n0min = npl0_arr[minloc]
tols = [1.10,1.05,1.01]
for tol in tols:
    print("Tolerance for "+str(tol)+":")
    result = PProp.Calc1DToleranceRange(bmag_image, npl0_arr[1]-npl0_arr[0], minloc, tol)
    print(" in +/- n:",result[0]/2)
    print(" rel tol :",result[0]/2/n0min)
    
dn = npl0_arr[1]-npl0_arr[0]
left = npl0_arr[0]+dn*result[1]
right = npl0_arr[0]+dn*result[2]

plt.plot(npl0_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([npl0_arr[0],npl0_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.scatter(n0min,min(bmag_image),color='k',label='min n0 = '+str(n0min))
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Flattop Density")
plt.xlabel(r'$ n_0 $ [cm^-3]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()

plt.semilogx(npl0_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([npl0_arr[0],npl0_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.scatter(n0min,min(bmag_image),color='k',label='min n0 = '+str(n0min))
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Flattop Density")
plt.xlabel(r'$ n_0 $ [cm^-3]')
plt.ylabel("B-mag")
plt.grid(); plt.legend(); plt.show()

plt.plot(npl0_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([npl0_arr[0],npl0_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.scatter(n0min,min(bmag_image),color='k',label='min n0 = '+str(n0min))
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Flattop Density")
plt.xlabel(r'$ n_0 $ [cm^-3]')
plt.ylabel("B-mag")
plt.grid(); plt.show()

plt.semilogx(npl0_arr, bmag_image,label='gamma='+str(gamma_set))
plt.plot([npl0_arr[0],npl0_arr[-1]],[tol,tol],ls='--',label='B-mag = 1.01')
plt.scatter(n0min,min(bmag_image),color='k',label='min n0 = '+str(n0min))
plt.axvspan(left, right, facecolor = 'g', alpha = 0.2,label = '+/- '+str(result[0]/2))
plt.title("B-mag vs Flattop Density")
plt.xlabel(r'$ n_0 $ [cm^-3]')
plt.ylabel("B-mag")
plt.grid(); plt.show()