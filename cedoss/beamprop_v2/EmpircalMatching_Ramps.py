#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:12:45 2019

Get empircal matching for various ramps.

At least that was the plan.  

@author: chris
"""

import BeamPropFuncs as PProp
import matplotlib.pyplot as plt
import numpy as np

"""Check to make sure energy gain is turned off in BeamPropFuncs"""

debug = 0
path = '/home/chris/Desktop/BeamProp/testGaussian'
gamma = PProp.def_gamma
n_0 = 0.5 #e17 cm^-3
"""Also check that n = this in BeamPropFuncs"""

z0=0
zvac = 0
betastar = 1/(1.33e-4*np.sqrt(n_0*1e17/gamma))

arrlen = 100
sigma_arr = np.linspace(5e-2*1e6, 50e-2*1e6, arrlen)
beta_arr = np.zeros(arrlen)
waist_arr = np.zeros(arrlen)
for i in range(len(sigma_arr)):
    sig_hw = sigma_arr[i] #um
    
    argon_params = PProp.ReturnDefaultPlasmaParams(path, scaledown = 10)
    argon_params['Z']=sig_hw * 5
    argon = PProp.GaussianExitRamp(argon_params, sig_hw, debug)
    
    n = argon.nez
    z = np.linspace(0,argon_params['Z'],len(n))/1e6
    
    beam_params = PProp.ReturnDefaultElectronParams(path, beta_star=betastar,
                                                   beta_offset=zvac, plasma_start=z0)
    
    betaf, alphaf, gammaf, gbf = PProp.Plot_CSEvo_FinalCompare(beam_params, n, z, z0, plot=debug)
    
    waistmat = alphaf[-1]*1/gammaf[-1] + z[-1]
    if debug == 1:
        print('Matching Beta:',1/gammaf[-1])
        print('Matching Waist:',waistmat)
    beta_arr[i]=1/gammaf[-1]
    waist_arr[i]=waistmat

sigma_arr=sigma_arr/1e6
waist_arr=waist_arr*-1
#plt.plot(beta_arr, sigma_arr)
#plt.plot(beta_arr, waist_arr)

#First is sigma vs beta
x_arr=np.linspace(beta_arr[0], beta_arr[-1], 200)
lfit1 = np.polyfit(beta_arr, sigma_arr, 1)
pfit1 = np.polyfit(beta_arr, sigma_arr, 2)
#lfit1=[1.65,-0.0170]

llabel1 = '{0:.2e}'.format(lfit1[0]) + r'$ \beta^* $'
llabel_simp1 = llabel1
llabel1 = llabel1 + " + " + '{0:.2e}'.format(lfit1[1])

plabel1 = '{0:.2e}'.format(pfit1[0]) + r'$ (\beta^*)^2 $'
plabel1 = plabel1 + " + " + '{0:.2e}'.format(pfit1[1]) + r'$ \beta^* $'
plabel1 = plabel1 + " + " + '{0:.2e}'.format(pfit1[2])

plt.scatter(beta_arr, sigma_arr, label="Contour Data")
plt.plot(x_arr, lfit1[0] * x_arr + lfit1[1], label=llabel1)
plt.plot(x_arr, pfit1[0] * np.square(x_arr) + pfit1[1] * x_arr + pfit1[2], label=plabel1)
plt.title("Matched ramp half-length vs beta")
plt.ylabel(r'$\sigma_{hw,m}$ [m]')
plt.xlabel(r'$\beta^*$[m]')
plt.grid(); plt.legend(); plt.show()

pfit = np.polyfit(beta_arr, waist_arr, 2)
#pfit=[-3.92,-4.67,0.067]
llabel2 = '{0:.2e}'.format(pfit[0]) + r'$ (\beta^*)^2 $'
llabel_simp2 = llabel2
llabel2 = llabel2 + " + " + '{0:.2e}'.format(pfit[1]) + r'$ \beta^* $'
llabel_simp3 = llabel2
llabel2 = llabel2 + " + " + '{0:.2e}'.format(pfit[2])

plt.scatter(beta_arr, waist_arr, label="Contour Data")
#plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr + pfit[2], label=llabel2)
plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr + pfit[2], label="all 3 terms")
plt.plot(x_arr, pfit[1] * x_arr + pfit[2], label="no quad")
plt.plot(x_arr, pfit[0] * np.square(x_arr), label="just quad")
plt.plot(x_arr, pfit[0] * np.square(x_arr) + pfit[1] * x_arr, label="no offset")
plt.title("Matched waist position vs beta")
plt.ylabel(r'$z_{\beta*,m}$ [m]')
plt.xlabel(r'$\beta^*$[m]')
plt.grid(); plt.legend(); plt.show()