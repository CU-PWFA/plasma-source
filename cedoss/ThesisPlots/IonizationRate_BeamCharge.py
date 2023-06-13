#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 11:26:33 2023

Making plots of ionization rate vs laser intensity and beam charge for my thesis

This version is for the ebeam

@author: chris
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "../../python")

from ionization import ionization
from ionization import adk

epsilon0 = 8.854e-12
lightspeed = 	299792458
echarge = 1.602e-19

#Want to interate over different three different taus and four
# different species

species_list=np.array([r'$H_2$',r'$Ar$',r'$Ar^+$',r'$He$'])
chi_list=np.array([15.4, 15.8, 24.6, 27.6])
color_list=np.array(['b','r','g','orange'])
text_height = np.array([.5, .5, .3, .3])
text_pos = np.array([3e17,5e17,3e17,5e17])#np.array([1.6e14,3e14,8e14,2e15])

sig_z_list = np.array([5,20,50])
radius_list = np.array([10,50,100])
ls_list = np.array(['solid','dashed','dotted'])

plt.figure(figsize=(8.5,3))

for m in range(len(species_list)):
    for n in range(len(ls_list)):
        
        #set up some necessary parameters
        chi = chi_list[m]
        sig_z = 20e-6 #m
        sig_xy = 5e-6
        delt_t = sig_z/lightspeed
        
        radius = radius_list[n]*1e-6
        
        #Following the function in ThreeDimensionAnalysis, which is
        # calling functions from Robert's python folder
        
        density_arr = np.logspace(17.4,19.6,1000)#cm-3
        charge_arr = np.linspace(0.5e-9,4e-9,50)#C
        ionfrac_arr = np.zeros(len(density_arr))
        
        for i in range(len(density_arr)):
            lambda0 = density_arr[i]*100**3*echarge*sig_xy**2*2*np.pi
            E = lambda0/(2*np.pi*epsilon0*radius)/1e9
            ionfrac_arr[i] = adk.gaussian_frac(chi,E,delt_t/1e-15,1,envelope=False)
        
        if m == 0:
            plt.semilogx(density_arr,ionfrac_arr,c=color_list[m],ls=ls_list[n],
                         label=r'$r=$'+str(radius_list[n])+"um")
        else:
            plt.semilogx(density_arr,ionfrac_arr,c=color_list[m],ls=ls_list[n])
        
        plt.text(text_pos[m], text_height[m] , species_list[m], color=color_list[m])
        
plt.legend(loc=0)
plt.xlabel("Peak Beam Density "+r'$\mathrm{(cm^{-3}})}$')
plt.ylabel("Ionization Fraction")
plt.show()
            
            
            