#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:50:30 2023

Making plots of ionization rate vs laser intensity and beam charge for my thesis

@author: chris
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "../../python")

from ionization import ionization
from ionization import adk

#Want to interate over different three different taus and four
# different species

species_list=np.array([r'$H_2$',r'$Ar$',r'$Ar^+$',r'$He$'])
chi_list=np.array([15.4, 15.8, 24.6, 27.6])
color_list=np.array(['b','r','g','orange'])
text_height = 0.5
text_pos = np.array([2e14,4e14,1e15,2.6e15])


delt_t_list = np.array([40,70,100])
ls_list = np.array(['solid','dashed','dotted'])

plt.figure(figsize=(8.5,3))

for m in range(len(species_list)):
    for n in range(len(delt_t_list)):
        
        #set up some necessary parameters
        chi = chi_list[m]
        delt_t = delt_t_list[n]*1e-15
        
        #Following the function in ThreeDimensionAnalysis, which is
        # calling functions from Robert's python folder
        
        intensity_arr = np.logspace(13.9,15.9,1000)
        ionfrac_arr = np.zeros(len(intensity_arr))
        
        for i in range(len(intensity_arr)):
            E = ionization.field_from_intensity(intensity_arr[i]/1e14)
            ionfrac_arr[i] = adk.gaussian_frac(chi,E,delt_t/1e-15,1)
        
        if m == 0:
            plt.semilogx(intensity_arr,ionfrac_arr,c=color_list[m],ls=ls_list[n],
                         label=r'$\tau_{FWHM}=$'+str(delt_t_list[n])+"fs")
        else:
            plt.semilogx(intensity_arr,ionfrac_arr,c=color_list[m],ls=ls_list[n])
        
        plt.text(text_pos[m], text_height , species_list[m], color=color_list[m])
        
plt.legend(loc=0)
plt.xlim([9e13, 5e15])
plt.xlabel("Peak Laser Intensity "+r'$\mathrm{(W/cm^2)}$')
plt.ylabel("Ionization Fraction")
plt.show()
            
            
            