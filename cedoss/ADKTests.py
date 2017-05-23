# -*- coding: utf-8 -*-
"""
Created on Wed May 17 17:49:06 2017

Tests the ADK_Combined.py code.  With print_table=1 information on the
ionization levels and their respective intensities will be printed out
in a table form.  A plot is also created at the end to overplot the
ionization curves for H, He, Ar, and Ar+2

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import ADK_Combined as adk

#Enter 1 to print out information on individual ionization levels
print_table=0

#Ionization energies [eV] for species
chi_H = adk.Get_chi_H()
chi_He = adk.Get_chi_He()
chi_Ar = adk.Get_chi_Ar1()
chi_Ar2 = adk.Get_chi_Ar2()

#We are using a square pulse of width 100fs
delt_t = 100e-15

#Setup the arrays for electric fields.  t1 is for cosmetics and t2 is more precise
start = 17
end = 117
step = (end-start)/40
t1 = np.arange(start, end, step)
t2 = np.arange(start, end, step/100)

#Get all the information on all the species if print_table=1
if print_table == 1:
    adk.get_Table(t2,delt_t,chi_H,1)
    adk.get_Table(t2,delt_t,chi_He,1)
    adk.get_Table(t2,delt_t,chi_Ar,1)
    adk.get_Table(t2,delt_t,chi_Ar2,2)

#Plots the ionization curves for all species
int_1 = adk.Intensity(t1)
int_2 = adk.Intensity(t2)
plt.figure(1)
plt.plot(int_1, adk.TotalProb(t1,delt_t,chi_H,1), 'bo', int_2, adk.TotalProb(t2,delt_t,chi_H,1), 'k',
         int_1, adk.TotalProb(t1,delt_t,chi_He,1), 'ro', int_2, adk.TotalProb(t2,delt_t,chi_He,1), 'k',
         int_1, adk.TotalProb(t1,delt_t,chi_Ar,1), 'go', int_2, adk.TotalProb(t2,delt_t,chi_Ar,1), 'k',
         int_1, adk.TotalProb(t1,delt_t,chi_Ar2,2), 'yo', int_2, adk.TotalProb(t2,delt_t,chi_Ar2,2), 'k')
plt.title("ADK Ionization of H(Blue), He(Red), Ar(Green 1st, Yellow 2nd)")
plt.ylabel("Ionization Fraction")
plt.xscale('log')
plt.xlabel("Laser Intensity (W/cm^2)")
plt.show()