#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:21:19 2020

@author: valentinalee
"""

#%%
import numpy as np
#%%
cross= 5e-19 #(m^2)
kb=1.38064852e-23  # J/K 
n0= 1e23 #m^-3
EV = 11604.5
Ti = 0.025*EV
Te = 2*EV
me= 9.11e-31
mi=4*1.67e-27
e= 1.6e-19
epsilon=8.85e-12
f_He_e= n0*cross*np.sqrt(kb*Te/me)
f_He_i= n0*cross*np.sqrt(kb*Ti/mi)
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me
f_e_e= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(me)*(kb*Te)**1.5))*10
f_i_i= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(mi)*(kb*Ti)**1.5))*10

print('Collision Frequency (1/s):')
print('CUPWFA_He_e='+"{:.2e}".format(f_He_e))
print('CUPWFA_He_i='+"{:.2e}".format(f_He_i))
print('CUPWFA_e_i='+"{:.2e}".format(f_e_i))
print('CUPWFA_e_e='+"{:.2e}".format(f_e_e))
print('CUPWFA_i_i='+"{:.2e}".format(f_i_i))

#%%Collision frequency HOFI
cross= 5e-19 #(m^2)
kb=1.38064852e-23  # J/K 
n0= 1e24 #m^-3
EV = 11604.5
Ti = 0.025*EV
Te = 13.7*EV
me= 9.11e-31
mi=1*1.67e-27
e= 1.6e-19
epsilon=8.85e-12
f_H_e= n0*cross*np.sqrt(kb*Te/me)
f_H_i= n0*cross*np.sqrt(kb*Ti/mi)
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me
f_e_e= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(me)*(kb*Te)**1.5))*10
f_i_i= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(mi)*(kb*Ti)**1.5))*10

print('HOFI_H_e='+"{:.2e}".format(f_H_e))
print('HOFI_H_i='+"{:.2e}".format(f_H_i))
print('HOFI_e_i='+"{:.2e}".format(f_e_i))
print('HOFI_e_e='+"{:.2e}".format(f_e_e))
print('HOFI_i_i='+"{:.2e}".format(f_i_i))

