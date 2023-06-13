#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 08:45:50 2022

@author: chris
"""
import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import scipy.io as sio
import scipy.constants as const

e = const.physical_constants['elementary charge'][0]

superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'
"""
#For lowres simulations
dataset_arr = np.array([
        'E308_02493',
        'E308_02512',
        'E308_02492',
        'E308_02498',
        'E308_02501',
        'E308_02506'])
pressure_arr = np.array([0,1,6,24,57.8,115.8])

toro_arr = np.array([
       'TORO_LI20_1988_0_TMIT',
       'TORO_LI20_2040_0_TMIT',
       'TORO_LI20_2452_0_TMIT',
       'TORO_LI20_3163_0_TMIT',
       'TORO_LI20_3255_0_TMIT'])

average_arr = np.zeros(  (len(dataset_arr) , len(toro_arr) )) 
var_arr = np.zeros(  (len(dataset_arr) , len(toro_arr) )) 
#For each dataset
for i in range(len(dataset_arr)):
    dataset = dataset_arr[i]

    path = superpath + day + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(path)
    
    data = mat['data_struct']
    #For each Toroid in S20 BSA  
    for j in range(len(toro_arr)):
        #Get Toroid data for each step
        scalars = data['scalars'][0][0]
        bsa = scalars['BSA_List_S20'][0][0]
        toro_name = toro_arr[j]
        toro = bsa[toro_name][0][0]
        toro = toro[:,0]*e*1e9
        average_arr[i][j] = np.average(toro[np.where(toro > 0)[0]])
        var_arr[i][j] = np.sqrt(np.var(toro[np.where(toro > 0)[0]]))
        
    #Plot average toroid signal vs toroid position
    """
colors = plt.cm.Set1(np.linspace(0, 1, len(toro_arr)))
plt.figure(figsize=(4.5,3.5))
for j in range(len(toro_arr)):
    plt.scatter(pressure_arr,average_arr[:,j],c=colors[j],label=toro_arr[j])
    plt.errorbar(pressure_arr,average_arr[:,j],yerr=var_arr[:,j],fmt="o",c=colors[j])
    if j%2==0:
        style='dashed'
    else:
        style='dashdot'
    plt.plot([pressure_arr[0]+.005*j,pressure_arr[-1]+200],[average_arr[0,j],average_arr[0,j]],ls=style,c=colors[j])
#plt.title("Toroid Signal vs Backing Pressure")
plt.xlabel("Backing Pressure (psi)")
plt.ylabel("Charge (nC)")
plt.xlim([0.5,170])
plt.xscale('log')

plt.legend(loc=2);plt.show()