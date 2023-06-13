#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 08:00:58 2022

Loop over all of the images in a single folder, and get a plot of sigma vs z_obj

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

superpath = '/home/chris/Desktop/SLACData/'
day = '20220812/'

#dataset = 'E308_02493' #no gas

#dataset = 'E308_02492' #6 psi
#dataset = 'E308_02497' #6 psi, highres

#dataset = 'E308_02498' #24 psi
#dataset = 'E308_02500' #24 psi, highres

#dataset = 'E308_02501' #5 bar over 0, 57.8 psi
#dataset = 'E308_02505' #5 bar over 0, 57.8 psi, highres

#dataset = 'E308_02506' #9 bar over 0, 115.8 psi
#dataset = 'E308_02511' #9 bar over 0, 115.8 psi, highres

dataset = 'E308_02512' #1 psi
dataset = 'E308_02513' #1 psi, highres

cut_pos = 130
y_roi = [40,200]
y_roi = [0,267]

path = superpath + day + dataset + '/' + dataset +'.h5'
f = h5.File(path,"r") 
data = f['data']
save_info = data['save_info']
params = data['params']
    
print("Comment:")
print(params['comment'][0][0])
print()
print("Scan Function: ",params['scanFuncs'][0][0])
func_start = params['startVals'][0][0]
func_end = params['stopVals'][0][0]
print(" ",func_start," to ",func_end)
print()

camerapath = 'DTOTR2'
path = superpath + day + dataset + '/images/' + camerapath + '/'
image_list = os.listdir(path)

n_shot = int(params['n_shot'][0][0])
n_step = int(len(image_list)/n_shot)

print("Shots per Step: ",n_shot)
print("Num Steps: ",n_step)
print()

f.close()

#Initialize a sigma array that is n_steps long with n_shots of data at each step
sigma_data = np.zeros((n_step,n_shot))

#Double for loop, going through all the steps and shots
for i in range(n_step):
    for j in range(n_shot):
        step = i+1
        nimg = j
        if step < 10:
            image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
        else:
            image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
        
        path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
        image = plt.imread(path)
        im_arr = np.array(image)
        
        slice_arr = im_arr[y_roi[0]:y_roi[1],cut_pos]
        slice_axs = np.arange(len(slice_arr))-int(.5*len(slice_arr))
        
        p0 = [100., 0., 3.]
        xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)
        
        sigma_data[i,j] = np.abs(xcoeff[2])

zob_arr = np.linspace(func_start,func_end,n_step)-func_start
sigma_mean = np.zeros((n_step))
sigma_var = np.zeros((n_step))
for i in range(n_step):
    sigma_mean[i]=np.average(sigma_data[i,:])
    sigma_var[i] =np.var(sigma_data[i,:])

plt.plot(zob_arr,sigma_mean,label="Sigma Average; Min = "+str(round(np.min(sigma_mean),2)))
plt.errorbar(zob_arr,sigma_mean,yerr=sigma_var,fmt="o",label="Variance per Step")
plt.title(dataset + ": "+camerapath+" Beam Size at Slice = "+str(cut_pos))
plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.ylabel(r'$\mathrm{\sigma}$'+" (pixels)")
plt.ylim([0.95*np.min(sigma_mean),1.05*np.max(sigma_mean)])
plt.legend();plt.show()










