#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 08:47:40 2022

In a Single Directory, find the slice with the smallest sigma

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
#dataset = 'E308_02498' #24 psi
dataset = 'E308_02500' #24 psi

cut_head = 80
cut_tail = 130
y_roi = [40,200]
#y_roi = [0,267]

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
cut_arr = np.arange(cut_head,cut_tail+1,1)
min_sigma = np.zeros((len(cut_arr)))
sigma_data = np.zeros((n_step,n_shot))

#Triple for loop, going through all the steps and shots
for k in range(len(cut_arr)):
    cut_pos = cut_arr[k]
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
            try:
                xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)
                sigma_data[i,j] = np.abs(xcoeff[2])
            except RuntimeError:
                print("Optimal not found for cut_po s=",cut_pos,", step =",step,", n =",nimg)
                sigma_data[i,j] = None
                pass
    
    zob_arr = np.linspace(func_start,func_end,n_step)-func_start
    sigma_mean = np.zeros((n_step))
    sigma_var = np.zeros((n_step))
    for i in range(n_step):
        sigma_mean[i]=np.average(sigma_data[i,:])
        sigma_var[i] =np.var(sigma_data[i,:])
    
    min_sigma[k] = np.min(sigma_mean)

plt.plot(cut_arr,min_sigma)
plt.ylabel("Minimum Spot Size (pixels)")
plt.xlabel("Energy Slice (pixel)")
plt.show()










