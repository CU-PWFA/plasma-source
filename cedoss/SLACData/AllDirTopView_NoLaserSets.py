#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:12:19 2023

Here I want to plot lineaout vs pressure

Unfortunately, had to use 1psi for background since 0psi has no laser
At 1psi, used a nolaser background and just manually cropped out the edges
Doens't really matter, 1psi was mostly noise and no real plasma signal

Editing to look at tghe no-laser images

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
#from scipy import optimize
from scipy.optimize import curve_fit
import scipy.io as sio
import matplotlib.colors as colors

"""
def FitDataSomething(data, axis, function, guess = [0.,0.,0.]):
    errfunc = lambda p, x, y: function(p, x) - y
    p0 = guess
    p1, success = optimize.leastsq(errfunc,p0[:], args=(axis, data))
    return p1
"""

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def SuperGauss(x, *p):
    A, mu, sigma, n = p
    return A*np.exp(-np.power( (x-mu)**2/(2.*sigma**2),n ))

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day_names = np.array(['20220812/','20220811/'])

"""
#For lowres simulations
dataset_arr = np.array([
        #'E308_02493',
        'E308_02484',#no beam here, later dataset at 1 psi had no topview :(
        'E308_02492',
        'E308_02498',
        'E308_02501',
        'E308_02506'])
pressure_arr = np.array([1,6,24,57.8,115.8])#0 is boring
"""

dataset_arr = np.array([
        'E308_02483',
        'E308_02478',])
pressure_arr = np.array([1, 150])
day_arr=np.array([0,1])

smoothBack = True
doNormalize = False
colors = plt.cm.brg(np.linspace(0, 1, 3*len(dataset_arr)))
plt.figure(figsize=(6.5,4))

for i in range(len(dataset_arr)):
    dataset = dataset_arr[i]
    day = day_names[day_arr[i]]
    path = superpath + day + dataset + '/' + dataset +'.h5'
    
    ###BACKGROUND FROM MAT FILE
    
    #print("Loading Background:")
    path = superpath + day + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(path)
    
    data = mat['data_struct']
    backs = data['backgrounds'][0][0]
    back_arr = backs['TopView'][0][0]
    back_arr = np.transpose(back_arr)
    """
    plt.set_cmap('gray')
    plt.imshow(im_arr)
    CB = plt.colorbar(orientation='horizontal')
    plt.title("image"); plt.show()
    """
    params = data['params'][0][0]
    n_shot = params['n_shot'][0][0][0][0]
    #print("Shots",n_shot)
    n_shot = 9
    #######################
    im_arr2 = None
    for j in range(n_shot):
        step = 1; nimg = j;
        camerapath = 'TopView'
        #print("Loading Step",step,"at n",nimg,":")
        if step < 10:
            image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
        else:
            image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
        path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
        image = plt.imread(path)
        im_arr_single = np.array(image)
        
        for x in range(im_arr_single.shape[0]):
            for y in range(im_arr_single.shape[1]):
                im_arr_single[x][y] = im_arr_single[x][y] - back_arr[x][y]
                if im_arr_single[x][y] > 65000:
                    im_arr_single[x][y] = 0
        if im_arr2 is None:
            im_arr2 = np.array(im_arr_single)
        else:
            im_arr2 = np.array(im_arr2+im_arr_single)
    im_arr2 = im_arr2/n_shot
        
    slice_arr = im_arr2[int(im_arr2.shape[0]/2),:]
    slice_axs = (np.arange(len(slice_arr))-int(.5*len(slice_arr)))*17.94e-6*100
    
    if smoothBack:
        window_width = 9
        cumsum_vec = np.cumsum(np.insert(slice_arr, 0, 0)) 
        ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
        slice_arr = ma_vec
        
        slice_axs = slice_axs[int(window_width/2):-int(window_width/2)]
        
    if doNormalize:
        plt.plot(slice_axs,slice_arr*pressure_arr[-1]/pressure_arr[i],c=colors[3*i],label=str(pressure_arr[i]) + " psi")
    else:
        plt.plot(slice_axs,slice_arr,c=colors[3*i+1],ls='solid',label=str(pressure_arr[i]) + " psi")#,label=" '' With Beam = "+str(pressure_arr[i]) + " psi")

plt.xlabel("Longitudinal Position (cm)")
if doNormalize:
    plt.ylabel("Normalized Intensity (arb. units)")
else:
    plt.ylabel("Camera Intensity (arb. units)")
plt.legend();plt.show()

