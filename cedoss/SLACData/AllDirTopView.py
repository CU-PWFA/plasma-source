#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 10:31:10 2022

Loading in images of TopView

Here I want to plot lineaout vs pressure

Unfortunately, had to use 1psi for background since 0psi has no laser
At 1psi, used a nolaser background and just manually cropped out the edges
Doens't really matter, 1psi was mostly noise and no real plasma signal

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
day = '20220812/'

#For lowres simulations
dataset_arr = np.array([
        #'E308_02493',
        'E308_02484',#no beam here, later dataset at 1 psi had no topview :(
        'E308_02492',
        'E308_02498',
        'E308_02501',
        'E308_02506'])
pressure_arr = np.array([1,6,24,57.8,115.8])#0 is boring

dataset = 'E308_02484'
print("Loading Background:")
path = superpath + day + dataset + '/' + dataset +'.mat'
mat = sio.loadmat(path)

data = mat['data_struct']
backs = data['backgrounds'][0][0]
topview_bg = backs['TopView'][0][0]
back_arr = np.transpose(topview_bg)
"""
plt.set_cmap('gray')
plt.imshow(back_arr)
CB = plt.colorbar(orientation='horizontal')
plt.title("background"); plt.show()
"""

dataset = 'E308_02493'
print("Loading Background:")
path = superpath + day + dataset + '/' + dataset +'.mat'
mat = sio.loadmat(path)

data = mat['data_struct']
backs = data['backgrounds'][0][0]
topview_bg = backs['TopView'][0][0]
back_arr2 = np.transpose(topview_bg)


doNormalize = False
compareBeam = False
smoothBack = True
tryFits = True
colors = plt.cm.brg(np.linspace(0, 1, 3*len(dataset_arr)))
if compareBeam:
    plt.figure(figsize=(6.5,5))
else:
    plt.figure(figsize=(6.5,4))

for i in range(len(dataset_arr)):
    if i != 0:
        dataset = dataset_arr[i]
        path = superpath + day + dataset + '/' + dataset +'.h5'
        
        ###BACKGROUND FROM MAT FILE
        
        #print("Loading Background:")
        path = superpath + day + dataset + '/' + dataset +'.mat'
        mat = sio.loadmat(path)
        
        data = mat['data_struct']
        backs = data['backgrounds'][0][0]
        im_arr = backs['TopView'][0][0]
        im_arr = np.transpose(im_arr)
        """
        plt.set_cmap('gray')
        plt.imshow(im_arr)
        CB = plt.colorbar(orientation='horizontal')
        plt.title("image"); plt.show()
        """
        params = data['params'][0][0]
        n_shot = params['n_shot'][0][0][0][0]
        #print("Shots",n_shot)
        n_shot = 1
        #######################
        if compareBeam:
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
                
                for x in range(im_arr.shape[0]):
                    for y in range(im_arr.shape[1]):
                        im_arr_single[x][y] = im_arr_single[x][y] - back_arr[x][y]
                        if im_arr_single[x][y] > 65000:
                            im_arr_single[x][y] = 0
                if im_arr2 is None:
                    im_arr2 = np.array(im_arr_single)
                else:
                    im_arr2 = np.array(im_arr2+im_arr_single)
            im_arr2 = im_arr2/n_shot

        if i != 0:
            for x in range(im_arr.shape[0]):
                for y in range(im_arr.shape[1]):
                    im_arr[x][y] = im_arr[x][y] - back_arr[x][y]
                    if im_arr[x][y] > 1000:#65000:  #Lowered to remove the nozzle edges
                        im_arr[x][y] = 0
        else:
            for x in range(im_arr.shape[0]):
                for y in range(im_arr.shape[1]):
                    im_arr[x][y] = im_arr[x][y] - back_arr2[x][y]
                    if im_arr[x][y] > 100:#65000:  #Lowered to remove the nozzle edges
                        im_arr[x][y] = 0
        
        slice_arr = im_arr[int(im_arr.shape[0]/2),:]
        #slice_axs = (np.arange(len(slice_arr))-int(.5*len(slice_arr))+550)/530
        slice_axs = (np.arange(len(slice_arr))-int(.5*len(slice_arr)))*17.94e-6*100
        
        if compareBeam:
            slice_arr2 = im_arr2[int(im_arr2.shape[0]/2),:]    
        if smoothBack:
            window_width = 9
            cumsum_vec = np.cumsum(np.insert(slice_arr, 0, 0)) 
            ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            #padding = np.zeros(int((window_width-1)/2))+.8*slice_arr[0]
            #slice_arr = np.append(np.append(padding,ma_vec),padding)
            slice_arr = ma_vec
            
            #window_width = 9
            cumsum_vec = np.cumsum(np.insert(slice_arr2, 0, 0)) 
            ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
            #padding = np.zeros(int((window_width-1)/2))+.8*slice_arr2[0]
            #slice_arr2 = np.append(np.append(padding,ma_vec),padding)
            slice_arr2 = ma_vec
            
            slice_axs = slice_axs[int(window_width/2):-int(window_width/2)]
            
        if doNormalize:
            if i != 0:
                plt.plot(slice_axs,slice_arr*pressure_arr[-1]/pressure_arr[i],c=colors[3*i],label=str(pressure_arr[i]) + " psi")
        else:
            plt.plot(slice_axs,slice_arr,c=colors[3*i],ls='solid',label=str(pressure_arr[i]) + " psi")
            if compareBeam:
                plt.plot(slice_axs,slice_arr2,c=colors[3*i+1],ls='solid')#,label=" '' With Beam = "+str(pressure_arr[i]) + " psi")

        if tryFits:
            p0=[max(slice_arr),0,0.5]
            xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)
            plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c=colors[3*i],ls='dashed')
            
            p0=[max(slice_arr),0,0.5,2]
            scoeff, var_matrix = curve_fit(SuperGauss, slice_axs, slice_arr, p0=p0)
            plt.plot(slice_axs,SuperGauss(slice_axs,*scoeff),c=colors[3*i],ls='dotted')
            
            print("Backing Pressure: ",pressure_arr[i])
            print(" Gaussian Fit (A,mu,sgima):")
            print(xcoeff)
            print(" Super-Gaussian Fit (A,mu,sgima,n):")
            print(scoeff)
            print()

plt.xlabel("Longitudinal Position (cm)")
if doNormalize:
    plt.ylim([-50,1400])
    plt.ylabel("Normalized Intensity (arb. units)")
else:
    if compareBeam:
        plt.ylim([-50,2100])
    else:
        plt.ylim([-50,1000])
    plt.ylabel("Camera Intensity (arb. units)")
plt.legend();plt.show()

