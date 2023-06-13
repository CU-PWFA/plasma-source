#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 19:59:12 2023

Loop over the dtotr2 images, and find what percentages of the beam are in
what energy bins.  then, over the full range, make a contour plot of how
this distriubiton changes when imaged.

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import scipy.io as sio
import MatrixCalc
import FACETPlotter
from functools import partial

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#emit_N = 5e-6 #Normalized Emittance m-rad
gam_L = 19569.5 #Beam Energy

def Betatron(x, *p, gam):
    #bstar, xstar = p
    bstar, xstar, emit_N = p
    return np.sqrt(emit_N/gam*(bstar+np.square(x-xstar)/bstar))

dx = 30.3e-6*1e3
#ybeam = 3.1
dnom = 60
ebend = 10

def ECaliInverse(yscreen,E_cent):
    return yscreen + dnom*(1 - ebend/E_cent)

def ECali(yscreen,ybeam):
    dy = dnom-ybeam
    return dnom*ebend/(yscreen + dy)

gasjetpos = 1993.3
dtotr2_cal = 30.3e-6 #m/pixel
#dtotr2_m = 9.8 #From a rough estimate, should calculate more rigorously

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'

#dataset = 'E308_02493';tset=0 #no gas

dataset = 'E308_02512';tset=1 #1 psig
#dataset = 'E308_02513';tset=1 #1 psi, highres

#dataset = 'E308_02492';tset=6 #6 psig
#dataset = 'E308_02497';tset=6 #6 psi, highres

#dataset = 'E308_02498';tset=24 #24 psi
#dataset = 'E308_02500';tset=24 #24 psi, highres

#dataset = 'E308_02501';tset=57.8 #5 bar over 0, 57.8 psi
#dataset = 'E308_02505';tset=57.8 #5 bar over 0, 57.8 psi, highres

#dataset = 'E308_02506';tset=115.8 #9 bar over 0, 115.8 psi
#dataset = 'E308_02511';tset=115.8 #9 bar over 0, 115.8 psi, highres, can't do 0.2 or 0.9

doLoop = True

if doLoop:
    
    #cut_pos = 130
    #y_roi = [40,200]
    #y_roi = [0,267]
    
    threshold = 20 #Threshold in DTOTR2 in which to ignore camera counts
    
    doPlot = False
    doReport = False
    
    ###############################################################################
    #Load Metadata from H5
    
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
    
    z_arr = FACETPlotter.getZarr(superpath,day,dataset)
    
    ###############################################################################
    #Load background image for dataset
    
    print("Loading Background:")
    path = superpath + day + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(path)
    
    data = mat['data_struct']
    backs = data['backgrounds'][0][0]
    background = backs['DTOTR2'][0][0]
    background = np.transpose(background)
    
    ###############################################################################
    #m_arr = m_arr * 2 #JUST TO GET SIGMAS OF WIRESCANNER HERE.  WHY THO
    
    #Initialize a sigma array that is n_steps long with n_shots of data at each step
    
    energy_slices = np.array([9.85, 9.9, 9.95, 10, 10.05, 10.1, 10.15])
    energy_slices_low = energy_slices-(0.05/2)
    energy_slices_high = energy_slices+(0.05/2)
    #energy_slices = np.array([9.6, 9.65, 9.7, 9.75, 9.8, 9.85])
    #energy_slices = np.array([9.9, 10, 10.1])
    
    percentage_data = np.zeros((len(energy_slices)+2,n_step,n_shot))
    
    badcount = 0
    bullshitcount = 0
    
    #Double for loop, going through all the steps and shots
    print("Doing the Loop:")
    for i in range(n_step):
        print(" "+str(i/n_step*100)+" %")
        for j in range(n_shot):
            #Load Image
            step = i+1
            nimg = j
            
            if nimg < 10:
                n_suffix = "_000" + str(nimg)
            else:
                n_suffix = "_00" + str(nimg)
            if step < 10:
                image_path = camerapath + "_data_step0" + str(step) + n_suffix + ".tif"
            else:
                image_path = camerapath + "_data_step" + str(step) + n_suffix + ".tif"
            
            path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
            image = plt.imread(path)
            im_arr = np.array(image)
            
            #Apply Background and Threshold
            for x in range(im_arr.shape[0]):
                for y in range(im_arr.shape[1]):
                    im_arr[x][y] = im_arr[x][y] - background[x][y]
                    if im_arr[x][y] > 65000:
                        im_arr[x][y] = 0
                    if im_arr[x][y] < threshold:
                        im_arr[x][y] = 0
            
            #Take Projection and Integrate
            en_arr = np.zeros(im_arr.shape[1])
            proj_arr = np.zeros(im_arr.shape[1])
            for x in range(im_arr.shape[1]):
                en_arr[x] = x
                proj_arr[x] = np.sum(im_arr[:,x])
            int_arr = np.zeros(len(proj_arr))
            for x in range(len(proj_arr)-1):
                int_arr[x+1] = proj_arr[x+1] + int_arr[x]
            total_counts = int_arr[-1]
            int_arr = int_arr/int_arr[-1]
            
            #If this is a BS dataset (see 02506_step14_n1 as an example) then discard
            test_ind = (np.where(int_arr<0.1)[0])[-1]
            if test_ind < 25:
                for k in range(len(energy_slices)):
                    percentage_data[k,i,j] = -100
                if doReport:
                    print("Bullshit Detected:")
                    print(" at Step",step,"N",nimg)
                bullshitcount = bullshitcount + 1
            else:
                #For each Beam Slice, find the energy index and fit sigma
                for k_tot in range(len(energy_slices)+2):
                    k=k_tot-1
                        
                    #Get an array of pixels and an array of energy
                    pixel_arr = np.arange(0,im_arr.shape[1]-1)
                    zval = z_arr[i]
                    ybeam = FACETPlotter.DTOTR2_ybeamFunc(zval)
                    energy_arr = FACETPlotter.DTOTR2_ECali(pixel_arr*dx, ybeam)
                    
                    #The pixel slice_ind is where we are just at the energy slice 
                    #slice_ind = (np.where(energy_arr>energy_slices[k])[0])[-1]
                    if k==-1:
                        highend = (np.where(energy_arr>energy_slices_low[0])[0])[-1]
                        lowend = im_arr.shape[1]-1
                    elif k==len(energy_slices):
                        highend = 0
                        lowend = (np.where(energy_arr>energy_slices_high[-1])[0])[-1]
                    else:
                        lowend = (np.where(energy_arr>energy_slices_low[k])[0])[-1]
                        highend = (np.where(energy_arr>energy_slices_high[k])[0])[-1]
                        
                    slice_arr = np.zeros(np.shape(im_arr)[0])
                    for im_i in range(np.shape(im_arr)[0]):
                        slice_arr[im_i] = np.sum(im_arr[im_i,highend:lowend])
                    binsum = np.sum(slice_arr)
                    frac = binsum/total_counts
                    
                    percentage_data[k+1,i,j] = frac

    print("Loop Completed!")
    print(" ",bullshitcount,"images thrown out from high background")
    print()

zob_arr = np.linspace(func_start,func_end,n_step)-gasjetpos#func_start

num_regions = len(energy_slices)+2
colors = plt.cm.brg(np.linspace(1, 0, num_regions))
frac_mean = np.zeros((num_regions,n_step))
frac_var = np.zeros((num_regions,n_step))

for k in range(num_regions):
    for i in range(n_step):
        frac_mean[k,i]=np.average(percentage_data[k,i,np.where((percentage_data[k,i,:]>0))[0]])
        frac_var[k,i] =np.sqrt(np.var(percentage_data[k,i,np.where((percentage_data[k,i,:]>0))[0]]))
        if k != 0:
            frac_mean[k,i]=frac_mean[k,i]+frac_mean[k-1,i]
       
plt.figure(figsize=(10*.85,7*.7))
        
for k in range(num_regions):
    if k != num_regions-1:
        plt.errorbar(zob_arr,frac_mean[k],yerr=np.sqrt(np.square(frac_var[k])+np.square(frac_var[k+1])),c="black")
    if k == 0:
        plt.fill_between(zob_arr, y1=frac_mean[k], y2=0, color = colors[k],label="<9.83")
    elif k == num_regions-1:
        plt.fill_between(zob_arr, y1=1.0, y2=frac_mean[k-1], color = colors[k],label=">10.18")
    else:
        plt.fill_between(zob_arr, y1=frac_mean[k], y2=frac_mean[k-1], color = colors[k], label="%.2f"%energy_slices[k-1])#+r'$\pm 0.025$')

plt.ylim([0,1])
plt.xlim([zob_arr[0],zob_arr[-1]])
plt.xlabel("Object Plane minus Gas Jet Location (m)")
plt.ylabel("Cumulative Beam Fraction")
plt.legend(title="Bin Energy (GeV)",framealpha=0.5,loc=1)







