#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:24:28 2022

Loop across all datasets given at a specific energy slice to see the change vs pressure

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit
import scipy.io as sio

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#emit_N = 5e-6 #Normalized Emittance m-rad
gam_L = 19569.5 #Beam Energy

def Betatron(x, *p):
    #bstar, xstar = p
    bstar, xstar, emit_N = p
    return np.sqrt(emit_N/gam_L*(bstar+np.square(x-xstar)/bstar))

gasjetpos = 1993.3
dtotr2_cal = 30.3e-6 #m/pixel
dtotr2_m = 9.8 #From a rough estimate, should calculate more rigorously

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'

#For lowres simulations
dataset_arr = np.array([
        'E308_02493',
        'E308_02512',
        'E308_02492',
        'E308_02498',
        'E308_02501',
        'E308_02506'])
pressure_arr = np.array([0,1,6,24,57.8,115.8])

beam_slice = 0.80

threshold = 20 #Threshold in DTOTR2 in which to ignore camera counts
massage_low = 0.8 
massage_high = 40 
sig_cent_guess = -25 #for when there are spikes fucking up my algo

doPlot = False
doReport = False

colors = plt.cm.brg(np.linspace(0, 1, len(dataset_arr)))
plt.figure(figsize=(10,7))
ymax = None
ymin = None

if doPlot:
    print("Can't really do this here, check the single dir version")
    sys.exit()

b_betastar = np.zeros(len(dataset_arr))
o_betastar = np.zeros(len(dataset_arr))
b_z0 = np.zeros(len(dataset_arr))
o_z0 = np.zeros(len(dataset_arr))
b_emitN = np.zeros(len(dataset_arr))
o_emitN = np.zeros(len(dataset_arr))

for m in range(len(dataset_arr)):
    print(" "+str(m/len(dataset_arr)*100)+" %")
    dataset = dataset_arr[m]
    
    ###############################################################################
    #Load Metadata from H5
    
    path = superpath + day + dataset + '/' + dataset +'.h5'
    f = h5.File(path,"r") 
    data = f['data']
    save_info = data['save_info']
    params = data['params']
        
    if doReport:
        print("Comment:")
        print(params['comment'][0][0])
        print()
        print("Scan Function: ",params['scanFuncs'][0][0])
    func_start = params['startVals'][0][0]
    func_end = params['stopVals'][0][0]
    if doReport:
        print(" ",func_start," to ",func_end)
        print()
    
    camerapath = 'DTOTR2'
    path = superpath + day + dataset + '/images/' + camerapath + '/'
    image_list = os.listdir(path)
    
    n_shot = int(params['n_shot'][0][0])
    n_step = int(len(image_list)/n_shot)
    
    if doReport:
        print("Shots per Step: ",n_shot)
        print("Num Steps: ",n_step)
        print()
    
    f.close()
    
    ###############################################################################
    #Load background image for dataset
    
    if doReport:
        print("Loading Background:")
    path = superpath + day + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(path)
    
    data = mat['data_struct']
    backs = data['backgrounds'][0][0]
    background = backs['DTOTR2'][0][0]
    background = np.transpose(background)
    
    ###############################################################################
    
    #Initialize a sigma array that is n_steps long with n_shots of data at each step
    sigma_data = np.zeros((n_step,n_shot))
    
    badcount = 0
    smallcount = 0
    largecount = 0
    bullshitcount = 0
    
    #Double for loop, going through all the steps and shots
    for i in range(n_step):
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
                sigma_data[i,j] = 0.5*massage_low
                if doReport:
                    print("Bullshit Detected:")
                    print(" at Step",step,"N",nimg)
                bullshitcount = bullshitcount + 1
            else:
                #For each Beam Slice, find the energy index and fit sigma
                slice_ind = (np.where(int_arr<beam_slice)[0])[-1]
                
                slice_arr = im_arr[:,slice_ind]
                slice_axs = np.arange(len(slice_arr))-int(.5*len(slice_arr))
                
                maxpos = np.argmax(slice_arr)
                p0 = [slice_arr[np.argmax(slice_arr)], slice_axs[np.argmax(slice_arr)], 3.]
                if (np.abs(slice_arr[maxpos] - slice_arr[maxpos-2])>0.5*slice_arr[maxpos]):
                    p0 = [slice_arr[np.argmax(slice_arr)], -25, 3.]
    
                #Detect Bullshit
                try:
                    xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)

                    if abs(xcoeff[2]) < massage_low:
                        if doReport:
                            print("Small sig:",xcoeff[2])
                            print(" at Slice",beam_slice,"index",slice_ind,"Step",step,"N",nimg)
                        if doPlot:
                            plt.plot(slice_axs,slice_arr,c='b')
                            plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                            plt.show()
                        smallcount = smallcount + 1
                    if abs(xcoeff[2]) > massage_high:
                        if doReport:
                            print("Large sig:",xcoeff[2])
                            print(" at Slice",beam_slice,"index",slice_ind,"Step",step,"N",nimg)
                        if doPlot:
                            plt.plot(slice_axs,slice_arr,c='g')
                            plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                            plt.show()
                        largecount = largecount + 1
                    sigma_data[i,j] = np.abs(xcoeff[2])
                
                except RuntimeError:
                    if doReport:
                        print("Couldn't Fit at Slice",beam_slice,"index",slice_ind,"Step",step,"N",nimg)
                    if doPlot:
                        plt.plot(slice_axs,slice_arr)
                        plt.show()
                    badcount = badcount + 1
                    sigma_data[i,j] = 0.5*massage_low
                    pass
    
    sigma_data = sigma_data * dtotr2_cal / dtotr2_m
    massage_low_um = massage_low * dtotr2_cal / dtotr2_m
    massage_high_um = massage_high * dtotr2_cal / dtotr2_m
    
    if doReport:
        print(" Dataset Loop Completed!")
        print("  ",bullshitcount,"images thrown out from high background")
        print("  ",badcount,"fits thrown out for runtime errors")
        print("  ",smallcount,"fits thrown out for too small sigmas")
        print("  ",largecount,"fits thrown out for too large sigmas")
        print()

    zob_arr = np.linspace(func_start,func_end,n_step)-gasjetpos#func_start

    sigma_mean_massaged = np.zeros((n_step))
    sigma_var_massaged = np.zeros((n_step))
    for i in range(n_step):
        sigma_mean_massaged[i]=np.average(sigma_data[i,np.where((sigma_data[i,:]>massage_low_um) & (sigma_data[i,:]<massage_high_um))[0]])
        sigma_var_massaged[i] =np.sqrt(np.var(sigma_data[i,np.where((sigma_data[i,:]>massage_low_um) & (sigma_data[i,:]<massage_high_um))[0]]))

    plt.plot(zob_arr,sigma_mean_massaged*1e6,c=colors[m],label="Pressure = "+str(pressure_arr[m])+" psi, Min = "+str(round(np.min(sigma_mean_massaged),2)))
    plt.errorbar(zob_arr,sigma_mean_massaged*1e6,yerr=sigma_var_massaged*1e6,fmt="o",c=colors[m])#,label="Variance per Step")

    if ymax is None:
        ymax = np.max(sigma_mean_massaged)
    else:
        if ymax < np.max(sigma_mean_massaged):
            ymax = np.max(sigma_mean_massaged)
    if ymin is None:
        ymin = np.min(sigma_mean_massaged)
    else:
        if ymin > np.min(sigma_mean_massaged):
            ymin = np.min(sigma_mean_massaged)
            
    #Make a fit of betatron for zob_arr vs sigma_mean_massaged and save the parameters
    doOrange = True
    try:
        fit_sigma_arr = sigma_mean_massaged#*30.3e-6/9.8
        fit_rms_arr = sigma_var_massaged#*30.3e-6/9.8
        
        pf = [0.1,0,20e-6]
        fcoeff_loop, var_matrix = curve_fit(Betatron, zob_arr, fit_sigma_arr, p0=pf)
        b_betastar[m], b_z0[m], b_emitN[m] = fcoeff_loop
        if doOrange:
            fit_sigma_arr_zoom = sigma_mean_massaged
            fit_sigma_arr_zoom = fit_sigma_arr_zoom[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
            zob_arr_zoom = zob_arr[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
            pfz = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
            fcoeffz_loop, var_matrix = curve_fit(Betatron, zob_arr_zoom, fit_sigma_arr_zoom, p0=pfz)
            o_betastar[m], o_z0[m], o_emitN[m] = fcoeffz_loop
    except RuntimeError:
        print("Couldn't betatron fit at pressure ",dataset_arr[m])


print("Full Loop Done!")

plt.title(camerapath+" Beam Size at Slice "+str(beam_slice)+"; Only " +str(massage_low_um)+"< "+r'$\sigma$'+" < "+str(massage_high_um))
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.xlabel("Imaged Distance From Gas Jet (m)")
#plt.ylabel(r'$\mathrm{\sigma}$'+" (pixels)")
plt.ylabel(r'$\mathrm{\sigma \ (\sim\mu m)}$')
plt.ylim([0.95*ymin*1e6,1.05*ymax*1e6])
plt.legend();plt.show()

plt.title(camerapath+" Beta Fit vs Pressure at Energy Slice "+str(beam_slice*100)+"%")
plt.scatter(pressure_arr,b_betastar,label="Full Fit")
if doOrange:
    plt.scatter(pressure_arr,o_betastar,label="Zoom Fit")
plt.xlabel("Pressure (psi)")
plt.ylabel("Fit Beta (m)")
plt.legend(); plt.show()

plt.title(camerapath+" Z* Fit vs Pressure at Energy Slice "+str(beam_slice*100)+"%")
plt.scatter(pressure_arr,b_z0,label="Full Fit")
if doOrange:
    plt.scatter(pressure_arr,o_z0,label="Zoom Fit")
plt.xlabel("Pressure (psi)")
plt.ylabel("Fit Beam Waist (m)")
plt.ylim([-2,1])
plt.legend(); plt.show()

plt.title(camerapath+" Emittance Fit vs Pressure at Energy Slice "+str(beam_slice*100)+"%")
plt.scatter(pressure_arr,b_emitN*1e6,label="Full Fit")
if doOrange:
    plt.scatter(pressure_arr,o_emitN*1e6,label="Zoom Fit")
plt.xlabel("Pressure (psi)")
plt.ylabel("Fit "+r'$\mathrm{\epsilon_N \ (\mu m-rad)}$')
plt.legend(); plt.show()




