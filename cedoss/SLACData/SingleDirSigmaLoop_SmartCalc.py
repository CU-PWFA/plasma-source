#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 12:08:32 2022


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

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#emit_N = 5e-6 #Normalized Emittance m-rad
gam_L = 19569.5 #Beam Energy

def Betatron(x, *p):
    #bstar, xstar = p
    bstar, xstar, emit_N = p
    return np.sqrt(emit_N/gam_L*(bstar+np.square(x-xstar)/bstar))

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

#dataset = 'E308_02493' #no gas

#dataset = 'E308_02512' #1 psig
#dataset = 'E308_02513' #1 psi, highres

dataset = 'E308_02492' #6 psig
#dataset = 'E308_02497' #6 psi, highres

#dataset = 'E308_02498' #24 psi
#dataset = 'E308_02500' #24 psi, highres

#dataset = 'E308_02501' #5 bar over 0, 57.8 psi
#dataset = 'E308_02505' #5 bar over 0, 57.8 psi, highres

#dataset = 'E308_02506' #9 bar over 0, 115.8 psi
#dataset = 'E308_02511' #9 bar over 0, 115.8 psi, highres, can't do 0.2 or 0.9

doLoop = True
fitZoom = False
staticMvG = True #Set to true to calc M using 10Gev, False to do relative position to beam centroid
fit_ind = 3

if doLoop:
    
    #cut_pos = 130
    #y_roi = [40,200]
    #y_roi = [0,267]
    
    threshold = 20 #Threshold in DTOTR2 in which to ignore camera counts
    massage_low = 0.8
    massage_high = 40
    sig_cent_guess = -25 #for when there are spikes fucking up my algo
    
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
    #Calculate the magnification for each step
    if staticMvG:
        gamma = 19569.5 #10 GeV
        m_arr = np.abs(MatrixCalc.getM11ForAllSteps(superpath,day,dataset,gamma))
        print(m_arr)
    #m_arr = m_arr * 2 #JUST TO GET SIGMAS OF WIRESCANNER HERE.  WHY THO
    
    #Initialize a sigma array that is n_steps long with n_shots of data at each step
    #beam_slices = np.array([0.1,.20,.30,.40,.50,.60,.70,.80,.90])
    beam_slices = np.array([.20,.30,.40,.50,.60,.70,.80])
    sigma_data = np.zeros((len(beam_slices),n_step,n_shot))
    
    badcount = 0
    smallcount = 0
    largecount = 0
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
                for k in range(len(beam_slices)):
                    sigma_data[k,i,j] = 0.5*massage_low
                if doReport:
                    print("Bullshit Detected:")
                    print(" at Step",step,"N",nimg)
                bullshitcount = bullshitcount + 1
            else:
                #For each Beam Slice, find the energy index and fit sigma
                for k in range(len(beam_slices)):
                    slice_ind = (np.where(int_arr<beam_slices[k])[0])[-1]
                    
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
                                print(" at Slice",beam_slices[k],"index",slice_ind,"Step",step,"N",nimg)
                            if doPlot:
                                plt.plot(slice_axs,slice_arr,c='b')
                                plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                                plt.show()
                            smallcount = smallcount + 1
                        if abs(xcoeff[2]) > massage_high:
                            if doReport:
                                print("Large sig:",xcoeff[2])
                                print(" at Slice",beam_slices[k],"index",slice_ind,"Step",step,"N",nimg)
                            if doPlot:
                                plt.plot(slice_axs,slice_arr,c='g')
                                plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                                plt.show()
                            largecount = largecount + 1
                        sigma_data[k,i,j] = np.abs(xcoeff[2])
                    
                    except RuntimeError:
                        if doReport:
                            print("Couldn't Fit at Slice",beam_slices[k],"Step",step,"N",nimg)
                        if doPlot:
                            plt.plot(slice_axs,slice_arr)
                            plt.show()
                        badcount = badcount + 1
                        sigma_data[k,i,j] = 0.5*massage_low
                        pass
                    
            # Now that we are done looping over beam slices, calc m for each slice
            if not staticMvG:
                center_ind = (np.where(int_arr<beam_slices[3])[0])[-1]
                E_cent = 10 #GeV
                ybeam = ECaliInverse(center_ind*dx,E_cent)
                for k in range(len(beam_slices)):
                    slice_ind = (np.where(int_arr<beam_slices[k])[0])[-1]
                    gam_slice = ECali(slice_ind*dx,ybeam) * gam_L / E_cent
                    m_arr = np.abs(MatrixCalc.getM11ForAllSteps(superpath,day,dataset,gam_slice))
                    sigma_data[k,i,j] = sigma_data[k,i,j] * dtotr2_cal / m_arr[i]
    if staticMvG:
        for i in range(n_step):
            sigma_data[:,i,:] = sigma_data[:,i,:] * dtotr2_cal / m_arr[i]
    
    print("Loop Completed!")
    print(" ",bullshitcount,"images thrown out from high background")
    print(" ",badcount,"fits thrown out for runtime errors")
    print(" ",smallcount,"fits thrown out for too small sigmas")
    print(" ",largecount,"fits thrown out for too large sigmas")
    print()

zob_arr = np.linspace(func_start,func_end,n_step)-gasjetpos#func_start

sigma_mean = np.zeros((len(beam_slices),n_step))
sigma_var = np.zeros((len(beam_slices),n_step))
for k in range(len(beam_slices)):
    for i in range(n_step):
        sigma_mean[k,i]=np.average(sigma_data[k,i,:])
        sigma_var[k,i] =np.sqrt(np.var(sigma_data[k,i,:]))

colors = plt.cm.brg(np.linspace(0, 1, len(beam_slices)))
"""
plt.figure(figsize=(10,7))
for k in range(len(beam_slices)):
    plt.plot(zob_arr,sigma_mean[k]*1e6,c=colors[k],label="E slice = "+str(beam_slices[k]*100)+"%, Min = "+str(round(np.min(sigma_mean[k]),2)))
    plt.errorbar(zob_arr,sigma_mean[k]*1e6,yerr=sigma_var[k]*1e6,fmt="o",c=colors[k])#,label="Variance per Step")
plt.title(dataset + ": "+camerapath+" Beam Size at Slices")
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.xlabel("Imaged Distance From Gas Jet (m)")
#plt.ylabel(r'$\mathrm{\sigma}$'+" (pixels)")
plt.ylabel(r'$\mathrm{\sigma \ (\sim\mu m)}$')
plt.ylim([0.95*np.min(sigma_mean)*1e6,1.05*np.max(sigma_mean)*1e6])
plt.legend();plt.show()
"""
sigma_mean_massaged = np.zeros((len(beam_slices),n_step))
sigma_var_massaged = np.zeros((len(beam_slices),n_step))
for k in range(len(beam_slices)):
    for i in range(n_step):
        sigma_mean_massaged[k,i]=np.average(sigma_data[k,i,np.where((sigma_data[k,i,:]>massage_low*dtotr2_cal/m_arr[i]) & (sigma_data[k,i,:]<massage_high*dtotr2_cal/m_arr[i]))[0]])
        sigma_var_massaged[k,i] =np.sqrt(np.var(sigma_data[k,i,np.where((sigma_data[k,i,:]>massage_low*dtotr2_cal/m_arr[i]) & (sigma_data[k,i,:]<massage_high*dtotr2_cal/m_arr[i]))[0]]))

plt.figure(figsize=(10,7))
for k in range(len(beam_slices)):
    plt.plot(zob_arr,sigma_mean_massaged[k]*1e6,c=colors[k],label="Proj. Slice = "+str(beam_slices[k]*100)+"%")#, Min = "+str(round(np.min(sigma_mean_massaged[k]),2)))
    plt.errorbar(zob_arr,sigma_mean_massaged[k]*1e6,yerr=sigma_var_massaged[k]*1e6,fmt="o",c=colors[k])#,label="Variance per Step")
#plt.title(dataset + ": "+camerapath+" Beam Size at Slices; Only " +str(massage_low)+"< "+r'$\sigma$'+" < "+str(massage_high))
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.xlabel("Object Plane minus Gas Jet Location (m)")
#plt.ylabel(r'$\mathrm{\sigma}$'+" (pixels)")
plt.ylabel(r'$\mathrm{\sigma \ (\mu m)}$')
plt.ylim([0.95*np.min(sigma_mean_massaged)*1e6,1.05*np.max(sigma_mean_massaged)*1e6])
plt.legend();plt.show()


#Doing a Single Betatron Fit


print("Fitting betatron at Energy Slice",beam_slices[fit_ind]*100,"%")
fit_sigma_arr = sigma_mean_massaged[fit_ind]#*30.3e-6/9.8
fit_rms_arr = sigma_var_massaged[fit_ind]#*30.3e-6/9.8
zob_arr_fine = np.linspace(zob_arr[0],zob_arr[-1],100)

pf = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
fcoeff, var_matrix = curve_fit(Betatron, zob_arr, fit_sigma_arr, p0=pf)

if fitZoom:
    fit_sigma_arr_zoom = sigma_mean_massaged[fit_ind]#*30.3e-6/9.8
    fit_sigma_arr_zoom = fit_sigma_arr_zoom[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
    zob_arr_zoom = zob_arr[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
    pfz = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
    fcoeffz, var_matrix = curve_fit(Betatron, zob_arr_zoom, fit_sigma_arr_zoom, p0=pfz)

plt.figure(figsize=(8,5))
#plt.title(dataset + ": "+camerapath+" Beam Size at Slices Vs Betatron Fit")
plt.plot(zob_arr,fit_sigma_arr*1e6,c=colors[fit_ind],label="E slice = "+str(beam_slices[fit_ind]*100)+"%, Min = "+str(round(np.min(fit_sigma_arr)*1e6,2))+r'$\mu m$')
plt.errorbar(zob_arr,fit_sigma_arr*1e6,yerr=fit_rms_arr*1e6,fmt="o",c=colors[fit_ind])#,label="Variance per Step")
plt.plot(zob_arr_fine,Betatron(zob_arr_fine,*fcoeff)*1e6,c='b',ls='--',label="Betatron Fit to Full Curve")
if fitZoom:
    plt.plot(zob_arr_fine,Betatron(zob_arr_fine,*fcoeffz)*1e6,c='r',ls='--',label="Betatron Fit to +/- 5 of Min.")
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.xlabel("Imaged Distance From Gas Jet (m)")
plt.ylabel(r'$\mathrm{\sigma \ (\mu m)}$')
plt.legend();plt.show()

print("Blue Curve:")
print(" BetaStar = ",fcoeff[0]*100,"cm")
print(" Z_0 = ",fcoeff[1],"m")
print(" EmitN = ",fcoeff[2]*1e6,"mm-mrad")
if fitZoom:
    print()
    print("Red Curve:")
    print(" BetaStar = ",fcoeffz[0]*100,"cm")
    print(" Z_0 = ",fcoeffz[1],"m")
    print(" EmitN = ",fcoeffz[2]*1e6,"mm-mrad")

sys.exit()

#Looping Over All Slices, and Plotting the Fit Parameters

doOrange = True

b_betastar = np.zeros(len(beam_slices))
o_betastar = np.zeros(len(beam_slices))
b_z0 = np.zeros(len(beam_slices))
o_z0 = np.zeros(len(beam_slices))
b_emitN = np.zeros(len(beam_slices))
o_emitN = np.zeros(len(beam_slices))
for i in range(len(beam_slices)):
    try:
        fit_ind = i
        fit_sigma_arr = sigma_mean_massaged[fit_ind]#*30.3e-6/9.8
        fit_rms_arr = sigma_var_massaged[fit_ind]#*30.3e-6/9.8
        zob_arr_fine = np.linspace(zob_arr[0],zob_arr[-1],100)
        
        pf = [0.1,0,20e-6]
        fcoeff_loop, var_matrix = curve_fit(Betatron, zob_arr, fit_sigma_arr, p0=pf)
        b_betastar[i], b_z0[i], b_emitN[i] = fcoeff_loop
        if doOrange:
            fit_sigma_arr_zoom = sigma_mean_massaged[fit_ind]#*30.3e-6/9.8
            fit_sigma_arr_zoom = fit_sigma_arr_zoom[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
            zob_arr_zoom = zob_arr[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
            pfz = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
            fcoeffz_loop, var_matrix = curve_fit(Betatron, zob_arr_zoom, fit_sigma_arr_zoom, p0=pfz)
            o_betastar[i], o_z0[i], o_emitN[i] = fcoeffz_loop
    except RuntimeError:
        print("Couldn't fit at slice ",beam_slices[i])
    
plt.title(dataset + ": "+camerapath+" Beta Fit vs Energy Slice")
plt.scatter(beam_slices*100,b_betastar,label="Full Fit")
if doOrange:
    plt.scatter(beam_slices*100,o_betastar,label="Zoom Fit")
plt.xlabel("Beam Energy Slices (%)")
plt.ylabel("Fit Beta (m)")
plt.legend(); plt.show()

plt.title(dataset + ": "+camerapath+" Z* Fit vs Energy Slice")
plt.scatter(beam_slices*100,b_z0,label="Full Fit")
if doOrange:
    plt.scatter(beam_slices*100,o_z0,label="Zoom Fit")
plt.xlabel("Beam Energy Slices (%)")
plt.ylabel("Fit Beam Waist (m)")
plt.ylim([-2,1])
plt.legend(); plt.show()

plt.title(dataset + ": "+camerapath+" Emittance Fit vs Energy Slice")
plt.scatter(beam_slices*100,b_emitN*1e6,label="Full Fit")
if doOrange:
    plt.scatter(beam_slices*100,o_emitN*1e6,label="Zoom Fit")
plt.xlabel("Beam Energy Slices (%)")
plt.ylabel("Fit "+r'$\mathrm{\epsilon_N \ (\mu m-rad)}$')
plt.legend(); plt.show()



