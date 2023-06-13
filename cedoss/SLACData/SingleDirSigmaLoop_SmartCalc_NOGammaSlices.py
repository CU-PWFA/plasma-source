#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:11:23 2023

This version of the loop skips energy slices, and just grabs the sigma from the total projection.

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

def BetatronSetEmittance(x, *p):
    emit_f = 5.17831716249e-05
    bstar, xstar = p
    return np.sqrt(emit_f/gam_L*(bstar+np.square(x-xstar)/bstar))

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

dataset = 'E308_02512' #1 psig
#dataset = 'E308_02513' #1 psi, highres

#dataset = 'E308_02492' #6 psig
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
    #beam_slices = np.array([.20,.30,.40,.50,.60,.70,.80])
    #sigma_data = np.zeros((len(beam_slices),n_step,n_shot))
    
    sigma_proj_data = np.zeros((n_step,n_shot))
    
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
            
            ##Here I just want the projection of the y
            sig_axis_arr = np.zeros(im_arr.shape[0])
            sig_proj_arr = np.zeros(im_arr.shape[0])
            for x in range(im_arr.shape[0]):
                sig_axis_arr[x] = x
                sig_proj_arr[x] = np.sum(im_arr[x,:])
            #plt.plot(sig_axis_arr,sig_proj_arr);plt.show()
            
            #Take Energy Projection and Integrate
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
                sigma_proj_data[i,j] = 0.5*massage_low
                if doReport:
                    print("Bullshit Detected:")
                    print(" at Step",step,"N",nimg)
                bullshitcount = bullshitcount + 1
                
            else:
                #Now, take the projection and fit a sigma to it
                
                    maxpos = np.argmax(sig_proj_arr)
                    p0 = [sig_proj_arr[maxpos], sig_axis_arr[maxpos], 10.]
                    
                    #if (np.abs(slice_arr[maxpos] - slice_arr[maxpos-2])>0.5*slice_arr[maxpos]):
                    #    p0 = [slice_arr[np.argmax(slice_arr)], -25, 3.]
        
                    #Detect Bullshit
                    try:
                        xcoeff, var_matrix = curve_fit(Gauss, sig_axis_arr, sig_proj_arr, p0=p0)
    
                        if abs(xcoeff[2]) < massage_low:
                            if doReport:
                                print("Small sig:",xcoeff[2])
                                print(" at index",sig_proj_arr,"Step",step,"N",nimg)
                            if doPlot:
                                plt.plot(sig_axis_arr,sig_proj_arr,c='b')
                                plt.plot(sig_axis_arr,Gauss(sig_axis_arr,*xcoeff),c='r')
                                plt.show()
                            smallcount = smallcount + 1
                        if abs(xcoeff[2]) > massage_high:
                            if doReport:
                                print("Large sig:",xcoeff[2])
                                print(" at index",sig_proj_arr,"Step",step,"N",nimg)
                            if doPlot:
                                plt.plot(sig_axis_arr,sig_proj_arr,c='g')
                                plt.plot(sig_axis_arr,Gauss(sig_axis_arr,*xcoeff),c='r')
                                plt.show()
                            largecount = largecount + 1
                        sigma_proj_data[i,j] = np.abs(xcoeff[2])
                    
                    except RuntimeError:
                        if doReport:
                            print("Couldn't Fit at Step",step,"N",nimg)
                        if doPlot:
                            plt.plot(sig_axis_arr,sig_proj_arr)
                            plt.show()
                        badcount = badcount + 1
                        sigma_proj_data[i,j] = 0.5*massage_low
                        pass
                    
            """        
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
            """
    if staticMvG:
        for i in range(n_step):
            sigma_proj_data[i,:] = sigma_proj_data[i,:] * dtotr2_cal / m_arr[i]
    
    print("Loop Completed!")
    print(" ",bullshitcount,"images thrown out from high background")
    print(" ",badcount,"fits thrown out for runtime errors")
    print(" ",smallcount,"fits thrown out for too small sigmas")
    print(" ",largecount,"fits thrown out for too large sigmas")
    print()

zob_arr = np.linspace(func_start,func_end,n_step)-gasjetpos#func_start

colors = plt.cm.brg(np.linspace(1, 0, 6))
plotcolor = colors[1]

sigmafull_mean = np.zeros(n_step)
sigmafull_var = np.zeros(n_step)
for i in range(n_step):
    sigmafull_mean[i]=np.average(sigma_proj_data[i,np.where((sigma_proj_data[i,:]>massage_low*dtotr2_cal/m_arr[i]) & (sigma_proj_data[i,:]<massage_high*dtotr2_cal/m_arr[i]))[0]])
    sigmafull_var[i] =np.sqrt(np.var(sigma_proj_data[i,np.where((sigma_proj_data[i,:]>massage_low*dtotr2_cal/m_arr[i]) & (sigma_proj_data[i,:]<massage_high*dtotr2_cal/m_arr[i]))[0]]))
    sigmafull_var[i] = np.sqrt(np.square(sigmafull_var[i]) + np.square(dtotr2_cal/m_arr[i]))

plt.figure(figsize=(6,4))
plt.plot(zob_arr,sigmafull_mean*1e6,c="k",label="Full Projection")#, Min = "+str(round(np.min(sigma_mean_massaged[k]),2)))
plt.errorbar(zob_arr,sigmafull_mean*1e6,yerr=sigmafull_var*1e6,fmt="o",c="k")#,label="Variance per Step")

#plt.title(dataset + ": "+camerapath+" Beam Size at Slices; Only " +str(massage_low)+"< "+r'$\sigma$'+" < "+str(massage_high))
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")
plt.xlabel("Object Plane minus Gas Jet Location (m)")
#plt.ylabel(r'$\mathrm{\sigma}$'+" (pixels)")
plt.ylabel(r'$\mathrm{\sigma_{obj} \ (\mu m)}$')
plt.ylim([0.75*np.min(sigmafull_mean)*1e6,1.25*np.max(sigmafull_mean)*1e6])
plt.legend();plt.show()


#Doing a Betatron Fit

setEmittance = False
theoryCurve = True


print("Fitting betatron")
fit_sigma_arr = sigmafull_mean
fit_rms_arr = sigmafull_var
zob_arr_fine = np.linspace(zob_arr[0],zob_arr[-1],1000)

pf = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
fcoeff, var_matrix = curve_fit(Betatron, zob_arr, fit_sigma_arr, p0=pf)

if fitZoom:
    fit_sigma_arr_zoom = sigmafull_mean
    fit_sigma_arr_zoom = fit_sigma_arr_zoom[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
    zob_arr_zoom = zob_arr[np.argmin(fit_sigma_arr)-5:np.argmin(fit_sigma_arr)+5]
    pfz = [0.08,zob_arr[np.argmin(fit_sigma_arr)],5e-3]
    fcoeffz, var_matrix = curve_fit(Betatron, zob_arr_zoom, fit_sigma_arr_zoom, p0=pfz)

min_rms=fit_rms_arr[np.argmin(fit_sigma_arr)]*1e6

plt.figure(figsize=(8,5))
#plt.title(dataset + ": "+camerapath+" Beam Size at Slices Vs Betatron Fit")
plt.semilogy(zob_arr,fit_sigma_arr*1e6,c=plotcolor,label="Min = "+str(round(np.min(fit_sigma_arr)*1e6,2))+r'$\mu m \ \pm$'+str(np.round(min_rms,2))+r'$\mu m$')
plt.errorbar(zob_arr,fit_sigma_arr*1e6,yerr=fit_rms_arr*1e6,fmt="o",c=plotcolor)#,label="Variance per Step")
plt.plot(zob_arr_fine,Betatron(zob_arr_fine,*fcoeff)*1e6,c='b',ls='--',label="Betatron Fit to Full Curve")
if fitZoom:
    plt.plot(zob_arr_fine,Betatron(zob_arr_fine,*fcoeffz)*1e6,c='orange',ls='--',label="Betatron Fit to +/- 5 of Min.")
#plt.xlabel("Spec z_obj + "+str(func_start)+" (m)")

if setEmittance:
    pf = [0.08,zob_arr[np.argmin(fit_sigma_arr)]]
    scoeff, var_matrix = curve_fit(BetatronSetEmittance, zob_arr, fit_sigma_arr, p0=pf)
    plt.plot(zob_arr_fine,BetatronSetEmittance(zob_arr_fine,*scoeff)*1e6,c='r',ls='--',label="Betatron Fit with Theoretical Emittance")

if theoryCurve:
    plt.plot(zob_arr_fine,Betatron(zob_arr_fine,*[0.0012039, 0.037, 5.1783e-5])*1e6,c='black',ls='--',label="Betatron with Theoretical Values")
    plt.plot(zob_arr,dtotr2_cal/m_arr*1e6,c='r',ls='dotted',label="Pixel Resolution Limit of DTOTR2")
    plt.ylim([1.6,200])
    

plt.xlabel("Imaged Distance From Gas Jet (m)")
plt.ylabel(r'$\mathrm{\sigma \ (\mu m)}$')
plt.legend();plt.show()

print("Blue Curve:")
print(" BetaStar = ",fcoeff[0]*100,"cm")
print(" Z_0 = ",fcoeff[1],"m")
print(" EmitN = ",fcoeff[2]*1e6,"mm-mrad")
if fitZoom:
    print()
    print("Orange Curve:")
    print(" BetaStar = ",fcoeffz[0]*100,"cm")
    print(" Z_0 = ",fcoeffz[1],"m")
    print(" EmitN = ",fcoeffz[2]*1e6,"mm-mrad")
    
if setEmittance:
    print()
    print("Red Curve:")
    print(" BetaStar = ",scoeff[0]*100,"cm")
    print(" Z_0 = ",scoeff[1],"m")
