#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 11:15:38 2022

Module that hosts scripts to generate good looking plots of FACET data

@author: chris
"""
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import os as os

#dtotr2_dx = 30.3e-6*1e3
#ybeam = 3.1

#This was because in my shift the spectrometer wasn't aligned correctly,
# and when there was no plasma and the object plane was scanned the
# energy centroid of the beam would shift.
#Using DTOTR2_ybeamcalibration.py for the values at 10 GeV
def DTOTR2_ybeamFunc(zval):
    a0 = 1022.04328885
    a1 = -0.51128705777
    return zval*a1 + a0

#Requires the usual 3 strings for the dataset
#Returns a 1d array for the z position of each step
def getZarr(superpath, day, dataset):
    path = superpath + day + '/' + dataset + '/' + dataset +'.h5'
    
    f = h5.File(path,"r") 
    data = f['data']
    params = data['params']
    func_start = params['startVals'][0][0]
    func_end = params['stopVals'][0][0]
    camerapath = 'DTOTR2'
    path = superpath + day + '/' + dataset + '/images/' + camerapath + '/'
    image_list = os.listdir(path)
    n_shot = int(params['n_shot'][0][0])
    n_step = int(len(image_list)/n_shot)
    f.close()
    return np.linspace(func_start,func_end,n_step)

#Requires the y axis value (yscreen) and the ybeam from DTOTR2_ybeamFunc
#Other inputs are constants from Doug S.
#Returns the approximate Energy in GeV
def DTOTR2_ECali(yscreen, ybeam, dtotr2_dnom = 60, dtotr2_ebend = 10):
    dy = dtotr2_dnom-ybeam
    return dtotr2_dnom*dtotr2_ebend/(yscreen + dy)

#Requires the usual 3 strings for the dataset, as well as the specific step and shot
#Threshold is an optional parameter for cutting out low noise
#Returns the background-subtracted 2d array for the image, as well as the pixel-um calibration
def loadDTOTR2(superpath, day, dataset, step, nimg, threshold = 20):
    mat_path = superpath + day + '/' + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(mat_path)
    
    data = mat['data_struct']
    backs = data['backgrounds'][0][0]
    background = backs['DTOTR2'][0][0]
    background = np.transpose(background)

    camerapath = 'DTOTR2'
    if step < 10:
        image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
    else:
        image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
    
    path = superpath + day + '/' + dataset + '/images/' + camerapath + '/' + image_path
    image = plt.imread(path)
    im_arr = np.array(image)
    
    for i in range(im_arr.shape[0]):
        for j in range(im_arr.shape[1]):
            im_arr[i][j] = im_arr[i][j] - background[i][j]
            if im_arr[i][j] > 65000:
                im_arr[i][j] = 0
            if im_arr[i][j] < threshold:
                im_arr[i][j] = 0
    
    metadata = data['metadata'][0][0]
    metadtotr2 = metadata['DTOTR2'][0][0]
    pixelcal = metadtotr2['RESOLUTION'][0][0][0][0]       
    
    return im_arr, pixelcal

#Requires the 2d image array, the pixel-to-mm calibration, and ybeam from DTOTR2_ybeamcalibration
def plotDTOTR2_withProjection(im_arr, dx, ybeam, annotation="Test"):
    hori_range = np.array([0,im_arr.shape[1]])
    vert_range = (np.array([0,im_arr.shape[0]])-1/2*im_arr.shape[0])*dx
    
    hori_arr = DTOTR2_ECali(np.arange(0,im_arr.shape[1])*dx, ybeam)
    aspect_ratio = 0.5*1/np.abs((vert_range[1]-vert_range[0])/(hori_range[1]-hori_range[0]))
    
    plt.set_cmap('OrRd')
    fig, ax = plt.subplots(1,1)
    ax.imshow(np.log(im_arr)*150,aspect=aspect_ratio,vmin=450,vmax=1000,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    #ax.imshow(im_arr,aspect=aspect_ratio,vmin=0,vmax=1000,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    ax.set_xlabel("Energy (GeV)")
    ax.set_ylabel("Position on Screen (mm)")

    #energy_ticks = np.arange(9.2,10.61,0.2)
    energy_ticks = np.arange(1.0,20.1,0.2)
    energy_ticks = energy_ticks[np.where(energy_ticks<np.max(hori_arr))[0]]
    energy_ticks = energy_ticks[np.where(energy_ticks>np.min(hori_arr))[0]]

    pixel_ticks = np.zeros(len(energy_ticks))
    for i in range(len(hori_arr)):
        for j in range(len(pixel_ticks)):
            if hori_arr[i] > energy_ticks[j]:
                pixel_ticks[j] = int(i)

    label_ticks = np.zeros(len(energy_ticks))
    for i in range(len(energy_ticks)):
        label_ticks[i] = str(energy_ticks[i])
    
    ax.set_xticks(pixel_ticks)
    ax.set_xticklabels(label_ticks)
    
    en_arr = np.linspace(hori_range[0],hori_range[1]-1,im_arr.shape[1])
    proj_arr = np.zeros(im_arr.shape[1])
    for i in range(im_arr.shape[1]):
        proj_arr[i] = (np.sum(im_arr[:,i])*(1/7000))+vert_range[0]
        
    plt.plot(en_arr,proj_arr)
    
    s_arr = np.linspace(vert_range[0],vert_range[1],im_arr.shape[0])
    proj_arr2 = np.zeros(im_arr.shape[0])
    for i in range(im_arr.shape[0]):
        proj_arr2[i] = (np.sum(im_arr[i,:])*(1/700))+hori_range[0]+1.5
    plt.plot(np.flip(proj_arr2,0),s_arr)
    
    plt.text(170,2.3,annotation,bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
    plt.grid()
    plt.show()
    return

#Requires the 2d image array, the pixel-to-mm calibration, and ybeam from DTOTR2_ybeamcalibration
def plotDTOTR2_rawData(im_arr, dx, ybeam, annotation="Test"):
    im_arr = np.transpose(im_arr)
    
    hori_range = np.array([0,im_arr.shape[1]])*dx
    vert_range = np.array([0,im_arr.shape[0]])*dx
    
    plt.set_cmap('gray')
    fig, ax = plt.subplots(1,1)
    ax.imshow(im_arr,aspect=1,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    #ax.imshow(im_arr,aspect=aspect_ratio,vmin=0,vmax=1000,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    ax.set_xlabel("Screen x (mm)")
    ax.set_ylabel("Screen y (mm)")

    plt.show()
    return

if __name__ == '__main__':
    #testing:
    tsuperpath = '/media/chris/New Volume/SLACData/'
    tday = '20220812'
    tdataset = 'E308_02493';tpsi=0 #no gas
    tdataset = 'E308_02512';tpsi=1 #1 psig
    #tdataset = 'E308_02513';tpsi=1 #1psig with high res scan
    #tdataset = 'E308_02492';tpsi=6 #6 psig
    #tdataset = 'E308_02498';tpsi=24 #24 psi
    #tdataset = 'E308_02501';tpsi=57.8 #5 bar over 0, 57.8 psi
    #tdataset = 'E308_02506';tpsi=115.8 #9 bar over 0, 115.8 psi
    
    tstep = 1#+tstep
    tnimg = 1# + tnimg
    
    im_arr, pixelcal = loadDTOTR2(tsuperpath, tday, tdataset, tstep, tnimg)
    
    pixelcal = pixelcal * 1e-6*1e3 #to get to mm
    z_arr = getZarr(tsuperpath, tday, tdataset)
    ybeam = DTOTR2_ybeamFunc(z_arr[tstep-1])
    
    gasjet = 1993.27
    annotation = "Backing Pressure: "+str(tpsi)+" psi"+'\n'+"Object Plane: "+str(z_arr[tstep-1]-1993.27)+" m"
    annotation = "Backing Pressure: "+str(tpsi)+" psi"+'\n'+"Object Plane: "+("%.2f" % (z_arr[tstep-1]-1993.27))+" m"
    
    
    plotDTOTR2_withProjection(im_arr, pixelcal, ybeam, annotation)
    plotDTOTR2_rawData(im_arr, pixelcal, ybeam, annotation)
    #print(DTOTR2_ECali(1.0, ybeam))
