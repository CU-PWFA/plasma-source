#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:07:44 2022

Looking at the 0 psi dataset, find ybeam vs z.  Also calculate the error in gamma

Loops over all datasets

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
#from scipy import optimize
from scipy.optimize import curve_fit
import scipy.io as sio
import matplotlib.colors as colors
import FACETPlotter

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

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'

dataset_arr = np.array(['E308_02493','E308_02512','E308_02492','E308_02498','E308_02501','E308_02506'])
label_arr =   np.array(['0 psi',     '1 psi',     '6 psi',     '24 psi',    '57.8 psi',  '115.8 psi',])

dx = 30.3e-6*1e3
gasjetpos = 1993.3

doReport = True
bullshitcount = 0

#Doing it the lazy way, either plot the ybeam or plot the charge, not both
choice = 1 # 0 for ybeam, 1 for qbeam, 2 for ebeam, 3 for both ybeam and ebeam

if choice == 3:
    fig, ax1 = plt.subplots(figsize=(11.5,3.5),nrows=1, ncols=2)

for d in range(len(dataset_arr)):
    dataset = dataset_arr[d]
    path = superpath + day + dataset + '/' + dataset +'.h5'
    f = h5.File(path,"r") 
    data = f['data']
    save_info = data['save_info']
    params = data['params']
    func_start = params['startVals'][0][0]
    func_end = params['stopVals'][0][0]
    camerapath = 'DTOTR2'
    path = superpath + day + dataset + '/images/' + camerapath + '/'
    image_list = os.listdir(path)
    
    n_shot = int(params['n_shot'][0][0])
    n_step = int(len(image_list)/n_shot)
    f.close()
    #n_shot=1
    ybeam_arr = np.zeros((n_step,n_shot))
    qbeam_arr = np.zeros((n_step,n_shot))
    Ebeam_arr = np.zeros((n_step,n_shot))
    #step = 1
    for s in range(n_step):
        step = s+1
        for n in range(n_shot):
            nimg = n
            
            path = superpath + day + dataset + '/' + dataset +'.h5'
            ###BACKGROUND FROM MAT FILE
            
            path = superpath + day + dataset + '/' + dataset +'.mat'
            mat = sio.loadmat(path)
            
            data = mat['data_struct']
            backs = data['backgrounds'][0][0]
            dtotr2_back = backs['DTOTR2'][0][0]
            dtotr2_back = np.transpose(dtotr2_back)
            #######################
            camerapath = 'DTOTR2'
            if step < 10:
                image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
            else:
                image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
            
            path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
            
            image = plt.imread(path)
            
            im_arr = np.array(image)
            
            threshold = 0
            for i in range(im_arr.shape[0]):
                for j in range(im_arr.shape[1]):
                    im_arr[i][j] = im_arr[i][j] - dtotr2_back[i][j]
                    if im_arr[i][j] > 65000:
                        im_arr[i][j] = 0
                    if im_arr[i][j] < threshold:
                        im_arr[i][j] = 0
            
            hori_range = np.array([0,im_arr.shape[1]])*dx
            en_arr = np.linspace(hori_range[0],hori_range[1],im_arr.shape[1])
            proj_arr = np.zeros(im_arr.shape[1])
            for i in range(im_arr.shape[1]):
                proj_arr[i] = np.sum(im_arr[:,i])
            int_arr = np.zeros(len(proj_arr))
            for i in range(len(proj_arr)):
                if i == 0:
                    int_arr[0] = proj_arr[0]
                else:
                    int_arr[i] = int_arr[i-1]+proj_arr[i]
                    
            #If this is a BS dataset (see 02506_step14_n1 as an example) then discard
            intfrac_arr = int_arr/int_arr[-1]
            test_ind = (np.where(intfrac_arr<0.1)[0])[-1]
            if test_ind < 25:
                qbeam_arr[s,n] = -101
                ybeam_arr[s,n] = -101
                Ebeam_arr[s,n] = -101
                if doReport:
                    print("Bullshit Detected:")
                    print(" at Step",step,"N",nimg)
                bullshitcount = bullshitcount + 1
            
            else:
                qbeam_arr[s,n]=int_arr[-1]
                int_arr = int_arr/int_arr[-1]
                cent_ind = (np.where(int_arr<0.5)[0])[-1]
                ybeam_arr[s,n]=en_arr[cent_ind]
                if choice > 1:
                    zvals = FACETPlotter.getZarr(superpath, day, dataset)
                    ybeam_cali = FACETPlotter.DTOTR2_ybeamFunc(zvals[s])
                    Ebeam_arr[s,n] = FACETPlotter.DTOTR2_ECali(ybeam_arr[s,n],ybeam_cali)
            
    z_arr = np.linspace(func_start, func_end, n_step)-gasjetpos
    
    if choice == 0:
        ybeam_ave = np.zeros(n_step)
        ybeam_std = np.zeros(n_step)
        for i in range(len(ybeam_ave)):
            ybeam_set = ybeam_arr[i,:]
            ybeam_trim = ybeam_set[np.where(ybeam_set>0)[0]]
            ybeam_ave[i] = np.average(ybeam_trim)
            ybeam_std[i] = np.sqrt(np.var(ybeam_trim))
        
        plt.errorbar(z_arr,ybeam_ave,ybeam_std,label=label_arr[d])
    
    if choice == 1:
        qbeam_ave = np.zeros(n_step)
        qbeam_std = np.zeros(n_step)
        for i in range(len(qbeam_ave)):
            qbeam_set = qbeam_arr[i,:]
            qbeam_trim = qbeam_set[np.where(qbeam_set>0)[0]]
            qbeam_ave[i] = np.average(qbeam_trim)
            qbeam_std[i] = np.sqrt(np.var(qbeam_trim))
        
        plt.errorbar(z_arr,qbeam_ave,qbeam_std,label=label_arr[d])

    if choice == 2:
        Ebeam_ave = np.zeros(n_step)
        Ebeam_std = np.zeros(n_step)
        for i in range(len(Ebeam_ave)):
            Ebeam_set = Ebeam_arr[i,:]
            Ebeam_trim = Ebeam_set[np.where(Ebeam_set>0)[0]]
            Ebeam_ave[i] = np.average(Ebeam_trim)
            Ebeam_std[i] = np.sqrt(np.var(Ebeam_trim))
        
        plt.errorbar(z_arr,Ebeam_ave,Ebeam_std,label=label_arr[d])
        
    if choice == 3:
        ybeam_ave = np.zeros(n_step)
        ybeam_std = np.zeros(n_step)
        for i in range(len(ybeam_ave)):
            ybeam_set = ybeam_arr[i,:]
            ybeam_trim = ybeam_set[np.where(ybeam_set>0)[0]]
            ybeam_ave[i] = np.average(ybeam_trim)
            ybeam_std[i] = np.sqrt(np.var(ybeam_trim))
            
        Ebeam_ave = np.zeros(n_step)
        Ebeam_std = np.zeros(n_step)
        for i in range(len(Ebeam_ave)):
            Ebeam_set = Ebeam_arr[i,:]
            Ebeam_trim = Ebeam_set[np.where(Ebeam_set>0)[0]]
            Ebeam_ave[i] = np.average(Ebeam_trim)
            Ebeam_std[i] = np.sqrt(np.var(Ebeam_trim))
            
        ax1[0].errorbar(z_arr,ybeam_ave,ybeam_std,label=label_arr[d])
        ax1[1].errorbar(z_arr,Ebeam_ave,Ebeam_std,label=label_arr[d])

if choice != 3:
    plt.xlabel("Object Plane Location (m)")
    if choice == 0:
        plt.ylabel("Energy Centroid on Screen (mm)")
    if choice == 1:
        plt.ylabel("Counts on DTOTR2 Camera (arb. units)")
        plt.ylim([550000, 550000*2.3])
    if choice == 2:
        plt.ylabel("Energy Centroid on Screen (GeV)")
    plt.legend()
    #plt.grid(); 
    plt.show()
else:
    ax1[0].set_xlabel("Object Plane Location (m)")
    ax1[1].set_xlabel("Object Plane Location (m)")
    ax1[0].set_ylabel("Energy Centroid on Screen (mm)")
    ax1[1].set_ylabel("Energy Centroid on Screen (GeV)")
    ax1[0].text(-1.12,3.55,"(a)")
    ax1[1].text(-1.12,10.154,"(b)")
    ax1[0].legend(); ax1[1].legend()
    #ax1[0].grid(); ax1[1].grid()
    plt.tight_layout(); plt.show()