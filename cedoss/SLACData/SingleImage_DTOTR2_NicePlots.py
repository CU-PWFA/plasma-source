#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 13:08:24 2022

Loading images of DTOTR2, trying to make a nice plot

OUTDATED: USE FACETPlotter INSTEAD!!

This was just to quickly make some movies back when I was first making these scripts,
 and I didn't correctly make the energy lines yet.  They are still square and bad

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

sys.exit()

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
#dataset = 'E308_02493' #no gas
dataset = 'E308_02512' #1 psi
#dataset = 'E308_02492' #6psi, lowres
#dataset = 'E308_02506' #9 bar over 0, 115.8 psi
#dataset = 'E308_02497' #6 psi, highres

dx = 30.3e-6*1e3
ybeam = 3.1
dnom = 60
ebend = 10

def ECali(yscreen):
    dy = dnom-ybeam
    return dnom*ebend/(yscreen + dy)

#5 and 14 for at high divergence point
#3 and 8 for minimum at 90%

step = 1
for i in range(10):
    nimg = i
    
    path = superpath + day + dataset + '/' + dataset +'.h5'
    """
    f = h5.File(path,"r")
    #datasetNames = [n for n in f.keys()]
    #for n in datasetNames:
    #    print(n)
        
    data = f['data']
    #datasetNames = [n for n in data.keys()]
    #for n in datasetNames:
    #    print(n)
        
    save_info = data['save_info']
    datasetNames = [n for n in save_info.keys()]
    for n in datasetNames:
        print(n)
    print()
    params = data['params']
    datasetNames = [n for n in params.keys()]
    for n in datasetNames:
        print(n)
        
    print()
    print("Comment:")
    print(params['comment'][0][0])
    
    print()
    print("Scan Function: ",params['scanFuncs'][0][0])
    print(" ",params['startVals'][0][0]," to ",params['stopVals'][0][0])
    
    camerapath = 'TopView'
    path = superpath + day + dataset + '/images/' + camerapath + '/'
    
    image_list = os.listdir(path)
    n_shot = int(params['n_shot'][0][0])
    print("Shots per Step: ",n_shot)
    n_step = int(len(image_list)/n_shot)
    print("Num Steps: ",n_step)
    
    f.close()
    print()
    """
    ###BACKGROUND FROM MAT FILE
    
    print("Loading Background:")
    path = superpath + day + dataset + '/' + dataset +'.mat'
    mat = sio.loadmat(path)
    
    data = mat['data_struct']
    backs = data['backgrounds'][0][0]
    dtotr2_back = backs['DTOTR2'][0][0]
    dtotr2_back = np.transpose(dtotr2_back)
    """
    plt.set_cmap('gray')
    plt.imshow(dtotr2_back,vmin = 0,vmax = 2000)
    CB = plt.colorbar(orientation='horizontal')
    plt.show()
    """
    #######################
    camerapath = 'DTOTR2'
    print("Loading Step",step,"at n",nimg,":")
    if step < 10:
        image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
    else:
        image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
    
    path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
    
    image = plt.imread(path)
    
    im_arr = np.array(image)
    
    
    threshold = 20
    for i in range(im_arr.shape[0]):
        for j in range(im_arr.shape[1]):
            im_arr[i][j] = im_arr[i][j] - dtotr2_back[i][j]
            if im_arr[i][j] > 65000:
                im_arr[i][j] = 0
            if im_arr[i][j] < threshold:
                im_arr[i][j] = 0
    
    
    #im_arr = np.transpose(im_arr)
    #im_arr = im_arr[int(1/4*im_arr.shape[0]):int(3/4*im_arr.shape[0]),:]
    
    if ybeam is None:
        hori_range = (np.array([0,im_arr.shape[1]]))*dx
        print("set ybeam!")
    else:
        hori_range = ECali((np.array([0,im_arr.shape[1]]))*dx)
    vert_range = (np.array([0,im_arr.shape[0]])-1/2*im_arr.shape[0])*dx
    
    aspect_ratio = 0.5*1/np.abs((vert_range[1]-vert_range[0])/(hori_range[1]-hori_range[0]))
    
    plt.set_cmap('OrRd')
    plt.imshow(im_arr,aspect=aspect_ratio,vmin=0,vmax=1000,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    #CB = plt.colorbar(orientation='horizontal')
    #plt.clim(0,1.)
    plt.xlabel("Energy (GeV)")
    plt.ylabel("Position on Screen (mm)")
    #plt.show()
    """
    im_arr2 = im_arr+0.5
    plt.set_cmap('gist_rainbow')
    plt.imshow(im_arr2,aspect=aspect_ratio,norm=colors.LogNorm(vmin=im_arr2.min(), vmax=im_arr2.max()),vmin=0.5,vmax=1500,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    #CB = plt.colorbar(orientation='horizontal')
    #plt.clim(0,1.)
    #plt.show()
    """
    
    en_arr = np.linspace(hori_range[0],hori_range[1],im_arr.shape[1])
    proj_arr = np.zeros(im_arr.shape[1])
    for i in range(im_arr.shape[1]):
        proj_arr[i] = (np.sum(im_arr[:,i])*(1/7000))+vert_range[0]
        
    plt.plot(en_arr,proj_arr)
    #plt.xlabel("Energy Slice (pixel)")
    #plt.ylabel("Total Camera Counts")
    plt.grid();plt.show()