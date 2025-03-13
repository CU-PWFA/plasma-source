#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 10:59:28 2022

Loading in images of TopView

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
#dataset = 'E308_02484' #1 psi
#dataset = 'E308_02500' #24 psi
#dataset = 'E308_02492' #6psi, lowres
dataset = 'E308_02506' #9 bar over 0, 115.8 psi

back = False
annotation = "(c) 115.8 psi, Electron Beam On, Laser On"

dx = 54.3
x0 = 50

step = 1
for i in range(1):
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
    if back:
        print("Loading Background:")
        path = superpath + day + dataset + '/' + dataset +'.mat'
        mat = sio.loadmat(path)
        
        data = mat['data_struct']
        backs = data['backgrounds'][0][0]
        dtotr2_back = backs['TopView'][0][0]
        dtotr2_back = np.transpose(dtotr2_back)
        
        #plt.set_cmap('gray')
        #plt.imshow(dtotr2_back,vmin = 0,vmax = 2000)
        #CB = plt.colorbar(orientation='horizontal')
        #plt.show()
        im_arr = np.array(dtotr2_back)
    #######################
    else:
        camerapath = 'TopView'
        print("Loading Step",step,"at n",nimg,":")
        if step < 10:
            image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
        else:
            image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
        
        path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
        
        image = plt.imread(path)
        
        im_arr = np.array(image)
    
    """
    threshold = 20
    for i in range(im_arr.shape[0]):
        for j in range(im_arr.shape[1]):
            im_arr[i][j] = im_arr[i][j] - dtotr2_back[i][j]
            if im_arr[i][j] > 65000:
                im_arr[i][j] = 0
            if im_arr[i][j] < threshold:
                im_arr[i][j] = 0
    """
    
    #im_arr = np.flip(im_arr,1)
    
    hori_range = (np.array([0,im_arr.shape[1]])-x0)/dx
    vert_range = (np.array([0,im_arr.shape[0]])-1/2*im_arr.shape[0])/dx
    
    plt.set_cmap('magma')
    #plt.set_cmap('jet')
    plt.imshow(im_arr,vmin=0,vmax=2100,extent=[hori_range[0],hori_range[1],vert_range[0],vert_range[1]])
    #CB = plt.colorbar(orientation='horizontal')
    #plt.clim(0,1.)
    plt.xlabel("Z (mm)")
    plt.ylabel("X (mm)")
    #plt.arrow(0,0,2,0,color='white',arrowstyle='->')
    plt.annotate("",xy=(5,0),xytext=(-0.7,0),arrowprops=dict(color="white",arrowstyle="->"))
    plt.text(-0.2,1.9,annotation,color='white')
    plt.show()
sys.exit()
im_arr2 = im_arr+0.5
plt.set_cmap('gist_rainbow')
plt.imshow(im_arr2,norm=colors.LogNorm(vmin=im_arr2.min(), vmax=im_arr2.max()))
CB = plt.colorbar(orientation='horizontal')
#plt.clim(0,1.)
plt.show()