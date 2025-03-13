#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 14:54:09 2022

@author: chris
"""

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
import matplotlib.colors as colors


def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'

#dataset = 'E308_02493' #no gas

#dataset = 'E308_02512' #1 psi
#dataset = 'E308_02513' #1 psi, highres

#dataset = 'E308_02492' #6 psi
#dataset = 'E308_02497' #6 psi, highres

#dataset = 'E308_02498' #24 psi
#dataset = 'E308_02500' #24 psi, highres

#dataset = 'E308_02501' #5 bar over 0, 57.8 psi
#dataset = 'E308_02505' #5 bar over 0, 57.8 psi, highres

dataset = 'E308_02506' #9 bar over 0, 115.8 psi
#dataset = 'E308_02511' #9 bar over 0, 115.8 psi, highres, can't do 0.2 or 0.9




#cut_pos = 130
#y_roi = [40,200]
#y_roi = [0,267]

threshold = 20 #Threshold in DTOTR2 in which to ignore camera counts
massage_low = 0.8
massage_high = 40
sig_cent_guess = -25 #for when there are spikes fucking up my algo

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

#Initialize a sigma array that is n_steps long with n_shots of data at each step
beam_slices = np.array([0.1,.20,.30,.40,.50,.60,.70,.80,.90])
#beam_slices = np.array([.30,.40,.50,.60,.70,.80])
sigma_data = np.zeros((len(beam_slices),n_step,n_shot))

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
                except RuntimeError:
                    print("Couldn't Fit at Slice",beam_slices[k],"Step",step,"N",j)
                    plt.plot(slice_axs,slice_arr)
                    plt.show()
                    pass
                if abs(xcoeff[2]) < massage_low:
                    print("Small sig:",xcoeff[2])
                    print(" at Slice",beam_slices[k],"index",slice_ind,"Step",step,"N",nimg)
                    plt.plot(slice_axs,slice_arr,c='b')
                    plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                    plt.show()
                if abs(xcoeff[2]) > massage_high:
                    print("Large sig:",xcoeff[2])
                    print(" at Slice",beam_slices[k],"index",slice_ind,"Step",step,"N",nimg)
                    plt.plot(slice_axs,slice_arr,c='g')
                    plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),c='r')
                    plt.show()
                    plt.set_cmap('gray')
                    plt.imshow(im_arr)
                    #plt.plot([cut_pos,cut_pos],y_roi,c='r')
                    CB = plt.colorbar()
                    #plt.clim(0,1.)
                    plt.show()
                    im_arr2 = im_arr+0.5
                    plt.set_cmap('gist_rainbow')
                    plt.imshow(im_arr2,norm=colors.LogNorm(vmin=im_arr2.min(), vmax=im_arr2.max()))
                    #plt.plot([cut_pos,cut_pos],y_roi,c='r')
                    CB = plt.colorbar()
                    #plt.clim(0,1.)
                    plt.show()
            
                sigma_data[k,i,j] = np.abs(xcoeff[2])