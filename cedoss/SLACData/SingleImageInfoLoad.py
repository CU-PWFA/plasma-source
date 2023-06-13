#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:30:22 2022

Testing going through the SLAC data

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
    return A*np.exp(-(x-mu)**2/(2.*(sigma)**2))

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'
dataset = 'E308_02493' #no gas
#dataset = 'E308_02500' #24 psi
#dataset = 'E308_02492' #6psi, lowres
#dataset = 'E308_02506' #9 bar over 0, 115.8 psi


step = 3
nimg = 2
cut_pos = 27
y_roi = [40,200]
slice_pos = 0.3

path = superpath + day + dataset + '/' + dataset +'.h5'

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

camerapath = 'DTOTR2'
path = superpath + day + dataset + '/images/' + camerapath + '/'

image_list = os.listdir(path)
n_shot = int(params['n_shot'][0][0])
print("Shots per Step: ",n_shot)
n_step = int(len(image_list)/n_shot)
print("Num Steps: ",n_step)

#Don't need to, just get the total number of files in the folder
"""
step_list = []
for i in range(len(image_list)):
    name_split = image_list[i].split('_')
    stepno = name_split[2].split('p')
    step_list.append(int(stepno[1]))

step_list = np.sort(step_list)
print(step_list)
"""

f.close()
print()

###BACKGROUND FROM MAT FILE

print("Loading Background:")
path = superpath + day + dataset + '/' + dataset +'.mat'
mat = sio.loadmat(path)

data = mat['data_struct']
backs = data['backgrounds'][0][0]
dtotr2_back = backs[camerapath][0][0]
dtotr2_back = np.transpose(dtotr2_back)

plt.imshow(dtotr2_back)
CB = plt.colorbar()
plt.show()

#######################


print("Loading Step",step,"at n",nimg,":")
if step < 10:
    image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
else:
    image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"

path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path

image = plt.imread(path)

im_arr = np.array(image)#[140:200,540:600]
#max_tiffunit = np.iinfo(im_arr.dtype).max
#im_arr = im_arr/max_tiffunit
"""
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
"""
#im_arr = im_arr/16
threshold = 20
for i in range(im_arr.shape[0]):
    for j in range(im_arr.shape[1]):
        if im_arr[i][j]*1. - dtotr2_back[i][j]*1. > 0:
            im_arr[i][j] = im_arr[i][j] - dtotr2_back[i][j]
        else:
            im_arr[i][j] = 0
        #if im_arr[i][j] > 65000:
        #    im_arr[i][j] = 0
        if im_arr[i][j] < threshold:
            im_arr[i][j] = 0

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

slice_arr = im_arr[y_roi[0]:y_roi[1],cut_pos]
slice_axs = np.arange(len(slice_arr))-int(.5*len(slice_arr))

p0 = [1., 0., 1.]
xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)

plt.plot(slice_axs,slice_arr,label="Data at Slice = "+str(cut_pos))
plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),label="Gaussian Fit, "+r'$\sigma=$'+str(round(xcoeff[2],2)))
plt.xlabel("Beam Image in X (pixels)")
plt.ylabel("Intensity (counts)")
plt.legend();plt.show()



print("Entire Projection Along x")

x_arr = np.zeros(im_arr.shape[0])
xproj_arr = np.zeros(im_arr.shape[0])
for i in range(im_arr.shape[0]):
    x_arr[i] = i
    xproj_arr[i] = np.sum(im_arr[i,:])

p0 = [xproj_arr[np.argmax(xproj_arr)], x_arr[np.argmax(xproj_arr)], 1.]
xcoeff, var_matrix = curve_fit(Gauss, x_arr, xproj_arr, p0=p0)
plt.plot(x_arr,xproj_arr,label="Projection in X")
plt.plot(x_arr,Gauss(x_arr,*xcoeff),label="Gaussian Fit, "+r'$\sigma=$'+str(round(xcoeff[2],2)))
plt.legend();plt.show()

m_per_pixel = 30.3e-6
magnification = 2.43

print(" Sigma in pixels:    ",xcoeff[2])
print(" Sigma in um:        ",xcoeff[2]*m_per_pixel*1e6)
print(" Actual sigma in um: ",xcoeff[2]*m_per_pixel*1e6/magnification)
print()

print("Energy Projection:")

en_arr = np.zeros(im_arr.shape[1])
proj_arr = np.zeros(im_arr.shape[1])
for i in range(im_arr.shape[1]):
    en_arr[i] = i
    proj_arr[i] = np.sum(im_arr[:,i])
    
plt.plot(en_arr,proj_arr)
plt.xlabel("Energy Slice (pixel)")
plt.ylabel("Total Camera Counts")
plt.show()

int_arr = np.zeros(len(proj_arr))
for i in range(len(proj_arr)-1):
    int_arr[i+1] = proj_arr[i+1] + int_arr[i]
total_counts = int_arr[-1]
int_arr = int_arr/int_arr[-1]

plt.plot(en_arr,int_arr)
plt.xlabel("Energy Slice (pixel)")
plt.ylabel("Percentage of the Beam (%)")
plt.grid();plt.show()

cut_pos = (np.where(int_arr<slice_pos)[0])[-1]

slice_arr = im_arr[:,cut_pos]
slice_axs = np.arange(len(slice_arr))-int(.5*len(slice_arr))

p0 = [1., 0., 1.]
xcoeff, var_matrix = curve_fit(Gauss, slice_axs, slice_arr, p0=p0)

plt.plot(slice_axs,slice_arr,label="Data at Slice = "+str(cut_pos))
plt.plot(slice_axs,Gauss(slice_axs,*xcoeff),label="Gaussian Fit, "+r'$\sigma=$'+str(round(xcoeff[2],2)))
plt.xlabel("Beam Image in X (pixels)")
plt.ylabel("Intensity (counts)")
plt.legend();plt.show()

print("Total Camera Counts:",total_counts)

