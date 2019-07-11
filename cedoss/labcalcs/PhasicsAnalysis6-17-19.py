#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:42:16 2019

Analysis of diffraction patterns in Phasics data taken on June 17 2019.

@author: chris
"""

#from PIL import Image as im
import numpy as np
import matplotlib.pyplot as plt

folder = '/home/chris/Desktop/'

#imagefile = folder + 'June172019Data/20mbar June 17 2019/SID4 09h38m47s039ms.tif'  #20mbar
imagefile = folder + 'June172019Data/12mbar June 17 2019/SID4 11h26m27s058ms.tif'  #12mbar
#imagefile = folder + 'June172019Data/4mbar June 17 2019/SID4 11h13m38s018ms.tif'  #4mbar
image = plt.imread(imagefile)

im_arr = np.array(image)#[140:200,540:600]
#max_tiffunit = np.iinfo(im_arr.dtype).max
#im_arr = im_arr/max_tiffunit

im_arr = im_arr

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(im_arr)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.show()

plt.figure(figsize=(20,10))
plt.plot(im_arr[:,600])
plt.plot(im_arr[:,595])
plt.plot(im_arr[:,590])
plt.show()

nosamples = 30
samplefreq = 1
centerline = 600

dataline = im_arr[:,centerline]
samples = np.zeros(2*nosamples + 1)
for i in range(len(dataline)):
    samples[0] = dataline[i]
    for jj in range(nosamples):
        j = jj + 1
        samples[2*j-1] = im_arr[i,centerline - samplefreq*j]
        samples[2*j]=im_arr[i,centerline + samplefreq*j]
    dataline[i] = np.average(samples)
    
plt.figure(figsize=(20,10))
plt.plot(dataline)
plt.show()