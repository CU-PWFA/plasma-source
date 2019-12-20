#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:30:51 2019

@author: chris
"""

#from PIL import Image as im
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from scipy.optimize import curve_fit

folder = '/home/chris/Desktop/cyllens_test1/'

def getAverageSet(directory, case):
    setlist = [f for f in listdir(directory+case) if isfile(join(directory+case, f))]
    image = None
    for set_i in setlist:
        imagefile = directory + case + set_i
        if image == None:
            image = np.array(plt.imread(imagefile),dtype='double')
        else:
            image = np.array(plt.imread(imagefile),dtype='double') + image
    return image / len(setlist)

image = getAverageSet(folder, '1911060004/')

plt.set_cmap('gray')
plt.imshow(image)
plt.title("Averaged Image")
CB = plt.colorbar()
plt.show()

image=image-image[0,0]

maxtot = np.argmax(image)

ysize, xsize = np.shape(image)

row = int(maxtot/xsize)
col = maxtot-xsize*row

window = 70
pixsize = 3.29

widdat = image[row-window:row+window,col]
nardat = image[row,col-window:col+window]

axis = np.arange(-window,window,1)*3.29

plt.plot(axis,nardat)
plt.plot(axis,widdat)
plt.xlabel("ums")
plt.ylabel("Intensity")
plt.title("Intensity Lineouts from camera")
plt.show()

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


p0 = [1., 0., 1.]
        
narcoeff, var_matrix = curve_fit(gauss, axis, nardat, p0=p0) 
narfit = gauss(axis, *narcoeff)

widcoeff, var_matrix = curve_fit(gauss, axis, widdat, p0=p0) 
widfit = gauss(axis, *widcoeff)

print("Narrow Standard Deviation: ",narcoeff[2],"um")
print("Wide Standard Deviation:   ",widcoeff[2],"um")

plt.plot(axis,nardat)
plt.plot(axis,gauss(axis,*narcoeff))
plt.show()

plt.plot(axis,widdat)
plt.plot(axis,gauss(axis,*widcoeff))
plt.show()
