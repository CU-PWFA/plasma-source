#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:32:37 2018

Test loading .tiff images into python

If plt.imread does not work, first install Pillow with 'pip install Pillow'

.tiff images saved as unit16 will have a maximum value of 65535, given by
'np.iinfo(im_arr.dtype).max'.  Dividing the array by this value will instead
give an array of decimals ranging from 0 to 1.

@author: chris
"""

#from PIL import Image as im
import numpy as np
import matplotlib.pyplot as plt

folder = '/home/chris/Desktop/'

imagefile = folder + 'June172019Data/4mbar June 17 2019/SID4 11h13m09s018ms.tif'
#imagefile = folder + 'CameraPics/Cam01_1805170002_0001.tiff'
image = plt.imread(imagefile)

im_arr = np.array(image)#[140:200,540:600]
#max_tiffunit = np.iinfo(im_arr.dtype).max
#im_arr = im_arr/max_tiffunit

im_arr = im_arr/16

plt.set_cmap('gray')
plt.imshow(im_arr)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.show()