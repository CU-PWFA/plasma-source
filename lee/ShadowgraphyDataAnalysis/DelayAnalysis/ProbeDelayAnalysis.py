#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 13:53:41 2021

@author: valentinalee
"""

import numpy as np;
import matplotlib.pyplot as plt
from photutils.centroids import centroid_com
from PIL import Image
from scipy.ndimage.interpolation import shift
from scipy.ndimage import gaussian_filter
import cv2
from skimage.transform import hough_line, hough_line_peaks
from skimage.feature import canny
#%%
#load background data, find centroid 
BG= np.array(Image.open('2109170088/19423598_2109170088_0000.tiff'))
BG_cx, BG_cy = centroid_com(BG)

#move BG so that centroid goes to center
BG= shift(BG, ((300-BG_cy), (300-BG_cx)))

#load data, find centroid
Data= np.array(Image.open('2109170089/19423598_2109170089_0000.tiff'))
Data_cx, Data_cy = centroid_com(Data)
 #if it did work, use a low pass filter first 
#move data so that centroid goes to center
Data= shift(Data, ((300-Data_cy), (300-Data_cx)))

#sutraction 
Signal= gaussian_filter(Data-BG, sigma=2)

#%%
plt.figure(2)
plt.pcolormesh(Signal)

#%%
#image = Signal[300:365, 105:550]
image = Signal[250:320, 105:550]
tested_angles = np.linspace(np.deg2rad(85), np.deg2rad(90), 50, endpoint=False)
h, theta, d = hough_line(image, theta=tested_angles)

for _, angle, dist in zip(*hough_line_peaks(h, theta, d, threshold=450)):
    (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
    print((x0, y0), np.tan(angle + np.pi/2))
#%%
fig, axes = plt.subplots(1, 3, figsize=(15, 6))
ax = axes.ravel()

ax[2].imshow(image, cmap=cm.gray)
ax[2].set_ylim((image.shape[0], 0))
ax[2].set_axis_off()
ax[2].set_title('Detected lines')

for _, angle, dist in zip(*hough_line_peaks(h, theta, d, threshold=450)):
    (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
    ax[2].axline((x0, y0), slope=np.tan(angle + np.pi/2))

