#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 16:40:35 2021

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from photutils.centroids import centroid_com
from PIL import Image
from scipy.ndimage.interpolation import shift
from scipy.ndimage import gaussian_filter
#%%
#load background data, find centroid 
BG= np.array(Image.open('2109170116/19423598_2109170116_0000.tiff'))
BG_cx, BG_cy = centroid_com(BG)

#move BG so that centroid goes to center
BG= shift(BG, ((300-BG_cy), (300-BG_cx)))

#load data, find centroid
Data= np.array(Image.open('2109170117/19423598_2109170117_0000.tiff'))
Data_cx, Data_cy = centroid_com(Data)
 #if it did work, use a low pass filter first 
#move data so that centroid goes to center
Data= shift(Data, ((300-Data_cy), (300-Data_cx)))

#sutraction 
Signal= gaussian_filter(Data-BG, sigma=5)

#%%
BGList= []
#for loop load 100 shots background
for shot_n in range (0, 10):
    number= str(shot_n).zfill(4)  
    BG= np.array(Image.open('2109170116/19423598_2109170116_'+number+'.tiff'))
    BG_cx, BG_cy = centroid_com(BG)
#each shot, move and shfit it    
#move BG so that centroid goes to center
    BGList.append(shift(BG, ((300-BG_cy), (300-BG_cx))))
#calculate the average
AvgBG= np.mean(BGList, axis=0)
#for loop load 100 shots data

SignalRangeX= [180, 430]
SignalRangeY= [215, 405]
CentralIdx= (SignalRangeY[1]-SignalRangeY[0])/2

for shot_n in range (0, 10):
    number= str(shot_n).zfill(4)  
    Data= np.array(Image.open('2109170117/19423598_2109170117_'+number+'.tiff'))
#each shot, move it, subtract the background
    Data_cx, Data_cy = centroid_com(Data)
    Data= shift(Data, ((300-Data_cy), (300-Data_cx)))
    Signal= gaussian_filter(Data-AvgBG, sigma=5)
#each signal, find peaks
    PeakIdx1= []
    PeakIdx2= []
    ValleyIdx= []
    Peak1= []
    Peak2= []
    Valley= []
    for XIdx in range (SignalRangeX[0], SignalRangeX[1]):
        print(XIdx)
        x= Signal[SignalRangeY[0]: SignalRangeY[1], XIdx]
        peaks, _ = find_peaks(x, height=4000, distance=20)
    
        x2= -x
        valleys, _ = find_peaks(x2)
    
        try:
            valleyidx= valleys[np.argsort(abs(valleys-CentralIdx))[0]]
            valley= x[valleyidx]
        except:
            pass
        
        try:
            if peaks[np.argsort(abs(peaks-valleyidx))[0]] > peaks[np.argsort(abs(peaks-valleyidx))[1]]:
                peak1idx, peak2idx= peaks[np.argsort(abs(peaks-valleyidx))[0]], peaks[np.argsort(abs(peaks-valleyidx))[1]]
            else:
                peak2idx, peak1idx= peaks[np.argsort(abs(peaks-valleyidx))[0]], peaks[np.argsort(abs(peaks-valleyidx))[1]]
            peak1, peak2= x[peak1idx], x[peak2idx]
        except:
            pass
    
        PeakIdx1.append([peak1idx, XIdx])
        PeakIdx2.append([peak2idx, XIdx])
        ValleyIdx.append([valleyidx, XIdx])
        Peak1.append([peak1, XIdx])
        Peak2.append([peak2, XIdx])
        Valley.append([valley, XIdx])
    PeakIdx1= np.array(PeakIdx1)
    PeakIdx2= np.array(PeakIdx2)
    ValleyIdx= np.array(ValleyIdx)
    Peak1= np.array(Peak1)
    Peak2= np.array(Peak2)
    Valley= np.array(Valley)

#Do a fitting for the first line
#move the signal 
#%%
PeakWidth= np.mean(PeakIdx1[:, 0]-PeakIdx2[:, 0])
PeakToTrough= np.mean(((Peak1[:, 0]+Peak2[:, 0])/2-Valley[:, 0])/Valley[:, 0])
print('PeakWidth', PeakWidth)
print('PeakToTrough', PeakToTrough)
#move the peaks
#overlap each signal calculate the average 