#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:40:59 2021

@author: valentinalee
"""
#%%
import numpy as np
from scipy.signal import find_peaks

def FindShadowgraphyMrm (MidLineOutI, PixelSize):

    x= MidLineOutI
    peaks, _ = find_peaks(x)
    peak1idx, peak2idx= peaks[np.argsort(x[peaks])[-1]], peaks[np.argsort(x[peaks])[-2]]
    peak1, peak2= x[peak1idx], x[peak2idx]
    
    if (np.argsort(x[peaks])[-1])-(np.argsort(x[peaks])[-2])==2:
       
        valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                +min(peak1idx, peak2idx)
        valley1idx, valley2idx= valleys[np.argsort(x[valleys])[0:2]]
        valley1, valley2= x[np.array([valley1idx, valley2idx])]
    else:
        valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                +min(peak1idx, peak2idx)
        valley1idx= valleys[np.argsort(x[valleys])[0]]
        valley2idx= valley1idx
        valley1, valley2= x[np.array([valley1idx, valley2idx])]
        
    if (np.mean(np.array([peak1, peak2])))-(np.mean(np.array([valley1, valley2])))>0.5:
        
        FirstPeakWidth= abs(peak1idx-peak2idx)*PixelSize
        PTValue= 100*(np.mean(np.array([peak1, peak2]))-np.mean(np.array([valley1, valley2])))/ \
                    ((np.mean(np.array([peak1, peak2]))+np.mean(np.array([valley1, valley2])))/2)
        return FirstPeakWidth, PTValue