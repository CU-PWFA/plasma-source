#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 11:06:11 2021

@author: valentinalee
"""

import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def FindPeaksFromSeed(coords, Signal):
    '''
    Coords: Dictionay, {'x0':_,'y0':_,'x1':_,'y1':_}
    Signal: 2D array
    '''
    LineCoord= [coords['x0'], coords['y0'], coords['x1'], coords['y1']]
    
    
    if LineCoord[0]> LineCoord[2]:
        x1= LineCoord[0]
        y1= LineCoord[1]
        x0= LineCoord[2]
        y0= LineCoord[3]
        Start= LineCoord[2]
        Stop= LineCoord[0]
    else:
        Start= LineCoord[0]
        Stop= LineCoord[2]
        x0= LineCoord[0]
        y0= LineCoord[1]
        x1= LineCoord[2]
        y1= LineCoord[3]
    
    FirstPeaks= np.zeros((abs(Stop-Start)+1, 2))
    a= (y0-y1)/(x0-x1)
    b= y0-a*x0
    
    for x in range(Start, Stop+1):
        seed= a*x+b
        peaks, _= find_peaks(Signal[:, x])
        DisToSeed= [abs(p-seed) for p in peaks]
        sortedArray= sorted(range(len(DisToSeed)), key=lambda k: DisToSeed[k])
        peak= int(peaks[sortedArray[0]])
        FirstPeaks[int(x-Start), 0]= x
        FirstPeaks[int(x-Start), 1]= peak
    return FirstPeaks

#%%
#Signal= np.load('ExampleSignal.npy')
#Coords= {'x0': 108, 'y0': 283, 'x1': 371, 'y1': 280}
#FirstPeaks= FindPeaksFromSeed(Coords, Signal)
#plt.figure()
#plt.pcolormesh(Signal, cmap= 'bwr')
#plt.plot(FirstPeaks[:, 0], FirstPeaks[:, 1])
#%%
'''
When click confirmLine1, it calls the function to find peaks around that line.
plot the peaks found 
save those points in memery 
and when click confirmLine2, it does the same

then click "confirm peaks found" to save the array that contains the peaks coords for both lines
'''

'''
load signal and load peaks coords
calculate the slopes
calculate the opening angle
rotate the signal by whatever 
decide where to measure width
decide signal range that is going to be used for calculating the signal strength
calculate the avg signal strength
'''