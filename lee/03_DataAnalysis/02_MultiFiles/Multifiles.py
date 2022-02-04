#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:32:19 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from photutils.centroids import centroid_com
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from scipy import optimize
import os
import glob
from scipy.ndimage.interpolation import shift

#%%
PixelSize= 3e-6
#%%
DataDir= '/media/valentinalee/Elements/PlasmaDiagnosticsData/day_29/2110290047/' 
BGDir= '/media/valentinalee/Elements/PlasmaDiagnosticsData/day_29/2110290048/'

shotN= 5
#%%
os.chdir(os.path.expanduser(BGDir))
for name in glob.glob('17529184_??????????_????.tiff'):
    DataNumber= name[9:19]    
print(DataNumber)

for shot in range (0, shotN):
    strN= str(shot)
    FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
    print(FileName)
    if shot==0:
        BG= np.array(Image.open(FileName))
        TotalBG= np.zeros((BG.shape[0], BG.shape[1]))
        BG_CX0, BG_CY0= centroid_com(BG)
    
    BG= np.array(Image.open(FileName))    
    BG_cx, BG_cy= centroid_com(BG)
#move BG so that centroid goes to center
    BG= shift(BG, ((BG_CY0-BG_cy), (BG_CX0-BG_cx)))
    TotalBG= TotalBG+ BG

AvgBG= TotalBG/shotN
#%%
def shiftBG(s1, s2):
    ShiftedBG= shift(AvgBG, (s1, s2))
    return ShiftedBG

os.chdir(os.path.expanduser(DataDir))
for name in glob.glob('17529184_??????????_????.tiff'):
    DataNumber= name[9:19]    
print(DataNumber)

for shot in [10]:
#for shot in range (0, shotN):
    if shot== 0:
        strN= str(shot)
        FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
        print(FileName)
        
    #load the tiff image and change it to array
        data= np.array(Image.open(FileName))
        errorFn= lambda p1: np.ravel(shiftBG(*p1)- data)
        p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
        print(p)
        data= shift(data, (-p[0], -p[1]))
        signal= data-AvgBG
        signalFiltered= gaussian_filter(signal, sigma=6)
#TODO figure out how to set range and peak finding 
        sumSignal= np.sum(-signalFiltered, axis= 1)[200:400]
        peaks, _ = find_peaks(sumSignal, prominence= (np.amax(sumSignal)-np.amin(sumSignal))*0.01, width= (5, 40))
        if len(peaks)>1:
            PeakInt= [np.sum(-signalFiltered, axis= 1)[int(200+peaks[n])] for n in range (0, len(peaks))]
            sortedArray= sorted(range(len(PeakInt)), key=lambda k: PeakInt[k])
            Peak= int(peaks[sortedArray[-1]]+200)
        elif len(peaks)==1:
            Peak= peaks[0]+200
        else:
            print('No peak found in '+FileName)
            continue
        Peak0= Peak
        TotalData= np.zeros((data.shape[0], data.shape[1]))
        
    strN= str(shot)
    FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
    print(FileName)
    
#load the tiff image and change it to array
    data= np.array(Image.open(FileName))
    errorFn= lambda p1: np.ravel(shiftBG(*p1)- data)
    p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
    print(p)
    data= shift(data, (-p[0], -p[1]))
    signal= data-AvgBG
    signalFiltered= gaussian_filter(signal, sigma=6)
    sumSignal= np.sum(-signalFiltered, axis= 1)[200:400]
    peaks, _ = find_peaks(sumSignal, prominence= (np.amax(sumSignal)-np.amin(sumSignal))*0.01, width= (5, 40))
    if len(peaks)>1:
        PeakInt= [np.sum(-signalFiltered, axis= 1)[int(200+peaks[n])] for n in range (0, len(peaks))]
        sortedArray= sorted(range(len(PeakInt)), key=lambda k: PeakInt[k])
        Peak= int(peaks[sortedArray[-1]]+200)
    elif len(peaks)==1:
        Peak= peaks[0]+200
    else:
        print('No peak found in '+FileName)
        continue
    print(Peak)
    movedData= shift(data, (Peak0-Peak, 0))
    TotalData= TotalData+movedData

#%%
AvgData= TotalData/shotN
errorFn= lambda p1: np.ravel(shiftBG(*p1)- AvgData)
p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
print(p)
AvgData= shift(AvgData, (-p[0], -p[1]))
signal= AvgData-AvgBG
signalFiltered= gaussian_filter(signal, sigma=2)
np.save(os.path.join(os.path.expanduser(ResultDir_npy), DataNumber+'_Avg'), \
        (signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)))
#%%
BaseResultDir= '~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/'
ResultDir_npy= BaseResultDir+'02_Plots/'
plt.pcolormesh((signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)), \
               cmap= 'bwr')
plt.axis('scaled')
plt.colorbar()
plt.title(DataNumber+'_Avg')
plt.savefig(os.path.expanduser(ResultDir_npy)+DataNumber)
print('3')
plt.close()
#%%
ResultDir= '~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/'
np.save(os.path.join(os.path.expanduser(ResultDir), DataNumber+'_Avg'), (signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)))

#%%
np.save('ExampleSignal', (signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)))
#%%
plt.figure()
plt.pcolormesh(signalFiltered, cmap= 'bwr')
plt.axis('scaled')
plt.colorbar()
#plt.savefig('Example.png', bbox_inches='tight')
#%%
Start= 100
Stop= 341
FirstPeaks= np.zeros((Stop-Start, 3))
#for x in [225]:
for x in range(Start, Stop):
    peaks, _= find_peaks(signalFiltered[:, x])
    PeakInt= [signalFiltered[:, x][peaks[n]] for n in range (0, len(peaks))]
    sortedArray= sorted(range(len(PeakInt)), key=lambda k: PeakInt[k])
    peak1= int(peaks[sortedArray[-1]])
    peak2= int(peaks[sortedArray[-2]])
    FirstPeaks[int(x-Start), 0]= x
    FirstPeaks[int(x-Start), 1]= peak1
    FirstPeaks[int(x-Start), 2]= peak2

#%%
plt.figure()
#plt.pcolormesh(AvgData)
plt.pcolormesh(signalFiltered, cmap= 'bwr')
plt.plot(FirstPeaks[:, 0], FirstPeaks[:, 1])
plt.plot(FirstPeaks[:, 0], FirstPeaks[:, 2])
plt.colorbar()

#%%
plt.figure()
plt.plot(np.sum(-signalFiltered, axis= 1))
