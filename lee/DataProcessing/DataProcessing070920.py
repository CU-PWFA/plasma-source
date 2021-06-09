#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 18:40:09 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy import optimize
from sklearn.linear_model import LinearRegression
import statistics 
#%%10
test= np.array(Image.open('2006260010/19423601_2006260010_0000.tiff'))
BG= np.array(Image.open('2006260011/19423601_2006260011_0001.tiff'))
#%%12
test= np.array(Image.open('2006260012/19423601_2006260012_0002.tiff'))
BG= np.array(Image.open('2006260013/19423601_2006260013_0000.tiff'))
#%%14
test= np.array(Image.open('2006260014/19423601_2006260014_0001.tiff'))
BG= np.array(Image.open('2006260015/19423601_2006260015_0000.tiff'))
#%%20
test= np.array(Image.open('2006260020/19423601_2006260020_0001.tiff'))
BG= np.array(Image.open('2006260021/19423601_2006260021_0000.tiff'))
#%%24
test= np.array(Image.open('2006260024/19423601_2006260024_0006.tiff'))
BG= np.array(Image.open('2006260025/19423601_2006260025_0001.tiff'))
#%%26
test= np.array(Image.open('2006260026/19423601_2006260026_0001.tiff'))
BG= np.array(Image.open('2006260027/19423601_2006260027_0001.tiff'))
#%%28
test= np.array(Image.open('2006260028/19423601_2006260028_0007.tiff'))
BG= np.array(Image.open('2006260029/19423601_2006260029_0001.tiff'))
#%%32
test= np.array(Image.open('2006260032/19423601_2006260032_0001.tiff'))
BG= np.array(Image.open('2006260033/19423601_2006260033_0002.tiff'))
#%%38
test= np.array(Image.open('2006260038/19423601_2006260038_0001.tiff'))
BG= np.array(Image.open('2006260039/19423601_2006260039_0001.tiff'))
#%%42
test= np.array(Image.open('2006260042/19423601_2006260042_0001.tiff'))
BG= np.array(Image.open('2006260043/19423601_2006260043_0001.tiff'))
#%%44
test= np.array(Image.open('2006260044/19423601_2006260044_0002.tiff'))
BG= np.array(Image.open('2006260045/19423601_2006260045_0000.tiff'))
#%%48
test= np.array(Image.open('2006260048/19423601_2006260048_0001.tiff'))
BG= np.array(Image.open('2006260049/19423601_2006260049_0000.tiff'))
#%%58
test= np.array(Image.open('2006260058/19423601_2006260058_0001.tiff'))
BG= np.array(Image.open('2006260059/19423601_2006260059_0003.tiff'))
#%%
BGfreq= fftshift(fft2(BG))
testfreq= fftshift(fft2(test))
Newfreq= testfreq-BGfreq
NewS= abs(ifft2(fftshift(Newfreq)))
SignalArea= NewS[500:800, :]
Newsignal= gaussian_filter(SignalArea, sigma=5)
NewsignalBG= gaussian_filter(SignalArea, sigma=20)
NS= Newsignal-NewsignalBG
#%%
plt.figure(3)
plt.pcolormesh(NS)
plt.axis('scaled')
plt.colorbar()
#%%
IntCenter= 195
CIdx= IntCenter
count=0
IdxStart= 600
IdxEnd= 1400
DataLength= int((IdxEnd-IdxStart)/10)
PeakIdxArray_R= np.zeros((DataLength, 2))
PeakIdxArray_L= np.zeros((DataLength, 2))

for Xidx in range (IdxStart, IdxEnd, 10):
    LONew=NS[:, Xidx]
    for n in range (CIdx, int(CIdx+50)):
        if LONew[n+1]-LONew[n] <0:
            if LONew[n+2]-LONew[n] <0:
                if LONew[n+3]-LONew[n] <0:            
                    if LONew[n+4]-LONew[n] <0:            
                        if LONew[n+5]-LONew[n] <0:            
                            if LONew[n-1]-LONew[n] <0:            
                                if LONew[n-2]-LONew[n] <0:            
                                    if LONew[n-3]-LONew[n] <0:            
                                        if LONew[n-4]-LONew[n] <0:            
                                            if LONew[n-5]-LONew[n] <0:            
                                                print(n)
                                                PeakIdxR= n
                                                break
    for n in range (CIdx, int(CIdx-50), -1):
        if LONew[n+1]-LONew[n] <0:
            if LONew[n+2]-LONew[n] <0:
                if LONew[n+3]-LONew[n] <0:            
                    if LONew[n+4]-LONew[n] <0:            
                        if LONew[n+5]-LONew[n] <0:            
                            if LONew[n-1]-LONew[n] <0:            
                                if LONew[n-2]-LONew[n] <0:            
                                    if LONew[n-3]-LONew[n] <0:            
                                        if LONew[n-4]-LONew[n] <0:            
                                            if LONew[n-5]-LONew[n] <0:            
                                                print(n)
                                                PeakIdxL= n
                                                break
    CIdx= int((PeakIdxR+PeakIdxL)/2)
    PeakIdxArray_R[count][0]= Xidx
    PeakIdxArray_R[count][1]= PeakIdxR
    PeakIdxArray_L[count][0]= Xidx
    PeakIdxArray_L[count][1]= PeakIdxL
    count= count+1

#%%
plt.plot(PeakIdxArray_L[:, 0], PeakIdxArray_L[:, 1], '.')
plt.plot(PeakIdxArray_R[:, 0], PeakIdxArray_R[:, 1], '.')
plt.xlabel('x (pixel#)')
plt.ylabel('y (pixel#)')

#%%
PeakIdxArray= PeakIdxArray_L
#%%
NewArray1= np.zeros((1, 2))
count=0
for idx in range (int(PeakIdxArray.shape[0]-1)):
    if abs(PeakIdxArray[idx+1, 1]-PeakIdxArray[idx, 1])<5:
        NewArray1[count][0]= PeakIdxArray[idx+1, 0]
        NewArray1[count][1]= PeakIdxArray[idx+1, 1]
        if idx<int(PeakIdxArray.shape[0]-2):
            NewArray1= np.append(NewArray1, [[0, 0]], axis=0)
        count= count+1

#Calculate slopes
count=0
idx=0
model= LinearRegression()
Slope= np.zeros((1,2))   
Intercept= np.zeros((1,2))   
while idx<int(NewArray1.shape[0]-10):
    Set=NewArray1[idx:idx+5, :]
    xdata=np.reshape(Set[:, 0], (int(Set[:, 0].shape[0]), 1))
    ydata=np.reshape(Set[:, 1], (int(Set[:, 1].shape[0]), 1))
    model.fit(xdata, ydata)
    slope= model.coef_[0][0]
    Slope[count][0]= idx
    Slope[count][1]= slope
    intercept= model.intercept_[0]
    Intercept[count][0]= idx
    Intercept[count][1]= intercept
    if idx<int(NewArray1.shape[0]-10-1):
        Slope= np.append(Slope, [[0, 0]], axis=0)
        Intercept= np.append(Intercept, [[0, 0]], axis=0)
    count= count+1
    idx= idx+1

hist, binEdges = np.histogram(Slope[:, 1])
argmax= np.argmax(hist)
IdealSlope= (binEdges[argmax]+binEdges[argmax+1])/2

#Find killing index
avgSlope= np.mean(Slope[:, 1])
avgIntercept= np.mean(Intercept[:, 1])
stdSlope= statistics.stdev(Slope[:, 1])
DelIdxS= np.zeros(1)
count= 0
for idx in range(int(Slope.shape[0])):
    if abs(Slope[idx, 1]-IdealSlope)>stdSlope*0.8:
        DelIdxS[count]= idx
        DelIdxS= np.append(DelIdxS, [0])
        count=count+1
DelIdxS= np.delete(DelIdxS, int(len(DelIdxS)-1))
DelIdxS= DelIdxS.astype(int)

Slope2 = np.delete(Slope[:, 1], DelIdxS)
Intercept2 = np.delete(Intercept[:, 1], DelIdxS)

avgSlope= np.mean(Slope2[:])
avgIntercept= np.mean(Intercept2[:])

killIdx= np.zeros(1)
mark=1
count=0
NewArray2= np.zeros((1, 2))

for idx in range (int(NewArray1.shape[0])):
    m= NewArray1[idx, 0]
    n= NewArray1[idx, 1]
    distance= abs(-avgSlope*m+n-avgIntercept)/np.sqrt((-avgSlope)**2+1**2)
    
    if distance< 5:
        if idx== int(NewArray1.shape[0]-2):
            break
        NewArray2[count][0]= NewArray1[idx, 0]
        NewArray2[count][1]= NewArray1[idx, 1]
        NewArray2= np.append(NewArray2, [[0, 0]], axis=0)
        count= count+1
NewArray2= np.delete(NewArray2, int(len(NewArray2)-1), 0)

xdata=np.reshape(NewArray2[:, 0], (int(NewArray2[:, 0].shape[0]), 1))
ydata=np.reshape(NewArray2[:, 1], (int(NewArray2[:, 1].shape[0]), 1))
model.fit(xdata, ydata)
Finalslope= model.coef_[0][0]
Finalintercept= model.intercept_[0]

fitfn= lambda x, S, I: S*x+I
#%%
xNew= np.linspace(np.amin(xdata), np.amax(xdata), 1000)
plt.plot(PeakIdxArray_L[:, 0], PeakIdxArray_L[:, 1], '.')
plt.plot(PeakIdxArray_R[:, 0], PeakIdxArray_R[:, 1], '.')
plt.plot(xdata, ydata, '.', label= 'Data')
plt.plot(xNew, fitfn(xNew, Finalslope, Finalintercept), '-', label='fitting')
plt.xlabel('x (pixel#)')
plt.ylabel('y (pixel#)')

#%%
np.save('Data28_L', [Finalslope, Finalintercept])