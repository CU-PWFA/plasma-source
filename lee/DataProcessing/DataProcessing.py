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
#%%14
test= np.array(Image.open('2009210111/19423601_2009210111_0007.tiff'))
#%%
BG0= np.array(Image.open('2009210112/19423601_2009210112_0002.tiff'))
#%%32
test= np.array(Image.open('2006260032/19423601_2006260032_0001.tiff'))
#%%
BG0= np.array(Image.open('2006260033/19423601_2006260033_0002.tiff'))

#%%38
test= np.array(Image.open('2006260038/19423601_2006260038_0001.tiff'))
#%%
BG0= np.array(Image.open('2006260039/19423601_2006260039_0001.tiff'))

#%%
BG= BG0
BGfreq= fftshift(fft2(BG))
testfreq= fftshift(fft2(test))
Newfreq= testfreq-BGfreq
NewS= abs(ifft2(fftshift(Newfreq)))
#%%
#SignalArea= NewS[500:800, 0:1600]
SignalArea= NewS[200:500, 600:850]
Newsignal= gaussian_filter(SignalArea, sigma=5)
NewsignalBG= gaussian_filter(SignalArea, sigma=20)
NS= Newsignal-NewsignalBG
#%%
plt.figure(2)
plt.pcolormesh(NS)
plt.axis('scaled')
plt.colorbar()
#plt.title('BG')
#%%
#NS= BG0[200:500, 600:850]
NS= test[200:500, 600:850]
Newsignal= gaussian_filter(NS, sigma=2)

#%%
plt.plot(Newsignal[:, 120])
#%%
IntCenter= 200
CIdx= IntCenter
count=0
IdxStart= 500
IdxEnd= 1550
DataLength= int((IdxEnd-IdxStart)/10)
PeakIdxRArray= np.zeros((DataLength, 2))
PeakIdxLArray= np.zeros((DataLength, 2))

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
    PeakIdxRArray[count][0]= Xidx
    PeakIdxRArray[count][1]= PeakIdxR
    PeakIdxLArray[count][0]= Xidx
    PeakIdxLArray[count][1]= PeakIdxL
    count= count+1

#%%
plt.plot(PeakIdxLArray[:, 0], PeakIdxLArray[:, 1], '.')
plt.plot(PeakIdxRArray[:, 0], PeakIdxRArray[:, 1], '.')
plt.xlabel('x (pixel#)')
plt.ylabel('y (pixel#)')


#%%
NewArray1= np.zeros((1, 2))
count=0
for idx in range (int(PeakIdxArray.shape[0]-1)):
    if abs(PeakIdxArray[idx+1, 1]-PeakIdxArray[idx, 1])<5:
        NewArray1[count][0]= PeakIdxArray[idx+1, 0]
        NewArray1[count][1]= PeakIdxArray[idx+1, 1]
        NewArray1= np.append(NewArray1, [[0, 0]], axis=0)
        count= count+1

#%%
#Calculate slopes
count=0
idx=0
model= LinearRegression()
Slope= np.zeros((1,2))   
while idx<int(NewArray1.shape[0]-10):
    Set=NewArray1[idx:idx+10, :]
    xdata=np.reshape(Set[:, 0], (int(Set[:, 0].shape[0]), 1))
    ydata=np.reshape(Set[:, 1], (int(Set[:, 1].shape[0]), 1))
    model.fit(xdata, ydata)
    slope= model.coef_[0][0]
    Slope[count][0]= idx
    Slope[count][1]= slope
    Slope= np.append(Slope, [[0, 0]], axis=0)
    count= count+1
    idx= idx+1
#Find killing index
mean= np.mean(Slope[:, 1])
std= statistics.stdev(Slope[:, 1])
killIdx= np.zeros(1)
mark=1
for idx in range(Slope.shape[0]):
    if abs(Slope[idx, 1]-mean)>std*1:
        print(idx)
        for i in range (0, 10):
            for j in range (len(killIdx)):
                if idx+i==killIdx[j]:
                    mark= 0
            if mark==1:
                killIdx[int(len(killIdx)-1)]= idx+i
                killIdx= np.append(killIdx, [0])
            mark=1
        print(killIdx)
    
#kill those  index
NewArray2= np.zeros((1, 2))
count=0
killList= killIdx.tolist()
for idx in range (5):
    killList.pop(0)

for idx in range (int(NewArray1.shape[0]-1)):
    if idx!=killList[0]:
        NewArray2[count][0]= NewArray1[idx, 0]
        NewArray2[count][1]= NewArray1[idx, 1]
        while int(NewArray2.shape[0])< int(killIdx.shape[0]):
            NewArray2= np.append(NewArray2, [[0, 0]], axis=0)
            break
        count= count+1
    else:
        killList.pop(0)
        print(killList)
        print(idx)
#%%
#%%
#Calculate slopes
count=0
idx=0
model= LinearRegression()
Slope= np.zeros((1,2))   
Intercept= np.zeros((1,2))   
while idx<int(NewArray1.shape[0]-10):
    Set=NewArray1[idx:idx+10, :]
    xdata=np.reshape(Set[:, 0], (int(Set[:, 0].shape[0]), 1))
    ydata=np.reshape(Set[:, 1], (int(Set[:, 1].shape[0]), 1))
    model.fit(xdata, ydata)
    slope= model.coef_[0][0]
    Slope[count][0]= idx
    Slope[count][1]= slope
    Slope= np.append(Slope, [[0, 0]], axis=0)
    intercept= model.intercept_[0]
    Intercept[count][0]= idx
    Intercept[count][1]= intercept
    Intercept= np.append(Intercept, [[0, 0]], axis=0)
    count= count+1
    idx= idx+1
#Find killing index
avgSlope= np.mean(Slope[:, 1])
avgIntercept= np.mean(Intercept[:, 1])
stdSlope= statistics.stdev(Slope[:, 1])
DelIdxS= np.zeros(1)
count= 0
for idx in range(int(Slope.shape[0])):
    if abs(Slope[idx, 1]-avgSlope)>stdSlope*2:
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
    if distance< 10:
        NewArray2[count][0]= NewArray1[idx, 0]
        NewArray2[count][1]= NewArray1[idx, 1]
        NewArray2= np.append(NewArray2, [[0, 0]], axis=0)
        count= count+1
#%%
NewLArray= np.zeros((1, 2))
count=0
for idx in range (int(PeakIdxLArray.shape[0]-1)):
    if abs(PeakIdxLArray[idx+1, 1]-PeakIdxLArray[idx, 1])<5:
        NewLArray[count][0]= PeakIdxLArray[idx+1, 0]
        NewLArray[count][1]= PeakIdxLArray[idx+1, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1

#%%
PeakIdxLArray= NewLArray

#%%
#Calculate slopes
count=0
idx=0
model= LinearRegression()
Slope= np.zeros((1,2))   
while idx<int(PeakIdxLArray.shape[0]-10):
    Set=PeakIdxLArray[idx:idx+10, :]
    xdata=np.reshape(Set[:, 0], (int(Set[:, 0].shape[0]), 1))
    ydata=np.reshape(Set[:, 1], (int(Set[:, 1].shape[0]), 1))
    model.fit(xdata, ydata)
    slope= model.coef_[0][0]
    Slope[count][0]= idx
    Slope[count][1]= slope
    Slope= np.append(Slope, [[0, 0]], axis=0)
    count= count+1
    idx= idx+1
#Find killing index
mean= np.mean(Slope[:, 1])
std= statistics.stdev(Slope[:, 1])
killIdx= np.zeros(1)
mark=1
for idx in range(Slope.shape[0]):
    if abs(Slope[idx, 1]-mean)>std*1:
        print(idx)
        for i in range (0, 10):
            for j in range (len(killIdx)):
                if idx+i==killIdx[j]:
                    mark= 0
            if mark==1:
                killIdx[int(len(killIdx)-1)]= idx+i
                killIdx= np.append(killIdx, [0])
            mark=1
        print(killIdx)
    
#kill those  index
NewLArray= np.zeros((1, 2))
count=0
killList= killIdx.tolist()
for idx in range (5):
    killList.pop(0)

for idx in range (int(PeakIdxLArray.shape[0]-1)):
    if idx!=killList[0]:
        NewLArray[count][0]= PeakIdxLArray[idx, 0]
        NewLArray[count][1]= PeakIdxLArray[idx, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1
    else:
        killList.pop(0)
        print(killList)
        print(idx)
#%%
plt.plot(NewLArray[:,0], NewLArray[:,1], '.')

#%%
x= PeakIdxLArray[:, 0]
y= PeakIdxLArray[:, 1]
fitfn= lambda p, x: p[0]*x+p[1]

errfunc = lambda p, x, y: fitfn(p, x) - y;
p0= [-0.1, 10]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(x, y, '.', label= 'plasma data')
plt.plot(xNew, fitfn(p1, xNew), '-', label='fitting')
plt.legend()

#%%
NewLArray= np.zeros((1, 2))
count=0
idx=0
ywidth= PeakIdxLArray[1, 0]-PeakIdxLArray[0, 0]
while idx<int(PeakIdxLArray.shape[0]-1):    
    if abs((PeakIdxLArray[idx+1, 1]-PeakIdxLArray[idx, 1])/ywidth)-abs(p1[0])<0.05:
        NewLArray[count][0]= PeakIdxLArray[idx+1, 0]
        NewLArray[count][1]= PeakIdxLArray[idx+1, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1
        idx= idx+1
    else:
        fixIdx= idx
        for nextIdx in range (int(PeakIdxLArray.shape[0])):
            if abs((PeakIdxLArray[idx+nextIdx, 1]-PeakIdxLArray[fixIdx, 1])/ywidth)-abs(p1[0])>0.05:
                idx= idx+1
            else:
                break
            break
        idx=idx+1
    print(idx)
#%%
NewLArray= np.zeros((1, 2))
count=0
idx=0
    
while idx<int(PeakIdxLArray.shape[0]-1):    
    if abs(PeakIdxLArray[idx+1, 1]-PeakIdxLArray[idx, 1])<5:
        NewLArray[count][0]= PeakIdxLArray[idx+1, 0]
        NewLArray[count][1]= PeakIdxLArray[idx+1, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1
        idx= idx+1
    else:
        fixIdx= idx
        while idx<int(PeakIdxLArray.shape[0]-1):    
            if abs(PeakIdxLArray[idx+1, 1]-PeakIdxLArray[fixIdx, 1])>5:
                idx= idx+1
            else:
                break
        if idx>=int(PeakIdxLArray.shape[0]-1):
            break
        else:
            continue
        NewLArray[count][0]= PeakIdxLArray[idx+1, 0]
        NewLArray[count][1]= PeakIdxLArray[idx+1, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1
        idx=idx+1
    print(idx)