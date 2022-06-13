#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 18:13:26 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
from PIL import Image
from scipy.ndimage import gaussian_filter
#%%04
test= np.array(Image.open('19423601_2006260004_0001.tiff'))
#%%
BG2= np.array(Image.open('19423601_2006260005_0002.tiff'))
BG3= np.array(Image.open('19423601_2006260005_0003.tiff'))
BG4= np.array(Image.open('19423601_2006260005_0004.tiff'))
BG5= np.array(Image.open('19423601_2006260005_0005.tiff'))
BG6= np.array(Image.open('19423601_2006260005_0006.tiff'))
BG7= np.array(Image.open('19423601_2006260005_0007.tiff'))
BG8= np.array(Image.open('19423601_2006260005_0008.tiff'))
#%%14
test= np.array(Image.open('19423601_2006260014_0001.tiff'))
#%%
BG0= np.array(Image.open('19423601_2006260015_0000.tiff'))
BG1= np.array(Image.open('19423601_2006260015_0001.tiff'))
BG2= np.array(Image.open('19423601_2006260015_0002.tiff'))
BG3= np.array(Image.open('19423601_2006260015_0003.tiff'))
BG4= np.array(Image.open('19423601_2006260015_0004.tiff'))
BG5= np.array(Image.open('19423601_2006260015_0005.tiff'))

#%%10
test= np.array(Image.open('19423601_2006260010_0000.tiff'))
#%%
BG0= np.array(Image.open('19423601_2006260011_0001.tiff'))
#%%20
test= np.array(Image.open('19423601_2006260020_0001.tiff'))
#%%
BG0= np.array(Image.open('19423601_2006260021_0000.tiff'))

#%%24
test= np.array(Image.open('19423601_2006260024_0002.tiff'))
#%%
BG0= np.array(Image.open('19423601_2006260025_0000.tiff'))

#%%26
test= np.array(Image.open('19423601_2006260026_0001.tiff'))
#%%
BG0= np.array(Image.open('19423601_2006260027_0001.tiff'))
#%%32
test= np.array(Image.open('2006260032/19423601_2006260032_0001.tiff'))
#%%
BG0= np.array(Image.open('2006260033/19423601_2006260033_0002.tiff'))

#%%38
test= np.array(Image.open('2006260038/19423601_2006260038_0001.tiff'))
#%%
BG0= np.array(Image.open('2006260039/19423601_2006260039_0001.tiff'))

#%%
Energy=1
BG= BG0
#NorE=On_4e17_500_62*Energy/sum(sum(On_4e17_500_62))*1e4
Nortest=test*Energy/sum(sum(test))*1e4
NorBG=BG*Energy/sum(sum(BG))*1e4

#%%
Norsignal= Nortest-NorBG
#Norsignal= Nortest/NorBG

#%%
BGfreq= fftshift(fft2(BG))
testfreq= fftshift(fft2(test))
Newfreq= testfreq-BGfreq
NewS= abs(ifft2(fftshift(Newfreq)))
#%%
def us(x):
    return (np.sign(x)+1)/2

def rect(x, halfwidth):
    return us(x+halfwidth)-us(x-halfwidth)

def Rect2D(Q, xlength, ylength):
    x=np.linspace(-xlength/2, xlength/2, xlength)
    y=np.linspace(-ylength/2, ylength/2, ylength)
    X, Y=np.meshgrid(x, y)
    R=np.sqrt(X**2+Y**2)
    return rect(R, Q)

def gaus2d(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))

def Rect2DXY(QX, QY, xlength, ylength):
    x=np.linspace(-xlength/2, xlength/2, xlength)
    y=np.linspace(-ylength/2, ylength/2, ylength)
    x, y = np.meshgrid(x, y) # get 2D variables instead of 1D
    z = gaus2d(x, y, 0, 0, QX, QY)
    return z

def HighPassFilter (Matrix, Q):
    NewMatrix=abs(ifft2(fftshift(fftshift(fft2(Matrix))*abs(Rect2D(Q, Matrix.shape[1], Matrix.shape[0])-1))))
    return NewMatrix

def LowPassFilter (Matrix, Q):
    NewMatrix=abs(ifft2(fftshift(fftshift(fft2(Matrix))*Rect2D(Q, Matrix.shape[1], Matrix.shape[0]))))
    return NewMatrix

def LowPassXYFilter (Matrix, QX, QY):
    NewMatrix=abs(ifft2(fftshift(fftshift(fft2(Matrix))*Rect2DXY(QX, QY, Matrix.shape[1], Matrix.shape[0]))))
    return NewMatrix

def HighPassXYFilter (Matrix, QX, QY):
    NewMatrix=abs(ifft2(fftshift(fftshift(fft2(Matrix))*(-1*Rect2DXY(QX, QY, Matrix.shape[1], Matrix.shape[0])+1))))
    return NewMatrix

#%%
plt.plot(abs(fftshift(fft2(Norsignal)))[:, 1024])
#%%
SignalArea= NewS[500:800, 400:2000]
SAfreq= fftshift(fft2(SignalArea))
testLP= LowPassXYFilter(SignalArea, 3, 8)
testLP2= HighPassXYFilter(testLP, 2, 5)
#%%
#%%
#Newsignal= gaussian_filter(Norsignal, sigma=5)
#Newsignal= gaussian_filter(SignalArea, sigma=5)
Newsignal= gaussian_filter(NewS, sigma=5)
#Newsignal= LowPassFilter(NewS, 200)
#Newsignal= LowPassFilter(SignalArea, 200)
#Newsignal= HighPassFilter(Newsignal, 4)
#Newsignal= LowPassFilter(Norsignal, 120)
#%%
NewsignalBG= gaussian_filter(SignalArea, sigma=20)
#%%
NS= Newsignal-NewsignalBG
#%%
LO=Newsignal[:, 600]
LOfreq= fftshift(fft(LO))
x= np.linspace(-LO.shape[0]/2, LO.shape[0]/2, LO.shape[0])

LOfreq2= fftshift(fft(LO))*rect(x, 2)
LONew= abs(ifft(fftshift(fftshift(fft(LO))*(-1*rect(x, 2)+1))))
#%%
plt.figure(9)
plt.plot(NS[:, 580])
#plt.plot(LO)

#%%
LOsave= np.zeros((Newsignal.shape[0], Newsignal.shape[1]))
for Xidx in range(Newsignal.shape[1]):
    LO= Newsignal[:, Xidx]
    LONew= abs(ifft(fftshift(fftshift(fft(LO))*(-1*rect(x, 2)+1))))
    LOsave[:, Xidx]= LONew
    plt.plot(LONew)

#%%
plt.figure(1)
plt.pcolormesh(LOsave)
plt.axis('scaled')    
#%%
plt.figure(6)
#plt.pcolormesh(NewsignalBG)
#plt.pcolormesh(Newsignal)
#plt.pcolormesh(NS)

#plt.pcolormesh(SignalArea)
#plt.pcolormesh(NewS)
plt.pcolormesh(NorBG)
#plt.pcolormesh(abs(SAfreq))
#plt.pcolormesh(testLP2)
#plt.pcolormesh(BG)
plt.axis('scaled')
plt.colorbar()
#%%
plt.figure(1)
plt.plot(LO)
#%%
np.save('Data14', Newsignal)
np.save('data14', Newsignal[:, 500])