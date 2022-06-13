#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 21:05:53 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.interpolate import interp1d
from scipy import stats

#%%
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res
def FWHM(X,Y):
    half_max = max(Y) / 2.
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    return X[right_idx] - X[left_idx] 
#%%
Eng= [2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 3.375, 3.625, 3.875, 4.125, 4.375, 4.625]
#Eng= [3.5, 4, 4.5, 5, 5.5]
CenPix= 89
Gwidth= np.empty(len(Eng))
Gpower= np.empty(len(Eng))
fwhm= np.empty(len(Eng))
PlsEng= np.empty(len(Eng))
basePath= ''
count= 0
for energy in Eng:

    ne= np.load(basePath+'ne_'+str(energy)+'.npy')
    ne0= np.load(basePath+'ne0_'+str(energy)+'.npy')
    ext= np.load(basePath+'ext_'+str(energy)+'.npy')
    IntPlsE= np.load(basePath+'IntplsEn_'+str(energy)+'.npy')
    tgrid= np.linspace(ext[0], ext[1], ne.shape[0])
    #xgrid= np.linspace(ext[2], ext[3], ne.shape[1])
    xgrid= np.linspace(ext[2]-(2*ext[3])/ne.shape[1]/2, ext[3]-(2*ext[3])/ne.shape[1]/2, ne.shape[1])
    x= xgrid[1004: 1044]
    y= ne[72, 1004: 1044]
    
#    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**15)
    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**p[2])
    errfunc = lambda p, x, y: fitfn(p, x) - y
    
#    p0= [1e17, 10]
    p0= [1e17, 10, 2]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
    print(p1)
    
    xNew= np.linspace(np.amin(x), np.amax(x), 1000)
    plt.plot(x, y, '.')
    plt.plot(xNew, fitfn(p1, xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
    plt.xlabel('r ($\mu m$)')
    plt.ylabel('Density ($cm^{-3}$)')
    plt.legend()
    
    fwhm[count]= FWHM(xNew, fitfn(p1, xNew))
    print(fwhm)
    
    Gwidth[count]= p1[1]
    Gpower[count]= p1[2]
    PlsEng[count]= IntPlsE
    count= count+1
    
#%%
xth= np.zeros(100)+222.5
yth= np.linspace(20, 125, len(xth))
x= PlsEng
y= fwhm
y1= Gwidth*2

#fitfn= lambda p, x: p[0]*np.log(p[1]*x)+p[2]
fitfn= lambda p, x: p[0]*np.arctan(p[1]*x)+p[2]
errfunc = lambda p, x, y: fitfn(p, x) - y

p0= [3, 10, 2]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)

fitfn1= lambda p, x: p[0]*np.arctan(p[1]*x)+p[2]
errfunc1 = lambda p, x, y1: fitfn1(p, x) - y1

p2, success = optimize.leastsq(errfunc1, p0[:], args=(x, y1))#, epsfcn= 2.5e-34);
print(p2)
#%%
xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(PlsEng, y1, 'o', label='Gaussian Width')
plt.plot(PlsEng, y, 'o', label= 'FWHM')
plt.plot(xNew, fitfn(p1, xNew), '-')
plt.plot(xNew, fitfn1(p2, xNew), '-')
plt.xlabel('Initial Pluse Energy (mJ)')
plt.ylabel('Width ($\mu m$)')
plt.plot(xth, yth, label='Full inozation')
plt.legend()
#%%
plt.figure(3)
plt.plot(PlsEng, Gwidth*2, 'o', label='Gaussian Full Width')
plt.plot(PlsEng, fwhm, 'o', label= 'FWHM')
plt.plot(PlsEng, dirfwhm, 'o', label= 'dirFWHM')
plt.xlabel('Initial Pluse Energy (mJ)')
plt.ylabel('Width ($\mu m$)')
plt.legend()
#plt.title('power= 5')
#%%
#xth= np.zeros(100)+222.5
xth= np.zeros(20)+253
yth= np.linspace(0, 20, len(xth))
plt.figure(4)
plt.plot(PlsEng, Gpower*2, 'o', label='Gaussian Power')
plt.plot(xth, yth, label='Full inozation')
plt.xlabel('Initial Pluse Energy (mJ)')
plt.ylabel('Power')
plt.legend()

#%%
x= np.array(Eng)
y= PlsEng*1e-3

fitfn= lambda p, x: p[0]+p[1]*x+p[2]*x**2
errfunc = lambda p, x, y: fitfn(p, x) - y

p0= [1, 2, 3]
p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
print(p1)

xNew= np.linspace(np.amin(x), np.amax(x), 1000)
plt.plot(x, y, 'o')#, label= 'FWHM')
plt.plot(xNew, fitfn(p1, xNew), '-')
plt.xlabel('E Field (1e-9 V/m)')
plt.ylabel('Initial Energy (J)')
#plt.legend()
#%%
#Eng= [2.25]
Eng= [2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 3.375, 3.625, 3.875, 4.125, 4.375, 4.625]
#Eng= [3.5, 4, 4.5, 5, 5.5]
CenPix= 89
Gwidth= np.empty(len(Eng))
Gpower= np.empty(len(Eng))
fwhm= np.empty(len(Eng))
PlsEng= np.empty(len(Eng))
Rsqrt= np.empty(len(Eng))
KStest= np.empty(len(Eng))

basePath= ''
count= 0

for energy in Eng:

    ne= np.load(basePath+'ne_'+str(energy)+'.npy')
    ne0= np.load(basePath+'ne0_'+str(energy)+'.npy')
    ext= np.load(basePath+'ext_'+str(energy)+'.npy')
    IntPlsE= np.load(basePath+'IntplsEn_'+str(energy)+'.npy')
    tgrid= np.linspace(ext[0], ext[1], ne.shape[0])
    #xgrid= np.linspace(ext[2], ext[3], ne.shape[1])
    xgrid= np.linspace(ext[2]-(2*ext[3])/ne.shape[1]/2, ext[3]-(2*ext[3])/ne.shape[1]/2, ne.shape[1])
    x= xgrid[1014: 1034]
    y= ne[CenPix, 1014: 1034]
    
#    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**15)
    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**p[2])
    errfunc = lambda p, x, y: fitfn(p, x) - y
    
#    p0= [1e17, 10]
    p0= [1e17, 10, 2]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
    print(p1)
    
    xNew= np.linspace(np.amin(x), np.amax(x), 1000)
    
#    plt.plot(x, y, '.')
#    plt.plot(xNew, fitfn(p1, xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
#    plt.xlabel('r ($\mu m$)')
#    plt.ylabel('Density ($cm^{-3}$)')
#    plt.legend()
    yavg= sum(y)/len(y)
    SSt= 0
    SSr= 0
    for i in range (len(y)):
        fi= fitfn(p1, x[i])
        SSt= SSt+ (y[i]- yavg)**2
        SSr= SSr+ (y[i]- fi)**2
#    R2= (SSr)
    R2= 1-(SSr/SSt)
    print(R2)
    KS= stats.kstest(y, fitfn(p1, xNew))
    print(KS)

    fwhm[count]= FWHM(xNew, fitfn(p1, xNew))
    Gwidth[count]= p1[1]
    Gpower[count]= p1[2]
    PlsEng[count]= IntPlsE
    Rsqrt[count]= R2
    KStest[count]= KS [0]
    count= count+1
#%%
plt.figure(1)
plt.plot(PlsEng, KStest, '.')
plt.title('KStest')
#plt.plot(PlsEng, Rsqrt, '.')
#plt.title('R2')
#%%
Eng= [2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 3.375, 3.625, 3.875, 4.125, 4.375, 4.625]
#Eng= [3.875]
#Eng= [3.875, 4, 4.125]
CenPix= 89
Gwidth= np.empty(len(Eng))
Gpower= np.empty(len(Eng))
fwhm= np.empty(len(Eng))
PlsEng= np.empty(len(Eng))
basePath= ''
count= 0
dirfwhm= np.zeros(len(Eng))
for energy in Eng:

    ne= np.load(basePath+'ne_'+str(energy)+'.npy')
    ne0= np.load(basePath+'ne0_'+str(energy)+'.npy')
    ext= np.load(basePath+'ext_'+str(energy)+'.npy')
    IntPlsE= np.load(basePath+'IntplsEn_'+str(energy)+'.npy')
    tgrid= np.linspace(ext[0], ext[1], ne.shape[0])
    #xgrid= np.linspace(ext[2], ext[3], ne.shape[1])
    xgrid= np.linspace(ext[2]-(2*ext[3])/ne.shape[1]/2, ext[3]-(2*ext[3])/ne.shape[1]/2, ne.shape[1])
    x= xgrid[1014: 1034]
    y= ne[CenPix, 1014: 1034]
    
    f= interp1d(x, y)
#    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**5)
    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**p[2])
    errfunc = lambda p, x, y: fitfn(p, x) - y
    
#    p0= [1e17, 50]
    p0= [1e17, 10, 2]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
    print(p1)
    
    xNew= np.linspace(np.amin(x), np.amax(x), 1000)
    plt.plot(x, y, '.')
    plt.plot(xNew, fitfn(p1, xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
    plt.plot(xNew, f(xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
    plt.xlabel('r ($\mu m$)')
    plt.ylabel('Density ($cm^{-3}$)')
    plt.legend()
    
    fwhm[count]= FWHM(xNew, fitfn(p1, xNew))
    print(fwhm)
    
    dirfwhm[count]= FWHM(xNew, f(xNew))
    Gwidth[count]= p1[1]
#    Gpower[count]= p1[2]
    PlsEng[count]= IntPlsE
    count= count+1
#%%
FWHMError= np.zeros(39)
for power in range (1, 40):
    Eng= [2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 3.375, 3.625, 3.875, 4.125, 4.375, 4.625]
    CenPix= 89
    Gwidth= np.empty(len(Eng))
    Gpower= np.empty(len(Eng))
    fwhm= np.empty(len(Eng))
    PlsEng= np.empty(len(Eng))
    basePath= ''
    count= 0
    dirfwhm= np.zeros(len(Eng))
    for energy in Eng:
    
        ne= np.load(basePath+'ne_'+str(energy)+'.npy')
        ne0= np.load(basePath+'ne0_'+str(energy)+'.npy')
        ext= np.load(basePath+'ext_'+str(energy)+'.npy')
        IntPlsE= np.load(basePath+'IntplsEn_'+str(energy)+'.npy')
        tgrid= np.linspace(ext[0], ext[1], ne.shape[0])
        #xgrid= np.linspace(ext[2], ext[3], ne.shape[1])
        xgrid= np.linspace(ext[2]-(2*ext[3])/ne.shape[1]/2, ext[3]-(2*ext[3])/ne.shape[1]/2, ne.shape[1])
        x= xgrid[1014: 1034]
        y= ne[CenPix, 1014: 1034]
        
        f= interp1d(x, y)
        fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**power)
    #    fitfn= lambda p, x: p[0]*np.exp(-(x**2/(2*p[1]**2))**p[2])
        errfunc = lambda p, x, y: fitfn(p, x) - y
        
        p0= [1e17, 50]
    #    p0= [1e17, 10, 2]
        p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))#, epsfcn= 2.5e-34);
#        print(p1)
        
        xNew= np.linspace(np.amin(x), np.amax(x), 1000)
        plt.plot(x, y, '.')
        plt.plot(xNew, fitfn(p1, xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
        plt.plot(xNew, f(xNew), '-', label='PulseEnergy='+str(np.round(IntPlsE, 2)))
        plt.xlabel('r ($\mu m$)')
        plt.ylabel('Density ($cm^{-3}$)')
        plt.legend()
        
        fwhm[count]= FWHM(xNew, fitfn(p1, xNew))
#        print(fwhm)
        
        dirfwhm[count]= FWHM(xNew, f(xNew))
        Gwidth[count]= p1[1]
    #    Gpower[count]= p1[2]
        PlsEng[count]= IntPlsE
        count= count+1
    
    fwhmError= sum((dirfwhm-fwhm)**2)
    FWHMError[power]= fwhmError
    print(power, fwhmError)
#%%
plt.plot(np.linspace(1, 35, 35), FWHMError[0:35], 'o')
plt.xlabel('power')
plt.ylabel('fwhm Error')