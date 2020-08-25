#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:03:48 2020

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:49:21 2019

Plot 6-11-19 data for laser power fitting

Results:  While Gaussian is okay, DoubleTanh is the best...again

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import ThreeDimensionAnalysis as ThrDim

timing=[50 ,60 ,70 ,80 ,90 ,100 ,105 ,110 ,115 ,120 ,125 ,130 ,135 ,140 ,145 ,150 ,170 ,190 ,210 ,230 ,250]
energy=[518,512,496,467,426,376 ,347 ,313 ,281.2,241.8,212.7,176.7,147.5,119.7,95.6,77.2,33.9,20.73,19.99,15.99,16.37]

timing = [0,1,2,3,4,5,6,7,8,9,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,21,22,23,24,25,27,29,31,33,35]
energy = [538,540,540,540,538,538,530,522,508,499,479,467,451,442,409,391,379,355,325,312,287,255,226.5,206.5,189.9,180.4,150.1,136.3,110.4,96.7,79.0,60.1,49.6,39.3,29.5,24.2,18.4,15.7,14.8,14.59,14.74]
nominal = 551

timing=np.array(timing)
energy=np.array(energy)

pfit_power = 6
supress = True

def GaussianOffset_setcenter(p, x):
    return  max(energy)*(np.exp(-.5*np.square(x-min(timing))/np.square(p[0]))+p[1]/max(energy))

def SingleTanh(p,x):
    return (.5 - .5*np.tanh((x-p[0])/p[1]) + p[3]/p[2]) * p[2]

gfit = ThrDim.FitDataSomething(energy,timing,GaussianOffset_setcenter,guess=[50,16],log=True,supress=supress)
tfit = ThrDim.FitDataSomething(energy,timing,ThrDim.DoubleTanhOffset,guess=[50,100,500,12],log=True,supress=supress)
pfit = np.polyfit(timing,energy,pfit_power)
sfit = ThrDim.FitDataSomething(energy,timing,SingleTanh,guess=[50,100,500,12],log=True,supress=False)

tim_arr=np.linspace(min(timing),max(timing),200)
enp_arr=np.zeros(200)
for i in range(pfit_power+1):
    enp_arr = enp_arr + pfit[i]*np.power(tim_arr,pfit_power-i)

gfit_full=[max(energy),gfit[0],min(timing),gfit[1]]
"""
plt.scatter(timing,energy,label='Data')
plt.plot(tim_arr,enp_arr,label='Polyfit')
plt.plot(tim_arr,ThrDim.GaussianOffset(gfit_full,tim_arr),label='Gaussian')
plt.plot(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("All fits on linear scale")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel(r'$\mathrm{Timing \ Delay \ }[\mu s]$')
plt.grid(); plt.legend(); plt.show()

plt.scatter(timing,energy,label='Data')
plt.semilogy(tim_arr,enp_arr,label='Polyfit')
plt.plot(tim_arr,ThrDim.GaussianOffset(gfit_full,tim_arr),label='Gaussian')
plt.plot(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("All fits on log scale")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel(r'$\mathrm{Timing \ Delay \ }[\mu s]$')
plt.grid(); plt.legend(); plt.show()

plt.scatter(timing,energy,label='Data')
#plt.semilogy(tim_arr,enp_arr,label='Polyfit')
#plt.plot(tim_arr,ThrDim.GaussianOffset(gfit_full,tim_arr),label='Gaussian')
plt.semilogy(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("DoubleTanh on log scale")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel(r'$\mathrm{Timing \ Delay \ }[\mu s]$')
plt.grid(); plt.legend(); plt.show()

plt.scatter(timing,energy,label='Data')
#plt.semilogy(tim_arr,enp_arr,label='Polyfit')
plt.semilogy(tim_arr,ThrDim.GaussianOffset(gfit_full,tim_arr),label='Gaussian')
#plt.semilogy(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("Gaussian on log scale")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel(r'$\mathrm{Timing \ Delay \ }[\mu s]$')
plt.grid(); plt.legend(); plt.show()
"""
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)

fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(1, 1, 1)
plt.scatter(timing+nominal,energy,label='3-11-20 Measurements')
#plt.semilogy(tim_arr,enp_arr,label='Polyfit')
plt.plot(tim_arr+nominal,SingleTanh(sfit,tim_arr),label='Tanh Fit')
#plt.semilogy(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("Laser Energy vs Timing Delay")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel("'Sync Qsw #2' Timing Delay "+r'$[ns]$')
ax.grid(b=True, which='major', color='k', linestyle='-')

xminorticks = np.arange(0,36,1)+nominal
yminorticks = np.arange(0,575,25)
ax.set_xticks(xminorticks, minor=True)
ax.set_yticks(yminorticks, minor=True)
ax.grid(b=True, which='minor', color='r', linestyle=':')
plt.xlim([nominal-1,nominal+1+timing[-1]])
#plt.grid(); 
plt.legend(); plt.show()
"""
plt.scatter(timing,energy,label='Data')
#plt.semilogy(tim_arr,enp_arr,label='Polyfit')
plt.semilogy(tim_arr,SingleTanh(sfit,tim_arr),label='SingleTanh')
#plt.semilogy(tim_arr,ThrDim.DoubleTanhOffset(tfit,tim_arr),label='DoubleTanh')
plt.title("Single Tanh on log scale")
plt.ylabel(r'$\mathrm{Energy \ }[mJ]$')
plt.xlabel(r'$\mathrm{Timing \ Delay \ }[\mu s]$')
plt.grid(); plt.legend(); plt.show()
"""