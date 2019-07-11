#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:03:23 2019

Plot 6-11-19 data for vacuum chamber outgassing and window ablation, which was
suspected due to the long axicon hitting the exit window of Chamber B

Results:  No significant ablation, and outgassing rate is 

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
from modules import ThreeDimensionAnalysis as ThrDim

minute=[0,1,2,4,7,10]
pgas0=[4.20e-05,2.01e-04,3.34e-04,6.37e-04,1.19e-03,1.68e-03]#Just outgassing, no laser
pgas1=[5.54e-05,2.26e-04,3.67e-04,7.15e-04,1.27e-03,1.76e-03]#Outgassing and half-power laser
pgas2=[5.54e-05,2.21e-04,3.58e-04,7.15e-04,1.27e-03,1.76e-03]#Outgassing and full-power laser

minute=np.array(minute); pgas0=np.array(pgas0); pgas1=np.array(pgas1); pgas2=np.array(pgas2)
minarr=np.linspace(0,10,100)

def ExponentialIncreasingDecay(p, x):
    return p[0]*(1-np.exp(-p[1]*x))

def InvExpDec(p,pr):
    return -1/p[1]*np.log(1-pr/p[0])

def DoubleTime(p,pr):
    return 1/p[1]*np.log((p[0]-pr)/(p[0]-2*pr))

efit0 = ThrDim.FitDataSomething(pgas0,minute,ExponentialIncreasingDecay,guess=[pgas0[0],1e-01],log=True)

plt.semilogy(minarr,ExponentialIncreasingDecay(efit0,minarr),label='Outgas')
plt.scatter(minute,pgas0,label='Data')
plt.title("Outgas Rate")
plt.xlabel("Time [min]")
plt.ylabel("Pressure [mbar]")
plt.grid(); plt.legend(); plt.show()

print("Exponential Fit:  p(t)=a(1-e^(-b*t))")
print("  a = "+str(efit0[0]))
print("  b = "+str(efit0[1]))

efit1 = ThrDim.FitDataSomething(pgas1,minute,ExponentialIncreasingDecay,guess=[pgas1[0],-2],log=True,supress=True)
efit2 = ThrDim.FitDataSomething(pgas2,minute,ExponentialIncreasingDecay,guess=[pgas2[0],-2],log=True,supress=True)

plt.semilogy(minarr,ExponentialIncreasingDecay(efit0,minarr),label='Outgas')
plt.semilogy(minarr,ExponentialIncreasingDecay(efit1,minarr),label='Half Power')
plt.semilogy(minarr,ExponentialIncreasingDecay(efit2,minarr),label='Full Power')
plt.title("Ablation Comparisons")
plt.xlabel("Time [min]")
plt.ylabel("Pressure [mbar]")
plt.grid(); plt.legend(); plt.show()

print("Time to go from 1e-04 to 2e-04 mbar:")
print(DoubleTime(efit0,1e-04),"min")
print("Time to go from 1e-03 to 2e-03 mbar:")
print(DoubleTime(efit0,1e-03),"min")
print("Time to go from 1e-02 to 2e-02 mbar:")
print(DoubleTime(efit0,1e-02),"min")
print("Time to go from 1e-01 to 2e-01 mbar:")
print(DoubleTime(efit0,1e-01),"min")

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)

plt.figure(figsize=(20,10))
pr_arr = np.logspace(-5,-1,100)
plt.loglog(pr_arr, DoubleTime(efit0,pr_arr))
plt.title("Outgas Time for Double Pressure")
plt.xlabel("Starting Pressure [mbar]")
plt.ylabel("Time to Double [min]")
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='r', linestyle=':')
plt.minorticks_on()
plt.show()