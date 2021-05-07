#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 12:03:04 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

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
PlasmaAfterglow= np.zeros((1, 6))
count= 0
#Den= [2, 5, 8, 10, 40, 80]
Den= [2, 5, 8, 10, 20, 40, 50, 80]
Pow= [5]
PWidth= [30, 35, 40, 45, 50]

#%%
GPcount= 0
PWcount= 0

for density in Den:
    vars()['AgI_n'+str(density)]= np.zeros((len(Pow), len(PWidth)))
    vars()['AgP_n'+str(density)]= np.zeros((len(Pow), len(PWidth)))
    vars()['AgW_n'+str(density)]= np.zeros((len(Pow), len(PWidth)))
#    vars()['AgFWHM_n'+str(density)]= np.zeros((len(Pow), len(PWidth)))
    for power in Pow:
        for pw in PWidth:
            try:
                
                PhotonMap= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
                IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
                r= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
                xvar= r[int(len(r)/2):len(r)]
                yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
                fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])
    
                errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
                p0= [1e11, 20, 2];
                p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
                            
                vars()['AgI_n'+str(density)][GPcount][PWcount]= p1[0] #photon #
                vars()['AgP_n'+str(density)][GPcount][PWcount]= p1[2] #power
                vars()['AgW_n'+str(density)][GPcount][PWcount]= p1[1] #width
                
                PWcount= PWcount+1
            except:
                pass
        GPcount= GPcount+1
        PWcount= 0
    GPcount= 0
#%%
GPcount= 0
DENcount= 0

for pw in PWidth:
    vars()['AgI_pw'+str(pw)]= np.zeros((len(Pow), len(Den)))
    vars()['AgP_pw'+str(pw)]= np.zeros((len(Pow), len(Den)))
    vars()['AgW_pw'+str(pw)]= np.zeros((len(Pow), len(Den)))
    for power in Pow:
        for density in Den:
            
            PhotonMap= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
            IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
            r= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
            xvar= r[int(len(r)/2):len(r)]
            yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
            fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])

            errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
            p0= [1e11, 20, 2];
            p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
            print(power, GPcount)
            print(pw, DENcount)
            print(p1)
            vars()['AgI_pw'+str(pw)][GPcount][DENcount]= p1[0] #photon #
            vars()['AgP_pw'+str(pw)][GPcount][DENcount]= p1[2] #power
            vars()['AgW_pw'+str(pw)][GPcount][DENcount]= p1[1] #width
            
            DENcount= DENcount+1
        GPcount= GPcount+1
        DENcount= 0
    GPcount= 0
#%%
for den in Den:
    vars()['den'+str(den)]= np.zeros(len(Pow))
    vars()['den'+str(den)].fill(den)
#%%
FittingPara= []
for den in Den:
    
    xvar= np.array(PWidth)
    yvar= vars()['AgW_n'+str(den)][0, :]
    fitfn_All= lambda p, xvar: p[0]*xvar+p[1]
    
    errfunc = lambda p, xvar, yvar: fitfn_All(p, xvar)- yvar
    p0= [1, 1]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    print(p1)
    FittingPara.append([p1[0], p1[1]])

#NewX= np.linspace(np.amin(xvar), np.amax(xvar), 1000)
#plt.plot(xvar, yvar, 'o')
#plt.plot(NewX, fitfn(p1, NewX))
FittingPara= np.array(FittingPara)
#%%
xvar= np.array(Den[0:6])
yvar= FittingPara[0:6, 0]
fitfn= lambda p, xvar: p[0]*xvar+p[1]

errfunc = lambda p, xvar, yvar: fitfn(p, xvar)- yvar
p0= [1, 1]
p_slope, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
#%%
NewX= np.linspace(np.amin(Den), np.amax(Den), 1000)
plt.plot(Den, FittingPara[:, 0], 'o')
plt.plot(NewX, fitfn(p_slope, NewX))
plt.xlabel('Density $(1e16cm^{-3})$')
plt.ylabel('Slope')
plt.text(2, 1.2, 'Slope='+str(round(p_slope[0], 4))+'*density+'+str(round(p_slope[1], 4)), fontsize=15)
#plt.xscale('log')
#%%
xvar= np.array(Den[0:6])
yvar= FittingPara[0:6, 1]
fitfn= lambda p, xvar: p[0]*xvar+p[1]

errfunc = lambda p, xvar, yvar: fitfn(p, xvar)- yvar
p0= [1, 1]
p_offset, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))

#%%
NewX= np.linspace(np.amin(Den), np.amax(Den), 1000)
plt.plot(Den, FittingPara[:, 1], 'o')
plt.plot(NewX, fitfn(p_offset, NewX))
plt.xlabel('Density $(1e16cm^{-3})$')
plt.ylabel('Offset')
plt.text(2, -8, 'Offset='+str(round(p_offset[0], 4))+'*density+'+str(round(p_offset[1], 4)), fontsize=15)
plt.xscale('log')
#%%
xvar= np.array(Den)
yvar= FittingPara[:, 1]
fitfn= lambda p, xvar: p[0]+p[1]*xvar+p[2]*xvar**2+p[3]*xvar**3

errfunc = lambda p, xvar, yvar: fitfn(p, xvar)- yvar
p0= [1, 1, 1, 1]
p_offset, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))

#%%
NewX= np.linspace(np.amin(Den), np.amax(Den), 1000)
plt.plot(Den, FittingPara[:, 1], 'o')
plt.plot(NewX, fitfn(p_offset, NewX))
plt.xlabel('Density $(1e16cm^{-3})$')
plt.ylabel('Offset')
plt.xscale('log')
#plt.text(2, 3, 'Offset='+str(round(p_offset[0], 4))+str(round(p_offset[1], 4))+'*density+'\
#         +str(round(p_offset[2], 4))+'*density^2'+str(round(p_offset[3], 4))+'*density^3', fontsize=10)
#%%
plt.figure(4)
plt.plot(PWidth, AgW_n2[0, :], 'o-', label= '2e16')
plt.plot(PWidth, AgW_n5[0, :], 'o-', label= '5e16')
plt.plot(PWidth, AgW_n8[0, :], 'o-', label= '8e16')
plt.plot(PWidth, AgW_n10[0, :], 'o-', label= '10e16')
plt.plot(PWidth, AgW_n20[0, :], 'o-', label= '18e16')
plt.plot(PWidth, AgW_n40[0, :], 'o-', label= '40e16')
plt.plot(PWidth, AgW_n80[0, :], 'o-', label= '80e16')

plt.legend()
plt.xlabel('Plasma Half Width $(\mu m)$')
plt.ylabel('Afterglow Half Width $(\mu m)$')
#plt.title('n10 & n25 & n50')
#%%
plt.figure(5)
plt.plot(Den, AgW_pw30[0, :], 'o-', label= '30')
plt.plot(Den, AgW_pw35[0, :], 'o-', label= '35')
plt.plot(Den, AgW_pw40[0, :], 'o-', label= '40')
plt.plot(Den, AgW_pw45[0, :], 'o-', label= '45')
plt.plot(Den, AgW_pw50[0, :], 'o-', label= '50')
plt.xscale("log")
plt.legend()
plt.xlabel('Density $(1e16 cm^{-3})$')
plt.ylabel('Afterglow Half Width $(\mu m)$')
#plt.title('n10 & n25 & n50')

#%%
plt.plot(Den, FittingPara[:, 0], 'o')
#%%
plt.plot(Den, FittingPara[:, 1], 'o')
#%%
NewX= np.linspace(np.amin(PW), np.amax(PW), 1000)
plt.figure(4)
plt.plot(PWidth, AgW_n2[0, :], 'o', label= '2e16')
plt.plot(PWidth, AgW_n5[0, :], 'o', label= '5e16')
plt.plot(PWidth, AgW_n8[0, :], 'o', label= '8e16')
plt.plot(PWidth, AgW_n10[0, :], 'o', label= '10e16')
plt.plot(PWidth, AgW_n20[0, :], 'o', label= '18e16')
plt.plot(PWidth, AgW_n40[0, :], 'o', label= '40e16')
plt.plot(PWidth, AgW_n80[0, :], 'o', label= '80e16')
plt.plot(NewX, fitfn_All(FittingPara[0, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[1, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[2, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[3, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[4, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[5, :], NewX))
plt.plot(NewX, fitfn_All(FittingPara[6, :], NewX))


plt.legend()
plt.xlabel('Plasma Half Width $(\mu m)$')
plt.ylabel('Afterglow Half Width $(\mu m)$')
#plt.title('n10 & n25 & n50')
