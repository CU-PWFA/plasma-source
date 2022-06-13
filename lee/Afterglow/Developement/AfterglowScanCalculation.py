#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 17:03:39 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D

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
Den= [2, 5, 8, 10, 40, 80]
#Den= [2, 5, 8, 10, 20, 40, 80]
Pow= [5]
PWidth= [30, 35, 40, 45, 50]

#%%
for density in Den:
    print(density)
    for power in Pow:
        for pw in PWidth:
            
            PhotonMap= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
            IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
            r= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
            xvar= r[int(len(r)/2):len(r)]
            yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
            fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**p[2]/(2*p[1]**p[2])))

            errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
            p0= [1e11, 200, 2];
            p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
            PlasmaAfterglow[count][0]= density*1e16
            PlasmaAfterglow[count][1]= power
            PlasmaAfterglow[count][2]= pw*1e-6
            PlasmaAfterglow[count][3]= p1[0] #photon #
            PlasmaAfterglow[count][4]= p1[2] #power
            PlasmaAfterglow[count][5]= p1[1] #width
            PlasmaAfterglow= np.append(PlasmaAfterglow, [[0, 0, 0, 0, 0, 0]], axis=0)
            count= count+1
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
            
            PhotonMap= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
            IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
            r= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
            xvar= r[int(len(r)/2):len(r)]
            yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
            fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])

            errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
            p0= [1e11, 200, 2];
            p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
                        
            vars()['AgI_n'+str(density)][GPcount][PWcount]= p1[0] #photon #
            vars()['AgP_n'+str(density)][GPcount][PWcount]= p1[2] #power
            vars()['AgW_n'+str(density)][GPcount][PWcount]= p1[1] #width
#            vars()['AgFWHM_n'+str(density)][GPcount][PWcount]= fwhm #width
            
            PWcount= PWcount+1
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
            p0= [1e11, 200, 2];
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
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(Pow, den1, AgI_n1[:, 0], 'o-', label= 'PW= 30', color= 'r')
ax.plot(Pow, den1, AgI_n1[:, 1], 'o-', label= 'PW= 70', color= 'b')
ax.plot(Pow, den1, AgI_n1[:, 2], 'o-', label= 'PW= 110', color= 'g')
ax.plot(Pow, den1, AgI_n1[:, 3], 'o-', label= 'PW= 150', color= 'y')

ax.plot(Pow, den5, AgI_n5[:, 0], 'o-', color= 'r')
ax.plot(Pow, den5, AgI_n5[:, 1], 'o-', color= 'b')
ax.plot(Pow, den5, AgI_n5[:, 2], 'o-', color= 'g')
ax.plot(Pow, den5, AgI_n5[:, 3], 'o-', color= 'y')

ax.plot(Pow, den10, AgI_n10[:, 0], 'o-', color= 'r')
ax.plot(Pow, den10, AgI_n10[:, 1], 'o-', color= 'b')
ax.plot(Pow, den10, AgI_n10[:, 2], 'o-', color= 'g')
ax.plot(Pow, den10, AgI_n10[:, 3], 'o-', color= 'y')

ax.plot(Pow, den25, AgI_n25[:, 0], 'o-', color= 'r')
ax.plot(Pow, den25, AgI_n25[:, 1], 'o-', color= 'b')
ax.plot(Pow, den25, AgI_n25[:, 2], 'o-', color= 'g')
ax.plot(Pow, den25, AgI_n25[:, 3], 'o-', color= 'y')

ax.plot(Pow, den50, AgI_n50[:, 0], 'o-', color= 'r')
ax.plot(Pow, den50, AgI_n50[:, 1], 'o-', color= 'b')
ax.plot(Pow, den50, AgI_n50[:, 2], 'o-', color= 'g')
ax.plot(Pow, den50, AgI_n50[:, 3], 'o-', color= 'y')

ax.legend()
# make labels
ax.set_xlabel('Gaussian Power')
ax.set_ylabel('Density (1e16 $cm^{-3}$)')
ax.set_zlabel('Afterglow Photon #')

plt.show()

#%%
plt.figure(1)
plt.plot(Pow, AgI_n1[:, 0], 'o-', label= 'PW= 30', color= 'r')
plt.plot(Pow, AgI_n1[:, 1], 'o-', label= 'PW= 70', color= 'b')
plt.plot(Pow, AgI_n1[:, 2], 'o-', label= 'PW= 110', color= 'g')
plt.plot(Pow, AgI_n1[:, 3], 'o-', label= 'PW= 150', color= 'y')
plt.legend()
plt.xlabel('GP')
plt.ylabel('Afterglow Photon #')
plt.title('n1')

#plt.figure(2)
plt.plot(Pow, AgI_n5[:, 0], 'o--', color= 'r')
plt.plot(Pow, AgI_n5[:, 1], 'o--', color= 'b')
plt.plot(Pow, AgI_n5[:, 2], 'o--', color= 'g')
plt.plot(Pow, AgI_n5[:, 3], 'o--', color= 'y')
plt.legend()
#plt.xlabel('GP')
#plt.ylabel('Afterglow Photon #')
#plt.title('n5')

#plt.figure(3)
plt.plot(Pow, AgI_n10[:, 0], 'o:', color= 'r')
plt.plot(Pow, AgI_n10[:, 1], 'o:', color= 'b')
plt.plot(Pow, AgI_n10[:, 2], 'o:', color= 'g')
plt.plot(Pow, AgI_n10[:, 3], 'o:', color= 'y')
#plt.legend()
#plt.xlabel('GP')
#plt.ylabel('Afterglow Photon #')
#plt.title('n10')
#%%
#plt.figure(4)
plt.plot(Pow, AgI_n25[:, 0], 'o--', color= 'r')
plt.plot(Pow, AgI_n25[:, 1], 'o--', color= 'b')
plt.plot(Pow, AgI_n25[:, 2], 'o--', color= 'g')
plt.plot(Pow, AgI_n25[:, 3], 'o--', color= 'y')
#plt.legend()
#plt.xlabel('GP')
#plt.ylabel('Afterglow Photon #')
#plt.title('n25')

#plt.figure(5)
plt.plot(Pow, AgI_n50[:, 0], 'o-', label= 'PW= 30', color= 'r')
plt.plot(Pow, AgI_n50[:, 1], 'o-', label= 'PW= 70', color= 'b')
plt.plot(Pow, AgI_n50[:, 2], 'o-', label= 'PW= 110', color= 'g')
plt.plot(Pow, AgI_n50[:, 3], 'o-', label= 'PW= 150', color= 'y')
plt.legend()
plt.xlabel('GP')
plt.ylabel('Afterglow Photon #')
plt.title('n50')
#%%
plt.figure(1)
plt.plot(PlasmaAfterglow[:, 0]*1e-16, PlasmaAfterglow[:, 3], '.')
plt.xlabel('Plasma Density ($cm^{-3}$)')
plt.ylabel('Afterglow Photon #')

#%%
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x= Pow
y= PWidth
X, Y= np.meshgrid(x, y)
#Z1= AgI_n1
#Z5= AgI_n5
#Z10= AgI_n10
Z25= AgI_n25
Z50= AgI_n50
#ax.plot_wireframe(X, Y, Z1)
#ax.plot_wireframe(X, Y, Z5)
#ax.plot_wireframe(X, Y, Z10)
ax.plot_wireframe(X, Y, Z25, label= 'n25')
ax.plot_wireframe(X, Y, Z50, label= 'n50')
ax.set_xlabel('Gaussian Power')
ax.set_ylabel('Density (1e16 $cm^{-3}$)')
ax.set_zlabel('Afterglow Photon #')
ax.legend()
plt.show()


#%%
plt.figure(4)
#plt.plot(PWidth, AgI_n1[0, :], 'o-', label= 'GP= 3', color= 'r')
#plt.plot(PWidth, AgI_n1[1, :], 'o-', label= 'GP= 5', color= 'b')
#plt.plot(PWidth, AgI_n1[2, :], 'o-', label= 'GP= 7', color= 'g')
#plt.plot(PWidth, AgI_n1[3, :], 'o-', label= 'GP= 9', color= 'y')

#plt.plot(PWidth, AgI_n5[0, :], 'o-.', color= 'r')
#plt.plot(PWidth, AgI_n5[1, :], 'o-.', color= 'b')
#plt.plot(PWidth, AgI_n5[2, :], 'o-.', color= 'g')
#plt.plot(PWidth, AgI_n5[3, :], 'o-.', color= 'y')

plt.plot(PWidth, AgI_n10[0, :], 'o:', color= 'r')
plt.plot(PWidth, AgI_n10[1, :], 'o:', color= 'b')
plt.plot(PWidth, AgI_n10[2, :], 'o:', color= 'g')
plt.plot(PWidth, AgI_n10[3, :], 'o:', color= 'y')

plt.plot(PWidth, AgI_n25[0, :], 'o-', label= 'GP= 3', color= 'r')
plt.plot(PWidth, AgI_n25[1, :], 'o-', label= 'GP= 5', color= 'b')
plt.plot(PWidth, AgI_n25[2, :], 'o-', label= 'GP= 7', color= 'g')
plt.plot(PWidth, AgI_n25[3, :], 'o-', label= 'GP= 9', color= 'y')

plt.plot(PWidth, AgI_n50[0, :], 'o--', color= 'r')
plt.plot(PWidth, AgI_n50[1, :], 'o--', color= 'b')
plt.plot(PWidth, AgI_n50[2, :], 'o--', color= 'g')
plt.plot(PWidth, AgI_n50[3, :], 'o--', color= 'y')

plt.plot(PWidth, AgI_n10[0, :], 'o:', color= 'r', label= 'n10')
plt.plot(PWidth, AgI_n25[0, :], 'o-', color= 'r', label= 'n25')
plt.plot(PWidth, AgI_n50[0, :], 'o--', color= 'r', label= 'n50')

plt.legend()
plt.xlabel('PW')
plt.ylabel('Afterglow Photon #')
plt.title('n10 & n25 & n50')
#%%
plt.figure(4)
#plt.plot(Pow, AgI_n1[:, 0], 'o-', label= 'PW=30', color= 'r')
#plt.plot(Pow, AgI_n1[:, 1], 'o-', label= 'PW=70', color= 'b')
#plt.plot(Pow, AgI_n1[:, 2], 'o-', label= 'PW=110', color= 'g')
#plt.plot(Pow, AgI_n1[:, 3], 'o-', label= 'PW=150', color= 'y')

#plt.plot(Pow, AgI_n5[:, 0], 'o--', color= 'r')
#plt.plot(Pow, AgI_n5[:, 1], 'o--', color= 'b')
#plt.plot(Pow, AgI_n5[:, 2], 'o--', color= 'g')
#plt.plot(Pow, AgI_n5[:, 3], 'o--', color= 'y')

plt.plot(Pow, AgI_n10[:, 0], 'o:', color= 'r')
plt.plot(Pow, AgI_n10[:, 1], 'o:', color= 'b')
plt.plot(Pow, AgI_n10[:, 2], 'o:', color= 'g')
plt.plot(Pow, AgI_n10[:, 3], 'o:', color= 'y')

plt.plot(Pow, AgI_n25[:, 0], 'o-', label= 'PW=30', color= 'r')
plt.plot(Pow, AgI_n25[:, 1], 'o-', label= 'PW=70', color= 'b')
plt.plot(Pow, AgI_n25[:, 2], 'o-', label= 'PW=110', color= 'g')
plt.plot(Pow, AgI_n25[:, 3], 'o-', label= 'PW=150', color= 'y')

plt.plot(Pow, AgI_n50[:, 0], 'o-.', color= 'r')
plt.plot(Pow, AgI_n50[:, 1], 'o-.', color= 'b')
plt.plot(Pow, AgI_n50[:, 2], 'o-.', color= 'g')
plt.plot(Pow, AgI_n50[:, 3], 'o-.', color= 'y')

#plt.plot(Pow, AgI_n1[:, 0], 'o-', color= 'r', label= 'n1')
#plt.plot(Pow, AgI_n5[:, 0], 'o--', color= 'r', label= 'n5')
plt.plot(Pow, AgI_n10[:, 0], 'o:', color= 'r', label= 'n10')
plt.plot(Pow, AgI_n25[:, 0], 'o-', color= 'r', label= 'n25')
plt.plot(Pow, AgI_n50[:, 0], 'o-.', color= 'r', label= 'n50')

plt.legend()
plt.xlabel('GP')
plt.ylabel('Afterglow Photon #')
plt.title('nBig')
#%%
plt.figure(4)
plt.plot(PWidth, AgP_n1[0, :], 'o-', label= 'GP= 3', color= 'r')
plt.plot(PWidth, AgP_n1[1, :], 'o-', label= 'GP= 5', color= 'b')
plt.plot(PWidth, AgP_n1[2, :], 'o-', label= 'GP= 7', color= 'g')
plt.plot(PWidth, AgP_n1[3, :], 'o-', label= 'GP= 9', color= 'y')

plt.plot(PWidth, AgP_n5[0, :], 'o-.', color= 'r')
plt.plot(PWidth, AgP_n5[1, :], 'o-.', color= 'b')
plt.plot(PWidth, AgP_n5[2, :], 'o-.', color= 'g')
plt.plot(PWidth, AgP_n5[3, :], 'o-.', color= 'y')

plt.plot(PWidth, AgP_n10[0, :], 'o:', color= 'r')
plt.plot(PWidth, AgP_n10[1, :], 'o:', color= 'b')
plt.plot(PWidth, AgP_n10[2, :], 'o:', color= 'g')
plt.plot(PWidth, AgP_n10[3, :], 'o:', color= 'y')

#plt.plot(PWidth, AgP_n25[0, :], 'o-', label= 'GP= 3', color= 'r')
#plt.plot(PWidth, AgP_n25[1, :], 'o-', label= 'GP= 5', color= 'b')
#plt.plot(PWidth, AgP_n25[2, :], 'o-', label= 'GP= 7', color= 'g')
#plt.plot(PWidth, AgP_n25[3, :], 'o-', label= 'GP= 9', color= 'y')

#plt.plot(PWidth, AgP_n50[0, :], 'o--', color= 'r')
#plt.plot(PWidth, AgP_n50[1, :], 'o--', color= 'b')
#plt.plot(PWidth, AgP_n50[2, :], 'o--', color= 'g')
#plt.plot(PWidth, AgP_n50[3, :], 'o--', color= 'y')

plt.plot(PWidth, AgP_n1[0, :], 'o-', color= 'r', label= 'n1')
plt.plot(PWidth, AgP_n5[0, :], 'o-.', color= 'r', label= 'n5')
plt.plot(PWidth, AgP_n10[0, :], 'o:', color= 'r', label= 'n10')

plt.legend()
plt.xlabel('PW')
plt.ylabel('Afterglow Gaussian Power')
plt.title('n1 & n5 & n10')
#%%
plt.figure(4)
plt.plot(Pow, AgP_n1[:, 0], 'o-.', color= 'r')
plt.plot(Pow, AgP_n1[:, 1], 'o-.', color= 'b')
plt.plot(Pow, AgP_n1[:, 2], 'o-.', color= 'g')
plt.plot(Pow, AgP_n1[:, 3], 'o-.', color= 'y')

plt.plot(Pow, AgP_n5[:, 0], 'o--', color= 'r')
plt.plot(Pow, AgP_n5[:, 1], 'o--', color= 'b')
plt.plot(Pow, AgP_n5[:, 2], 'o--', color= 'g')
plt.plot(Pow, AgP_n5[:, 3], 'o--', color= 'y')

plt.plot(Pow, AgP_n10[:, 0], 'o:', color= 'r')
plt.plot(Pow, AgP_n10[:, 1], 'o:', color= 'b')
plt.plot(Pow, AgP_n10[:, 2], 'o:', color= 'g')
plt.plot(Pow, AgP_n10[:, 3], 'o:', color= 'y')

plt.plot(Pow, AgP_n25[:, 0], 'o-', label= 'PW=30', color= 'r')
plt.plot(Pow, AgP_n25[:, 1], 'o-', label= 'PW=70', color= 'b')
plt.plot(Pow, AgP_n25[:, 2], 'o-', label= 'PW=110', color= 'g')
plt.plot(Pow, AgP_n25[:, 3], 'o-', label= 'PW=150', color= 'y')

plt.plot(Pow, AgP_n50[:, 0], 'o-.', color= 'r')
plt.plot(Pow, AgP_n50[:, 1], 'o-.', color= 'b')
plt.plot(Pow, AgP_n50[:, 2], 'o-.', color= 'g')
plt.plot(Pow, AgP_n50[:, 3], 'o-.', color= 'y')

plt.plot(Pow, AgP_n1[:, 0], 'o-.', color= 'r', label= 'n1')
plt.plot(Pow, AgP_n5[:, 0], 'o--', color= 'r', label= 'n5')
plt.plot(Pow, AgP_n10[:, 0], 'o:', color= 'r', label= 'n10')
plt.plot(Pow, AgP_n25[:, 0], 'o-', color= 'r', label= 'n25')
plt.plot(Pow, AgP_n50[:, 0], 'o-.', color= 'r', label= 'n50')

plt.legend()
plt.xlabel('GP')
plt.ylabel('Afterglow Gaussian Power')
plt.title('nbig')

#%%
#%%
plt.figure(4)
plt.plot(PWidth, AgW_n2[0, :], 'o-', label= 'GP= 3', color= 'r')
#plt.plot(PWidth, AgW_n2[1, :], 'o-', label= 'GP= 5', color= 'b')
#plt.plot(PWidth, AgW_n2[2, :], 'o-', label= 'GP= 7', color= 'g')
#plt.plot(PWidth, AgW_n2[3, :], 'o-', label= 'GP= 9', color= 'y')

plt.plot(PWidth, AgW_n5[0, :], 'o-.', color= 'r')
#plt.plot(PWidth, AgW_n5[1, :], 'o-.', color= 'b')
#plt.plot(PWidth, AgW_n5[2, :], 'o-.', color= 'g')
#plt.plot(PWidth, AgW_n5[3, :], 'o-.', color= 'y')

plt.plot(PWidth, AgW_n8[0, :], 'o-.', color= 'r')
#plt.plot(PWidth, AgW_n8[1, :], 'o-.', color= 'b')
#plt.plot(PWidth, AgW_n8[2, :], 'o-.', color= 'g')
#plt.plot(PWidth, AgW_n8[3, :], 'o-.', color= 'y')

plt.plot(PWidth, AgW_n10[0, :], 'o:', color= 'r')
#plt.plot(PWidth, AgW_n10[1, :], 'o:', color= 'b')
#plt.plot(PWidth, AgW_n10[2, :], 'o:', color= 'g')
#plt.plot(PWidth, AgW_n10[3, :], 'o:', color= 'y')

#plt.plot(PWidth, AgW_n25[0, :], 'o-', label= 'GP= 3', color= 'r')
#plt.plot(PWidth, AgW_n25[1, :], 'o-', label= 'GP= 5', color= 'b')
#plt.plot(PWidth, AgW_n25[2, :], 'o-', label= 'GP= 7', color= 'g')
#plt.plot(PWidth, AgW_n25[3, :], 'o-', label= 'GP= 9', color= 'y')

#plt.plot(PWidth, AgW_n50[0, :], 'o--', color= 'r')
#plt.plot(PWidth, AgW_n50[1, :], 'o--', color= 'b')
#plt.plot(PWidth, AgW_n50[2, :], 'o--', color= 'g')
#plt.plot(PWidth, AgW_n50[3, :], 'o--', color= 'y')

#plt.plot(PWidth, AgW_n10[0, :], 'o:', color= 'r', label= 'n10')
#plt.plot(PWidth, AgW_n25[0, :], 'o-', color= 'r', label= 'n25')
#plt.plot(PWidth, AgW_n50[0, :], 'o--', color= 'r', label= 'n50')

plt.legend()
plt.xlabel('PW')
plt.ylabel('Afterglow Width')
plt.title('n10 & n25 & n50')
#%%
plt.figure(4)
plt.plot(Pow, AgW_n1[:, 0], 'o-.', color= 'r')
plt.plot(Pow, AgW_n1[:, 1], 'o-.', color= 'b')
plt.plot(Pow, AgW_n1[:, 2], 'o-.', color= 'g')
plt.plot(Pow, AgW_n1[:, 3], 'o-.', color= 'y')

plt.plot(Pow, AgW_n5[:, 0], 'o--', color= 'r')
plt.plot(Pow, AgW_n5[:, 1], 'o--', color= 'b')
plt.plot(Pow, AgW_n5[:, 2], 'o--', color= 'g')
plt.plot(Pow, AgW_n5[:, 3], 'o--', color= 'y')

plt.plot(Pow, AgW_n10[:, 0], 'o:', color= 'r')
plt.plot(Pow, AgW_n10[:, 1], 'o:', color= 'b')
plt.plot(Pow, AgW_n10[:, 2], 'o:', color= 'g')
plt.plot(Pow, AgW_n10[:, 3], 'o:', color= 'y')

plt.plot(Pow, AgW_n25[:, 0], 'o-', label= 'PW=30', color= 'r')
plt.plot(Pow, AgW_n25[:, 1], 'o-', label= 'PW=70', color= 'b')
plt.plot(Pow, AgW_n25[:, 2], 'o-', label= 'PW=110', color= 'g')
plt.plot(Pow, AgW_n25[:, 3], 'o-', label= 'PW=150', color= 'y')

plt.plot(Pow, AgW_n50[:, 0], 'o-.', color= 'r')
plt.plot(Pow, AgW_n50[:, 1], 'o-.', color= 'b')
plt.plot(Pow, AgW_n50[:, 2], 'o-.', color= 'g')
plt.plot(Pow, AgW_n50[:, 3], 'o-.', color= 'y')

plt.plot(Pow, AgW_n1[:, 0], 'o-.', color= 'r', label= 'n1')
plt.plot(Pow, AgW_n5[:, 0], 'o--', color= 'r', label= 'n5')
plt.plot(Pow, AgW_n10[:, 0], 'o:', color= 'r', label= 'n10')
plt.plot(Pow, AgW_n25[:, 0], 'o-', color= 'r', label= 'n25')
plt.plot(Pow, AgW_n50[:, 0], 'o-.', color= 'r', label= 'n50')

plt.legend()
plt.xlabel('GP')
plt.ylabel('Afterglow Width')
plt.title('nall')

#%%
#%%
plt.figure(4)
plt.plot(PWidth, AgW_n2[0, :], 'o-', label= '2e16')
plt.plot(PWidth, AgW_n5[0, :], 'o-', label= '5e16')
plt.plot(PWidth, AgW_n8[0, :], 'o-', label= '8e16')
plt.plot(PWidth, AgW_n10[0, :], 'o-', label= '10e16')
#plt.plot(PWidth, AgW_n20[0, :], 'o-', label= '18e16')
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