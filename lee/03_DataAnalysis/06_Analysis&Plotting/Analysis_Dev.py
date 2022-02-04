#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 14:36:27 2021

@author: valentinalee
"""

import numpy as np
import pandas as pd
from scipy import optimize
import matplotlib.pyplot as plt
from scipy import ndimage
from sympy import symbols, Eq, solve
from scipy.interpolate import griddata
import os
#%%
Path= '~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan'
DataNumber= '2110290053'

#%%
def y_valley(x):
    y_valley= ((-a2/np.sqrt(a2**2+1))-(a1/np.sqrt(a1**2+1)))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))*x+ \
        ((-b2/np.sqrt(a2**2+1))-(b1/(np.sqrt(a1**2+1))))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))
    return y_valley

def ResultProcess(x_valley, data, BG):
    '''
    x_valley: int, the x coord. of the valley line
    data: 2D array, the 2D array of data (or signal) that is being processed
    --------
    return: 
    ResultList: 2D List, ([x_valley, Intensity of Peak1, Intensity of Peak2, Intensity of Valley, Distance of Peak1 and Peak2])
    '''
    global a_bis_v, x_grid, y_grid, a1, b1, a2, b2, x, y, griddata
    y_current= y_valley(x_valley)
    b_bis_v_current= y_current- a_bis_v*x_valley
    
    eq0= Eq(a_bis_v*x- y+ b_bis_v_current)
    eq1= Eq(a1*x- y+ b1)
    eq2= Eq(a2*x- y+ b2)
    x1y1= solve((eq0,eq1), (x, y))
    x2y2= solve((eq0,eq2), (x, y))
    print(x1y1)
    print(x2y2)
    print(x_valley, y_current)
    
    points= np.array((x_grid.flatten(), y_grid.flatten())).T
    values= data.flatten()
    BGvalues= BG.flatten()
    z0= girddata(points, values, (x_valley, y_current))
    z1= griddata(points, values, (x1y1[x], x1y1[y]))
    z2= griddata(points, values, (x2y2[x], x2y2[y]))
    z3= griddata(points, BGvalues, (x_valley, y_current))
    dist= np.sqrt((float(x1y1[x])-float(x2y2[x]))**2+ (float(x1y1[y])-float(x2y2[y]))**2)
    return [x_valley, float(z1), float(z2), float(z0), dist, float(z3)]


Path= os.path.expanduser(Path)
#WhateverData= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_Signal.npy')
#WhateverData= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_Data.npy')
BG= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_BG.npy')
Coords= np.load(Path+'/03_PeaksFoundByClickyTool'+'/'+DataNumber+'_PeaksInfo.npy')


X= np.linspace(0, 499, 500)
Y= np.linspace(0, 499, 500)
x_grid, y_grid= np.meshgrid(X, Y)

#%%

#Orgnize the data and zeros first 
Peak1X= Coords[:, 0][Coords[:, 0]!=0]
Peak1Y= Coords[:, 1][Coords[:, 1]!=0]
Peak2X= Coords[:, 2][Coords[:, 2]!=0]
Peak2Y= Coords[:, 3][Coords[:, 3]!=0]

#use the coords to do a linear fit and calculate the slope (for both)
xvar= Peak1X
yvar= Peak1Y
fitfn= lambda p, xvar: p[0]*xvar+p[1]

errfunc= lambda p, xvar, yvar: fitfn(p, xvar)- yvar
p0= [-1, 300]
p1, success= optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
a1= p1[0]
b1= p1[1]

xvar= Peak2X
yvar= Peak2Y
fitfn= lambda p, xvar: p[0]*xvar+p[1]

errfunc= lambda p, xvar, yvar: fitfn(p, xvar)- yvar
p0= [-1, 0]
p1, success= optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
a2= p1[0]
b2= p1[1]

a_bis= ((-a2/np.sqrt(a2**2+1))-(a1/np.sqrt(a1**2+1)))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))

a_bis_v= -1/a_bis
x, y= symbols('x y')
ResultList= [ResultProcess(x, WhateverData, BG) for x in range(250, 251)]

#%%

PV= [(ResultList[x][1]+ResultList[x][2])/2-ResultList[x][3] for x in range(0, np.shape(ResultList)[0])]
PV_avg= sum(PV)/ len(PV)
Separation= ResultList[int(np.shape(ResultList)[0]/2)][4]
#%%
#rotate the data
RotatedData= ndimage.rotate(Avg_Signal, np.rad2deg(TiltedAng), reshape=False)
#
#%%
plt.figure()
#plt.pcolormesh(RotatedData)
#plt.pcolormesh(Avg_data)
plt.pcolormesh(Avg_Signal)
#plt.plot(Peak1X, Peak1Y, '.')
#plt.plot(Peak2X, Peak2Y, '.')
#plt.plot(x_grid, a1*x_grid+b1)
#plt.plot(x_grid, a2*x_grid+b2)
#plt.plot(x_grid, y_mid)
#plt.plot(x, -(2*a2-a1)*x-(2*b2-b1))
#plt.plot(x, (-a2-a1)/2*x+(-b2-b1)/2)
plt.colorbar()
