#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 09:02:51 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
PW30= np.load('W30.npy')
PW40= np.load('W40.npy')
PW50= np.load('W50.npy')
PW60= np.load('W60.npy')
PW70= np.load('W70.npy')
PW80= np.load('W80.npy')
PW90= np.load('W90.npy')
PW100= np.load('W100.npy')
PW120= np.load('W120.npy')
PW140= np.load('W140.npy')
PW160= np.load('W160.npy')
PW180= np.load('W180.npy')
PW200= np.load('W200.npy')

#%%
LONew= PW200
Physical_Y= 60000
Newgrid= 10000
newY= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)

#%%
CIdx= int(LONew.shape[0]/2)
MaxIdx= int(LONew.shape[0])
for n in range (CIdx, int(MaxIdx-250)):
    if LONew[n+1]-LONew[n] <0:
        if LONew[n+50]-LONew[n] <0:
            if LONew[n+100]-LONew[n] <0:            
                if LONew[n+150]-LONew[n] <0:            
                    if LONew[n+250]-LONew[n] <0:            
                        if LONew[n-1]-LONew[n] <0:            
                            if LONew[n-50]-LONew[n] <0:            
                                if LONew[n-100]-LONew[n] <0:            
                                    if LONew[n-150]-LONew[n] <0:            
                                        if LONew[n-250]-LONew[n] <0:            
                                            print(n)
                                            PeakIdxH= n
                                            break
for n in range (CIdx, 250, -1):
    if LONew[n+1]-LONew[n] <0:
        if LONew[n+50]-LONew[n] <0:
            if LONew[n+100]-LONew[n] <0:            
                if LONew[n+150]-LONew[n] <0:            
                    if LONew[n+250]-LONew[n] <0:            
                        if LONew[n-1]-LONew[n] <0:            
                            if LONew[n-50]-LONew[n] <0:            
                                if LONew[n-100]-LONew[n] <0:            
                                    if LONew[n-150]-LONew[n] <0:            
                                        if LONew[n-250]-LONew[n] <0:            
                                            print(n)
                                            PeakIdxL= n
                                            break
#%%
PatternL= PeakIdxH-PeakIdxL
pixel= Physical_Y/Newgrid
PeakSpacing= PatternL*pixel
PeakHRatio= LONew[PeakIdxH]/LONew[CIdx]
PeakHeigh= LONew[PeakIdxH]-LONew[CIdx]

#%%
#W30= PeakSpacing
#HR30= PeakHRatio
#H30= PeakHeigh
#W40= PeakSpacing
#HR40= PeakHRatio
#H40= PeakHeigh
#W50= PeakSpacing
#HR50= PeakHRatio
#H50= PeakHeigh
#W60= PeakSpacing
#HR60= PeakHRatio
#H60= PeakHeigh
#W70= PeakSpacing
#HR70= PeakHRatio
#H70= PeakHeigh
#W80= PeakSpacing
#HR80= PeakHRatio
#H80= PeakHeigh
#W90= PeakSpacing
#HR90= PeakHRatio
#H90= PeakHeigh
#W100= PeakSpacing
#HR100= PeakHRatio
#H100= PeakHeigh
#W120= PeakSpacing
#HR120= PeakHRatio
#H120= PeakHeigh
#W140= PeakSpacing
#HR140= PeakHRatio
#H140= PeakHeigh
#W160= PeakSpacing
#HR160= PeakHRatio
#H160= PeakHeigh
#W180= PeakSpacing
#HR180= PeakHRatio
#H180= PeakHeigh
W200= PeakSpacing
HR200= PeakHRatio
H200= PeakHeigh
#%%
WHMatrix= np.zeros((13, 3))
WHMatrix[0][0]= W30
WHMatrix[0][1]= HR30
WHMatrix[0][2]= H30
WHMatrix[1][0]= W40
WHMatrix[1][1]= HR40
WHMatrix[1][2]= H40
WHMatrix[2][0]= W50
WHMatrix[2][1]= HR50
WHMatrix[2][2]= H50
WHMatrix[3][0]= W60
WHMatrix[3][1]= HR60
WHMatrix[3][2]= H60
WHMatrix[4][0]= W70
WHMatrix[4][1]= HR70
WHMatrix[4][2]= H70
WHMatrix[5][0]= W80
WHMatrix[5][1]= HR80
WHMatrix[5][2]= H80
WHMatrix[6][0]= W90
WHMatrix[6][1]= HR90
WHMatrix[6][2]= H90
WHMatrix[7][0]= W100
WHMatrix[7][1]= HR100
WHMatrix[7][2]= H100
WHMatrix[8][0]= W120
WHMatrix[8][1]= HR120
WHMatrix[8][2]= H120
WHMatrix[9][0]= W140
WHMatrix[9][1]= HR140
WHMatrix[9][2]= H140
WHMatrix[10][0]= W160
WHMatrix[10][1]= HR160
WHMatrix[10][2]= H160
WHMatrix[11][0]= W180
WHMatrix[11][1]= HR180
WHMatrix[11][2]= H180
WHMatrix[12][0]= W200
WHMatrix[12][1]= HR200
WHMatrix[12][2]= H200

#%%
plt.plot(newY, PW40, label= 'PW40')
plt.plot(newY, PW60, label= 'PW60')
plt.plot(newY, PW80, label= 'PW80')
plt.plot(newY, PW100, label= 'PW100')
plt.plot(newY, PW120, label= 'PW120')
plt.plot(newY, PW140, label= 'PW140')
plt.plot(newY, PW160, label= 'PW160')
plt.plot(newY, PW180, label= 'PW180')
plt.plot(newY, PW200, label= 'PW200')
plt.legend()
plt.xlabel('Vertical Lineout ($\mu$m)')
plt.ylabel('Intensity (W/cm^2)')

#%%
PlasmaW=np.array([30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]) 
plt.plot(PlasmaW, WHMatrix[:, 0], 'o')
plt.xlabel('Plasma Width $\mu$m')
#plt.ylabel('Peak-Valley Ratio')
#plt.ylabel('Abs(Peak-Valley) W/cm^2')
plt.ylabel('FirstPeak Spacing $\mu$m')
