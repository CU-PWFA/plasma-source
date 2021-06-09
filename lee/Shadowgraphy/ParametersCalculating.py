#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 06:37:16 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
Physical_Y= 60000
Newgrid= 10000
newY= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)

PowerL= 9
PlasmaWidthL= 11
PWLine= np.linspace(20, 20+(10*PlasmaWidthL), (PlasmaWidthL+1))
GPLine= np.linspace(1, PowerL+1, PowerL+1)

PeakSpacingM= np.zeros((PowerL, PlasmaWidthL))
PeakSpacing2M= np.zeros((PowerL, PlasmaWidthL))
PeakHRatioM= np.zeros((PowerL, PlasmaWidthL))


#for power in range (9, 10):    
for power in range (1, int(PowerL+1)):    
#    for PlasmaWidthS in range (30, 40, 10):
    for PlasmaWidthS in range (20, 130, 10):
        LONew= np.load('PW'+str(PlasmaWidthS)+'_GP'+str(power)+'.npy')
        
        CIdx= int(LONew.shape[0]/2)
        MaxIdx= int(LONew.shape[0])
        for n in range (CIdx, int(MaxIdx-250)):
            if LONew[n+10]-LONew[n] <0:
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
        for n in range (int(CIdx-100), 250, -1):
            if LONew[n+10]-LONew[n] <0:
                if LONew[n+50]-LONew[n] <0:
                    if LONew[n+100]-LONew[n] <0:            
                        if LONew[n+150]-LONew[n] <0:            
                            if LONew[n+250]-LONew[n] <0:            
                                if LONew[n-1]-LONew[n] <0:            
                                    if LONew[n-50]-LONew[n] <0:            
                                        if LONew[n-100]-LONew[n] <0:            
                                            if LONew[n-150]-LONew[n] <0:            
                                                if LONew[n-250]-LONew[n] <0:            
                                                    print(-n)
                                                    PeakIdxL= n
                                                    break

        for n in range (int(PeakIdxH+100), int(MaxIdx-250)):
            if LONew[n+10]-LONew[n] <0:
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
                                                    PeakIdxH2= n
                                                    break
        for n in range (int(PeakIdxL-100), 250, -1):
            if LONew[n+10]-LONew[n] <0:
                if LONew[n+50]-LONew[n] <0:
                    if LONew[n+100]-LONew[n] <0:            
                        if LONew[n+150]-LONew[n] <0:            
                            if LONew[n+250]-LONew[n] <0:            
                                if LONew[n-1]-LONew[n] <0:            
                                    if LONew[n-50]-LONew[n] <0:            
                                        if LONew[n-100]-LONew[n] <0:            
                                            if LONew[n-150]-LONew[n] <0:            
                                                if LONew[n-250]-LONew[n] <0:            
                                                    print(-n)
                                                    PeakIdxL2= n
                                                    break


        PatternL= PeakIdxH-PeakIdxL
        PatternL2= PeakIdxH2-PeakIdxL2
        pixel= Physical_Y/Newgrid
        PeakSpacing= PatternL*pixel
        PeakSpacing2= PatternL2*pixel
        PeakHRatio= LONew[PeakIdxH]/LONew[CIdx]
        PhysicalGrid= np.linspace(-Physical_Y/2, Physical_Y/2, Newgrid)
        print(PeakSpacing)
        
        PeakSpacingM[int(power-1), int((PlasmaWidthS-20)/10)]= PeakSpacing
        PeakSpacing2M[int(power-1), int((PlasmaWidthS-20)/10)]= PeakSpacing2
        PeakHRatioM[int(power-1), int((PlasmaWidthS-20)/10)]= PeakHRatio
#%%
#plt.pcolormesh(PeakSpacingM)
plt.figure(2)
#plt.title('Plasma Refraction Experiment \n Peak-to-trough Height Ratio 2D Parameters Scan')
plt.title('Plasma Refraction Experiment \n Second PeakSpacing 2D Parameters Scan')
#plt.pcolormesh(PWLine, GPLine, PeakHRatioM, cmap= 'jet')
plt.pcolormesh(PWLine, GPLine, PeakSpacing2M, cmap= 'jet')
plt.xlabel('w ($\mu$m)')
plt.ylabel('p')
plt.colorbar(label= '$\mu$m')
#%%
PWLineOL= np.linspace(20, 20+(10*(PlasmaWidthL-1)), (PlasmaWidthL))
GPLineOL= np.linspace(1, PowerL, PowerL)

#%%
plt.plot(PWLineOL, PeakSpacingM[0, :], 'o-', label= 'p= 1')
plt.plot(PWLineOL, PeakSpacingM[1, :], 'o-', label= 'p= 2')
plt.plot(PWLineOL, PeakSpacingM[2, :], 'o-', label= 'p= 3')
plt.plot(PWLineOL, PeakSpacingM[3, :], 'o-', label= 'p= 4')
plt.plot(PWLineOL, PeakSpacingM[4, :], 'o-', label= 'p= 5')
plt.plot(PWLineOL, PeakSpacingM[5, :], 'o-', label= 'p= 6')
plt.plot(PWLineOL, PeakSpacingM[6, :], 'o-', label= 'p= 7')
plt.plot(PWLineOL, PeakSpacingM[7, :], 'o-', label= 'p= 8')
plt.plot(PWLineOL, PeakSpacingM[8, :], 'o-', label= 'p= 9')
plt.legend()
plt.xlabel('Plasma w ($\mu$m)')
plt.ylabel('First PeakSpacing ($\mu$m)')

#%%
plt.plot(PWLineOL, PeakSpacing2M[0, :], 'o-', label= 'p= 1')
plt.plot(PWLineOL, PeakSpacing2M[1, :], 'o-', label= 'p= 2')
plt.plot(PWLineOL, PeakSpacing2M[2, :], 'o-', label= 'p= 3')
plt.plot(PWLineOL, PeakSpacing2M[3, :], 'o-', label= 'p= 4')
plt.plot(PWLineOL, PeakSpacing2M[4, :], 'o-', label= 'p= 5')
plt.plot(PWLineOL, PeakSpacing2M[5, :], 'o-', label= 'p= 6')
plt.plot(PWLineOL, PeakSpacing2M[6, :], 'o-', label= 'p= 7')
plt.plot(PWLineOL, PeakSpacing2M[7, :], 'o-', label= 'p= 8')
plt.plot(PWLineOL, PeakSpacing2M[8, :], 'o-', label= 'p= 9')
plt.legend()
plt.xlabel('Plasma w ($\mu$m)')
plt.ylabel('Second PeakSpacing ($\mu$m)')

#%%
#plt.plot(PWLineOL, PeakHRatioM[0][:], 'o-', label= 'p= 1')
plt.plot(PWLineOL, PeakHRatioM[1, :], 'o-', label= 'p= 2')
plt.plot(PWLineOL, PeakHRatioM[2, :], 'o-', label= 'p= 3')
plt.plot(PWLineOL, PeakHRatioM[3, :], 'o-', label= 'p= 4')
plt.plot(PWLineOL, PeakHRatioM[4, :], 'o-', label= 'p= 5')
plt.plot(PWLineOL, PeakHRatioM[5, :], 'o-', label= 'p= 6')
plt.plot(PWLineOL, PeakHRatioM[6, :], 'o-', label= 'p= 7')
plt.plot(PWLineOL, PeakHRatioM[7, :], 'o-', label= 'p= 8')
plt.plot(PWLineOL, PeakHRatioM[8, :], 'o-', label= 'p= 9')
plt.legend()
plt.xlabel('Plasma w ($\mu$m)')
plt.ylabel('Peak-to-trough Height Ratio')

#%%
plt.plot(GPLineOL, PeakSpacingM[:, 0], 'o-', label= 'w= 20')
plt.plot(GPLineOL, PeakSpacingM[:, 1], 'o-', label= 'w= 30')
plt.plot(GPLineOL, PeakSpacingM[:, 2], 'o-', label= 'w= 40')
plt.plot(GPLineOL, PeakSpacingM[:, 3], 'o-', label= 'w= 50')
plt.plot(GPLineOL, PeakSpacingM[:, 4], 'o-', label= 'w= 60')
plt.plot(GPLineOL, PeakSpacingM[:, 5], 'o-', label= 'w= 70')
plt.plot(GPLineOL, PeakSpacingM[:, 6], 'o-', label= 'w= 80')
plt.plot(GPLineOL, PeakSpacingM[:, 7], 'o-', label= 'w= 90')
plt.plot(GPLineOL, PeakSpacingM[:, 8], 'o-', label= 'w= 100')
plt.plot(GPLineOL, PeakSpacingM[:, 9], 'o-', label= 'w= 110')
plt.plot(GPLineOL, PeakSpacingM[:, 10], 'o-', label= 'w= 120')
plt.legend()
plt.xlabel('Gaussian p')
plt.ylabel('First PeakSpacing ($\mu$m)')

#%%
plt.plot(GPLineOL, PeakSpacing2M[:, 0], 'o-', label= 'w= 20')
plt.plot(GPLineOL, PeakSpacing2M[:, 1], 'o-', label= 'w= 30')
plt.plot(GPLineOL, PeakSpacing2M[:, 2], 'o-', label= 'w= 40')
plt.plot(GPLineOL, PeakSpacing2M[:, 3], 'o-', label= 'w= 50')
plt.plot(GPLineOL, PeakSpacing2M[:, 4], 'o-', label= 'w= 60')
plt.plot(GPLineOL, PeakSpacing2M[:, 5], 'o-', label= 'w= 70')
plt.plot(GPLineOL, PeakSpacing2M[:, 6], 'o-', label= 'w= 80')
plt.plot(GPLineOL, PeakSpacing2M[:, 7], 'o-', label= 'w= 90')
plt.plot(GPLineOL, PeakSpacing2M[:, 8], 'o-', label= 'w= 100')
plt.plot(GPLineOL, PeakSpacing2M[:, 9], 'o-', label= 'w= 110')
plt.plot(GPLineOL, PeakSpacing2M[:, 10], 'o-', label= 'w= 120')
plt.legend()
plt.xlabel('Gaussian p')
plt.ylabel('Second PeakSpacing ($\mu$m)')

#%%
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 0], 'o-', label= 'w= 20')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 1], 'o-', label= 'w= 30')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 2], 'o-', label= 'w= 40')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 3], 'o-', label= 'w= 50')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 4], 'o-', label= 'w= 60')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 5], 'o-', label= 'w= 70')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 6], 'o-', label= 'w= 80')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 7], 'o-', label= 'w= 90')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 8], 'o-', label= 'w= 100')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 9], 'o-', label= 'w= 110')
plt.plot(GPLineOL[1:9], PeakHRatioM[1:9, 10], 'o-', label= 'w= 120')
plt.legend()
plt.xlabel('Gaussian p')
plt.ylabel('Peak-to-trough Height Ratio')

#%%
plt.title('Normalized Nominal Case \n Plasma w=30 n=1e17')
plt.plot(GPLineOL, PeakHRatioM[:, 1]/(np.min(PeakHRatioM[:, 1])), 'o-', label= 'Peak-to-trough Height Ratio')
plt.plot(GPLineOL, PeakSpacingM[:, 1]/np.min(PeakSpacingM[:, 1]), 'o-', label= 'Frist PeakSpacing')
plt.plot(GPLineOL, PeakSpacing2M[:, 1]/(np.min(PeakSpacing2M[:, 1])), 'o-', label= 'Second PeakSpacing')
plt.legend()
plt.xlabel('Gaussian p')
plt.ylabel('Sensitivities (Normalized Parameters)')

#%%
plt.title('Normalized Nominal Case \n Gaussian p=5 n=1e17')
plt.plot(PWLineOL, PeakHRatioM[4, :]/(np.min(PeakHRatioM[4, :])), 'o-', label= 'Peak-to-trough Height Ratio')
plt.plot(PWLineOL, PeakSpacingM[4, :]/(np.min(PeakSpacingM[4, :])), 'o-', label= 'Frist PeakSpacing')
plt.plot(PWLineOL, PeakSpacing2M[4, :]/np.min(PeakSpacing2M[4, :]), 'o-', label= 'Second PeakSpacing')
plt.legend()
plt.xlabel('Plasma w ($\mu$m)')
plt.ylabel('Sensitivities (Normalized Parameters)')

