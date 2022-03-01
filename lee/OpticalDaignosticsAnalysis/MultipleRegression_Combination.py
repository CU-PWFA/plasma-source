#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 17:44:35 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import find_peaks
from sklearn.linear_model import LinearRegression
from itertools import combinations
#%%
Den= [2, 5, 8, 10, 20, 40, 50, 80]
Pow= [5]
PWidth= [30, 35, 40, 45, 50]

#%%
PixelSize= 10

count= 0

XMatrix= []
BasePath_AG= '/media/valentinalee/TRANSCEND/AfterglowFinalPowerFix/ScanResults/AfterglowScan/'
BasePath_Sh= '/media/valentinalee/TRANSCEND/ShadowgraphyResults/FFResults2/'
#%%
for density in Den:
    for power in Pow:
        for pw in PWidth:
            try:
                
                PhotonMap= np.load(BasePath_AG+'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
                IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
                r= np.load(BasePath_AG+'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
                xvar= r[int(len(r)/2):len(r)]
                yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
                fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])
    
                errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
                p0= [1e11, 20, 2];
                p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))

                dataE= np.load(BasePath_Sh+'PW'+str(pw)+'_GP'+str(power)+'_n'+str(density)+'.npy')
                MidLineOutI= abs(dataE[:, int(dataE.shape[1]/2)])**2

                x= MidLineOutI
                peaks, _ = find_peaks(x)
                peak1idx, peak2idx= peaks[np.argsort(x[peaks])[-1]], peaks[np.argsort(x[peaks])[-2]]
                peak1, peak2= x[peak1idx], x[peak2idx]
                
                if (np.argsort(x[peaks])[-1])-(np.argsort(x[peaks])[-2])==2:
                   
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx, valley2idx= valleys[np.argsort(x[valleys])[0:2]]
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                else:
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx= valleys[np.argsort(x[valleys])[0]]
                    valley2idx= valley1idx
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                    
                if (np.mean(np.array([peak1, peak2])))-(np.mean(np.array([valley1, valley2])))>0.5:
                    
                    plt.plot(x)
                    plt.plot(np.array([peak1idx, peak2idx]), x[np.array([peak1idx, peak2idx])], 'x')
                    plt.show()     
                    plt.plot(np.array([valley1idx, valley2idx]), x[np.array([valley1idx, valley2idx])], 'x')
    
                    FirstPeakWidth= abs(peak1idx-peak2idx)*PixelSize
                    PTValue= 100*(np.mean(np.array([peak1, peak2]))-np.mean(np.array([valley1, valley2])))/ \
                                ((np.mean(np.array([peak1, peak2]))+np.mean(np.array([valley1, valley2])))/2)
                    XMatrix.append([density, pw, p1[1], p1[0], FirstPeakWidth, PTValue])
                    
                    count= count+1

            except:
                print(density, power, pw)
                print('Oh no!')
                pass
                
XMatrix= np.array(XMatrix)            
#%%first order terms and second order temrs
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
NewXMatrix[:, 0]= XMatrix[:, 2]/np.amax(XMatrix[:, 2])
NewXMatrix[:, 1]= XMatrix[:, 3]/np.amax(XMatrix[:, 3])
NewXMatrix[:, 2]= XMatrix[:, 4]/np.amax(XMatrix[:, 4])
NewXMatrix[:, 3]= XMatrix[:, 5]/np.amax(XMatrix[:, 5])
NewXMatrix[:, 4]= XMatrix[:, 2]*XMatrix[:, 3]/np.amax(XMatrix[:, 2]*XMatrix[:, 3])
NewXMatrix[:, 5]= XMatrix[:, 2]*XMatrix[:, 4]/np.amax(XMatrix[:, 2]*XMatrix[:, 4])
NewXMatrix[:, 6]= XMatrix[:, 2]*XMatrix[:, 5]/np.amax(XMatrix[:, 2]*XMatrix[:, 5])
NewXMatrix[:, 7]= XMatrix[:, 3]*XMatrix[:, 4]/np.amax(XMatrix[:, 3]*XMatrix[:, 4])
NewXMatrix[:, 8]= XMatrix[:, 3]*XMatrix[:, 5]/np.amax(XMatrix[:, 3]*XMatrix[:, 5])
NewXMatrix[:, 9]= XMatrix[:, 4]*XMatrix[:, 5]/np.amax(XMatrix[:, 4]*XMatrix[:, 5])
NewXMatrix[:, 10]= XMatrix[:, 2]**2/np.amax(XMatrix[:, 2]**2)
NewXMatrix[:, 11]= XMatrix[:, 3]**2/np.amax(XMatrix[:, 3]**2)
NewXMatrix[:, 12]= XMatrix[:, 4]**2/np.amax(XMatrix[:, 4]**2)
NewXMatrix[:, 13]= XMatrix[:, 5]**2/np.amax(XMatrix[:, 5]**2)

#%%for n
comb= combinations([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 3)
XM= np.zeros((XMatrix.shape[0],  3))
Score_for_n= []
for c in list(comb):
    XM[:, 0]= NewXMatrix[:, c[0]]
    XM[:, 1]= NewXMatrix[:, c[1]]
    XM[:, 2]= NewXMatrix[:, c[2]]
#    XM[:, 3]= NewXMatrix[:, c[3]]
#    XM4[:, 4]= NewXMatrix[:, c[4]]    
    n_reg= LinearRegression().fit(XM, XMatrix[:, 0])
    n_coef= n_reg.coef_
    n_intercept= n_reg.intercept_
    n_score= n_reg.score(XM, XMatrix[:, 0])
    Score_for_n.append([n_score, c])
    n_predict= n_intercept+ n_coef[0]*XM[:, 0]+ n_coef[1]*XM[:, 1]+ n_coef[2]*XM[:, 2]#+ \
#           n_coef[3]*XM[:, 3]#+ n_coef[4]*NewXMatrix[:, 4]+ n_coef[5]*NewXMatrix[:, 5]+ \
#           n_coef[6]*NewXMatrix[:, 6]
    errorM= np.zeros((n_predict.shape[0], 4))
    errorM[:, 0]= XMatrix[:, 0]
    errorM[:, 1]= XMatrix[:, 1]
    errorM[:, 2]= abs(XMatrix[:, 0]-n_predict)/XMatrix[:, 0]
#    errorM[:, 3]= abs(XMatrix[:, 1]-w_predict)/XMatrix[:, 1]
    
    if np.mean(sorted(errorM[:, 2])[0:-1]) < 0.1:
        if np.amax(sorted(errorM[:, 2])[0:-1]) < 0.15:
            plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 2], s= 500)
            plt.colorbar()
            plt.xscale('log')
            plt.title('Density Error'+str(c))
            plt.xlabel('True Density (1e16 cm$^{-3}$)')
            plt.ylabel('True Plasma Width ($\mu$m)')
            plt.savefig('n_'+str(c))
            plt.close()

#%%f w
comb= combinations([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], 2)
XM= np.zeros((XMatrix.shape[0],  2))
Score_for_w= []
for c in list(comb):
    XM[:, 0]= NewXMatrix[:, c[0]]
    XM[:, 1]= NewXMatrix[:, c[1]]
#    XM[:, 2]= NewXMatrix[:, c[2]]
#    XM[:, 3]= NewXMatrix[:, c[3]]
#    XM4[:, 4]= NewXMatrix[:, c[4]]    
    w_reg= LinearRegression().fit(XM, XMatrix[:, 1])
    w_coef= w_reg.coef_
    w_intercept= w_reg.intercept_
    w_score= w_reg.score(XM, XMatrix[:, 1])
    Score_for_w.append([w_score, c])
    w_predict= w_intercept+ w_coef[0]*XM[:, 0]+ w_coef[1]*XM[:, 1]#+ w_coef[2]*XM[:, 2]#+ \
#           w_coef[3]*XM[:, 3]#+ w_coef[4]*NewXMatrix[:, 4]+ w_coef[5]*NewXMatrix[:, 5]+ \
#           n_coef[6]*NewXMatrix[:, 6]
    errorM= np.zeros((n_predict.shape[0], 4))
    errorM[:, 0]= XMatrix[:, 0]
    errorM[:, 1]= XMatrix[:, 1]
    errorM[:, 3]= abs(XMatrix[:, 1]-w_predict)/XMatrix[:, 1]
    
#    if np.mean(sorted(errorM[:, 3])[0:-1]) < 0.05:
    if np.mean(errorM[:, 3]) < 0.03:
        if np.amax(errorM[:, 3]) < 0.08:
#        if np.amax(sorted(errorM[:, 3])[0:-1]) < 0.08:
            plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 3], s= 500)
            plt.colorbar()
            plt.xscale('log')
            plt.title('Width Error'+str(c))
            plt.xlabel('True Density (1e16 cm$^{-3}$)')
            plt.ylabel('True Plasma Width ($\mu$m)')
            plt.savefig('w_'+str(c))
            plt.close()