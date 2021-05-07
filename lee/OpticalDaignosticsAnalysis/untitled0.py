#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:13:21 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import find_peaks
from scipy import interpolate
#%%
def xyFINDz (Xaxis, Yaxis, CorrZ, GivenX, GivenY):
    '''
    Args:
        Xaxis: 1D array, the x axis of the 2D array
        Yaxis: 1D array, the y axis of the 2D array
        CorrZ: 2D array, z values
        GivenX: float, the given x where you want to find its z(x, y)
        GivenY: float, the given y where you want to find its z(x, y)
    
    Returns:
        ZValue: float, z(GivenX, GivenY)
    
    '''
    idxX= (np.abs(Xaxis-GivenX)).argmin()
    idxY= (np.abs(Yaxis-GivenY)).argmin()
    return(CorrZ[idxY, idxX])
    
def FindShadowgraphyMrm (MidLineOutI, PixelSize):

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
        
        FirstPeakWidth= abs(peak1idx-peak2idx)*PixelSize
        PTValue= 100*(np.mean(np.array([peak1, peak2]))-np.mean(np.array([valley1, valley2])))/ \
                    ((np.mean(np.array([peak1, peak2]))+np.mean(np.array([valley1, valley2])))/2)
        return FirstPeakWidth, PTValue

def AllPointAtThisZ (z, Zpoint, error= 0.1):
    '''
    Args:
        z: 2D array, z value, z(x, y)
        Zpoint: float, given point Zpoint, this fn will return points that are close to Zpoint
        error= 0.01: acceptance error, the function will return all the points where abs(z-Zpoint)<error
    
    Returns:
        x_acpt: 1D array, x value of the accepted points
        y_acpt: 1D array, y value of the accepted points
        z_acpt: 1D array, z value of the accepted points
    
    '''
    x_acpt= []
    y_acpt= []
    z_acpt= []
    for xidx in range (z.shape[1]):
        for yidx in range (z.shape[0]):
            if abs(z[yidx, xidx]-Zpoint)/Zpoint< error:
                x_acpt.append(xidx)
                y_acpt.append(yidx)
                z_acpt.append(z[yidx, xidx])
    return np.array(x_acpt), np.array(y_acpt), np.array(z_acpt)
#%%

#Density= [10]
#PWidth= [50]
Density= [2, 5, 8, 10, 20, 50, 80]
PWidth= [30, 35, 40, 45, 50]
DenError= np.zeros((len(Density), len(PWidth)))
PWError= np.zeros((len(Density), len(PWidth)))

DenC= 0
PWC= 0

for den in Density:
    for pw in PWidth:
        try:
            
            #import shadowgraphy fitting parameter
            #FirstPeakWidth= FirstPeakPara[0]*ln(FirstPeakPara[1]*Den)
            FirstPeakPara= np.load('FirstPeakWidthVSDen_Para.npy')
            #%
            #import afterglow fitting parameter
            #The loaded one den is in 1e16
            #Slope= AGSlopePara[0]*Den+AGSlopePara[1]
            AGSlopePara= np.load('p_slope.npy')
            #Offset= AGOffsetPara[0]*Den+AGOffsetPara[1]
            AGOffsetPara= np.load('p_offset.npy')
        
            #%
            #import shadowgraphy data
            ShaData= np.load('/media/valentinalee/TRANSCEND/ShadowgraphyResults/FFResults2/PW'\
                             +str(pw)+'_GP5_n'+str(den)+'.npy')
            #pixel
            x_pixel= 15000/2**14*1e-6
            y_pixel= 3200/2**8*1e-6
            #import afterglow data
            AGData= np.load('/media/valentinalee/TRANSCEND/AfterglowFinalPowerFix/ScanResults/AfterglowScan/n'\
                            +str(den)+'_GP5_PW'+str(pw)+'_map.npy')
            AGData_x= np.load('/media/valentinalee/TRANSCEND/AfterglowFinalPowerFix/ScanResults/AfterglowScan/n'\
                              +str(den)+'_GP5_PW'+str(pw)+'_x.npy')
            
            #%
            #measure the first peak width
            FirstPeakWidthM, PTValueM= FindShadowgraphyMrm(abs(ShaData[:, int(ShaData.shape[1]/2)])**2, y_pixel)
            #measure the afterglow width
            IntAG= AGData[:, int(AGData.shape[1]-1)]
            xvar= AGData_x[int(len(AGData_x)/2):len(AGData_x)]
            yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
            fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])
            
            errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
            p0= [1e11, 20, 2];
            p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
            AGWidthM= p1[1]
            #%
            #creat 2D array of shadowgraphy peak-trough results
            PT2D= np.load('PTValue_Results2D.npy')
            PT_x= np.load('PTValue_Results2D_x.npy')
            PT_y= np.load('PTValue_Results2D_y.npy')
            
            #fill with nan with zeros
            PT2D= np.nan_to_num(PT2D)
            
            #interpolation to high grid
            f= interpolate.interp2d(PT_x, PT_y, PT2D, kind='cubic')
            newPT_x= np.linspace(min(PT_x), max(PT_x), 2000)
            newPT_y= np.linspace(min(PT_y), max(PT_y), 2000)
            
            newPT2D= f(newPT_x, newPT_y)
            #%
            #find initial den from shadowgraphy
            intden= np.exp(FirstPeakWidthM/FirstPeakPara[0])/FirstPeakPara[1]*10
            #find initial pw from afterglow & den
            intSlp= AGSlopePara[0]*intden+AGSlopePara[1]
            intOfst= AGOffsetPara[0]*intden+AGOffsetPara[1]
            intpw= (AGWidthM-intOfst)/intSlp
            
            #use xyFindz to find z find PT value
            x_acptidx, y_acptidx, z_acpt= AllPointAtThisZ(newPT2D, PTValueM, error= 0.05)
            GoodIdx= np.argmin(np.sqrt((newPT_x[[x_acptidx]]-intpw)**2+(newPT_y[[y_acptidx]]-intden)**2))
            
            denH= newPT_y[y_acptidx[GoodIdx]]
            pwH= newPT_x[x_acptidx[GoodIdx]]
            #make a leastsq fitting where at that den & pw, 
            #sqrt(first peak witdh error)* sqrt(PT value error) * (afterglow error) is the smallest 
            #return den and pw
            #%
            avgDen= (intden+denH)/2
#            avgSlp= AGSlopePara[0]+AGSlopePara[1]*avgDen+AGSlopePara[2]*avgDen**2+AGSlopePara[3]*avgDen**3
#            avgOfst= AGOffsetPara[0]+AGOffsetPara[1]*avgDen+AGOffsetPara[2]*avgDen**2+AGOffsetPara[3]*avgDen**3
#            avgPW= (AGWidthM-avgOfst)/avgSlp
            avgPW= (intpw+pwH)/2
            
#            DenE= abs(den-denH)/den
#            PWE= abs(pw-pwH)/pw
            DenE= abs(den-avgDen)/den
            PWE= abs(pw-avgPW)/pw
            DenError[DenC, PWC]= DenE
            PWError[DenC, PWC]= PWE
            PWC= PWC+1
            print(avgDen)
            print(avgPW)
            print(denH)
            print(pwH)

        except:
            pass
    DenC= DenC+1
    PWC= 0
    
#%%
Density= ['0', '2', '5', '8', '10', '20', '50', '80']
PWidth= ['0', '30', '35', '40', '45', '50']
#plt.pcolormesh(DenError)
plt.pcolormesh(PWidth, Density, DenError)
plt.xlabel('PlasmaWidth')
plt.ylabel('Density (1e16)')
plt.colorbar()
#%%
plt.figure(2)
plt.pcolormesh(PWidth, Density, PWError)
plt.xlabel('PlasmaWidth')
plt.ylabel('Density')
plt.colorbar()
#%%
def fitfn (var):
    den= var[0]
    pw= var[1]
    PTvalueC= xyFINDz(newPT_x, newPT_y, newPT2D, intpw, intden)
    FirstPeakWidthC= FirstPeakPara[0]*np.log(FirstPeakPara[1]*den)
    Slope= AGSlopePara[0]*den+AGSlopePara[1]
    Offset= AGOffsetPara[0]*den+AGOffsetPara[1]
    AGWidthC= Slope*pw+ Offset
#    print(np.sqrt(abs(FirstPeakWidthC-FirstPeakWidthM)))
#    print(np.sqrt(abs(PTvalueC-PTvalueM)))
#    print(abs(AGWidthC-AGWidthM))
    
    return (abs(FirstPeakWidthC-FirstPeakWidthM)/2)+(abs(PTvalueC-PTvalueM)/2)+abs(AGWidthC-AGWidthM)
#%%
p0= [den, pw]
#p0= np.array([den, pw])
errorFn= lambda p: fitfn(p)
p1, success = optimize.leastsq(errorFn, p0[:], epsfcn=2e-6)
print(p1)
