#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 17:41:18 2020

@author: valentinalee
"""
#Run this in dir 'ScanResults/'
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy import optimize
from os import path

#%%
NewTimeGrid= 10000
alpha= 2.36e-20 #m^3/s

for density in range (1, 51):
    for power in range (1, 11):
        for pw in range (30, 160, 5):
            
            CheckRun= 'AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy'
            if path.exists(CheckRun) == False:
            
                File= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy'
                
                if path.exists(File)==True:
                    print(File)
                    PlasmaExpan= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy') 
                    Ext= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_SimExt.npy')
                    Box= Ext[0] 
                    GridNum= int(Ext[1]) 
                    SimTime= Ext[2] 
                    FrameNum= int(Ext[3]) 
                    FnlTimeStepNum= int(FrameNum-1) 
                    
                    x= np.linspace(-Box, Box, GridNum) 
                    y= np.linspace(-Box, Box, GridNum) 
                    r= x[int(GridNum/2):GridNum-1] #r is the positive x 
                    X, Y= np.meshgrid(x, y)  
                    R= np.sqrt(X**2+Y**2)
                    
                    LineoutTime= PlasmaExpan[int(PlasmaExpan.shape[0]/2):int(PlasmaExpan.shape[0]-1), :] #density in positve r-t
                    time= np.linspace(0, SimTime, FrameNum) #time grid
                    
                    f= interp2d(time, r, LineoutTime, kind='cubic')
                    Time= np.linspace(0, SimTime, NewTimeGrid)
                    LineOutTime= f(Time, r)
        

                    PhotonNum= np.zeros((int(GridNum/2-1), int(NewTimeGrid-1)))
                    dt= np.amax(Time)/(NewTimeGrid-1)
                    for f in range(LineOutTime.shape[0]):
                        den= LineOutTime[f, :]
                        n= den[0]
                        DN= np.zeros(LineOutTime.shape[1])
                        for t in range (LineOutTime.shape[1]):        
                            dn= -alpha*n**2*dt
                            n= den[t]+dn
                            DN[t]= dn
                        PhotonNum[f, :]=-DN[0:int(len(DN)-1)]



               
                    AfterglowTime= np.empty((int(GridNum), int(FnlTimeStepNum)))
                    for sec in range (0, int(FnlTimeStepNum)):
                        IntPhoNum= np.empty(int(GridNum/2-1))
                        for f in range (0, int(GridNum/2-1)):
                            frame= PhotonNum[f, :] 
                            IntPhoNum[f]= np.sum(frame[0:int(NewTimeGrid/FnlTimeStepNum)+sec*int(NewTimeGrid/FnlTimeStepNum)])
                        xvar= r 
                        yvar= IntPhoNum
                        fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**p[2]/(2*p[1]**p[2])))
                    
                        errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
                        p0= [5e7, 20, 2]
                        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
                    
                        AfterglowMap= p1[0]*np.exp(-(R**p1[2]/(2*p1[1]**p1[2])))
                    
                        AfterglowSide= np.sum(AfterglowMap, axis= 0) 
                        AfterglowTime[:, sec]= AfterglowSide 
                    Timeplt= np.linspace(0, SimTime*1e9, int(FnlTimeStepNum))
                    np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x', x)
                    np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_t', Timeplt)
                    np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map', AfterglowTime)
    
                    plt.pcolormesh(Timeplt, x, AfterglowTime)
                    plt.xlabel('Integration time (ns)')
                    plt.ylabel('Lineout ($\mu$m)')
                    plt.colorbar(label= '# of photons')
                    plt.title('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw))
                    plt.savefig('AfterglowPlots/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_AfterglowMap')
                    plt.close()
