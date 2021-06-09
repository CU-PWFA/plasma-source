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
alpha= 3e-16#5e-18 #m^3/s

for density in range (1, 10, 4):
    for power in range (1, 11):
        for pw in range (30, 160, 10):
            
            CheckRun= 'AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy'
            if path.exists(CheckRun) == False:
            
                File= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy'
                print(File)
                if path.exists(File)==True:
                    PlasmaExpan= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy')
                    Ext= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_SimExt.npy')
                    Box= Ext[0]
                    GridSize= Ext[1]
                    SimTime= Ext[2]
                    Frame= Ext[3]
                    FnlTimeStepNum= int(Frame-1)
        
                    GridNum= PlasmaExpan.shape[0]
                    FrameNum= PlasmaExpan.shape[1]
        
                    x= np.linspace(-Box, Box, GridNum)
                    y= np.linspace(-Box, Box, GridNum)
                    r= x[int(GridNum/2):GridNum-1]
                    X, Y= np.meshgrid(x, y)
                    R= np.sqrt(X**2+Y**2)
        
                    LineoutTime= PlasmaExpan[int(PlasmaExpan.shape[0]/2):int(PlasmaExpan.shape[0]-1), :]
                    time= np.linspace(0, SimTime, FrameNum)
                    f= interp2d(time, r, LineoutTime, kind='cubic')
                    Time= np.linspace(0, SimTime, NewTimeGrid)
                    LineOutTime= f(Time, r)
        
                    Ntimestep= int(Time.shape[0])
                    PhotonNum= np.zeros((int(GridNum/2-1), int(NewTimeGrid-1)))
                    for f in range (0, int(GridNum/2-1)):
                        den= LineOutTime[f, :]
                        n0=den[0]
                        dndt=-alpha*n0**2
                        dt= np.amax(Time)/(Ntimestep-1)
                        dn= dndt*dt
                        nNew=np.zeros(Ntimestep-1)
                        pn= np.zeros(Ntimestep)
                        Den=n0
                        for t in range (0, Ntimestep-1):
                            pn[t]=dn 
                            n= Den+(den[t+1]-den[t])+dn
                            dndt=-alpha*n**2
                            dn= dndt*dt
                            if n>0 :
                                nNew[t]=n
                                Den= n
                            else:
                                nNew[t]=0
                                Den= 0
                
                        PhotonNum[f, :]= -pn[1:int(NewTimeGrid)]
               
                    AfterglowTime= np.zeros((int(GridNum), int(FnlTimeStepNum)))
                    for sec in range (0, int(FnlTimeStepNum)):
                        IntPhoNum= np.zeros(int(GridNum/2-1))
                        for f in range (0, int(GridNum/2-1)):
                            frame= PhotonNum[f, :]
                            IntPhoNum[f]= np.sum(frame[0:int(NewTimeGrid/FnlTimeStepNum)+sec*int(NewTimeGrid/FnlTimeStepNum)])
        
                        xvar= r
                        yvar= IntPhoNum
                        fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)
        
                        errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
                        p0= [10e20, 200, 0];
                        p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
                        AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
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