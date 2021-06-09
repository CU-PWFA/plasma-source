#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 12:52:26 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy import optimize

#%%
NewTimeGrid= 10000
alpha= 2.36e-20 #m^3/s


density= 18 
power= 5
pw= 50
            
PlasmaExpan= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy') #plasma density in r-t
Ext= np.load('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_SimExt.npy')
Box= Ext[0] #The Simulation goes from -Box to Box micron
GridNum= int(Ext[1]) #Number of dots in Box
SimTime= Ext[2] #Simualtion time in s
FrameNum= int(Ext[3]) #number of time frame in SimTime
FnlTimeStepNum= int(FrameNum-1) #Time frame-1

x= np.linspace(-Box, Box, GridNum) #x in space
y= np.linspace(-Box, Box, GridNum) #y in space
r= x[int(GridNum/2):GridNum-1] #r is the positive x 
#I don't know why r has to miss one grid but it has something with even or odd points
X, Y= np.meshgrid(x, y)  
R= np.sqrt(X**2+Y**2)

LineoutTime= PlasmaExpan[int(PlasmaExpan.shape[0]/2):int(PlasmaExpan.shape[0]-1), :] #density in positve r-t

time= np.linspace(0, SimTime, FrameNum) #time grid

f= interp2d(time, r, LineoutTime, kind='cubic')
Time= np.linspace(0, SimTime, NewTimeGrid)
LineOutTime= f(Time, r)
#The above three lines are increasing the time grid density 

#%%
PhotonNum= np.zeros((int(GridNum/2-1), int(NewTimeGrid-1)))#same grid size as LineOutTime but time dimention is one grid less
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
#%%  
AfterglowTime= np.empty((int(GridNum), int(FnlTimeStepNum)))
for sec in range (0, int(FnlTimeStepNum)):
    IntPhoNum= np.empty(int(GridNum/2-1))
    for f in range (0, int(GridNum/2-1)):
        frame= PhotonNum[f, :] #photon-time at given location (r)
        IntPhoNum[f]= np.sum(frame[0:int(NewTimeGrid/FnlTimeStepNum)+sec*int(NewTimeGrid/FnlTimeStepNum)])
        #Intergrate up to sec, this time

#The next few line until "Mycroft" do fitting of IntPhoNum at each time frame
#do a fitting first of a half 1d Gaussian
#Then take the variable to plot at R which is in 2D polar cordinate 
    xvar= r 
    yvar= IntPhoNum
    fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**p[2]/(2*p[1]**p[2])))

    errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
    p0= [5e7, 10, 2]
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))

    AfterglowMap= p1[0]*np.exp(-(R**p1[2]/(2*p1[1]**p1[2])))
#Here, AfterglowMap is a 2D afterglow map looking transversly at each time frame
#"Mycroft"

    AfterglowSide= np.sum(AfterglowMap, axis= 0) # Do a side integration to get a side view
    AfterglowTime[:, sec]= AfterglowSide #Save to AfterglowTime which is a r-integrated time plot

Timeplt= np.linspace(0, SimTime*1e9, int(FnlTimeStepNum))
#np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x', x)
#np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_t', Timeplt)
#np.save('AfterglowScan/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map', AfterglowTime)

plt.pcolormesh(Timeplt, x, AfterglowTime)
plt.xlabel('Integration time (ns)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= '# of photons')
plt.title('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw))
plt.savefig('AfterglowPlots/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_AfterglowMap')
#plt.close()

