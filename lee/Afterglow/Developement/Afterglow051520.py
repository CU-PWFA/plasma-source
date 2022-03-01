vvvv#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 17:41:18 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp2d
from scipy import optimize

#%%
ELECMASS = 9.10938356e-31     # kG
PROTMASS = 1.672621898e-27    # kG
#%% Call data
def qReader(FileName):

    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    q=np.array(Flu.get('q'))
    Den=q[:, 0]
    Eng=q[:, 4]
    return Den, Eng
#%%
BasePath= 'AfterglowStudy/DensityStudy/1e22/'
BaseFileName= 'cap_2D_'
#%%
GridNum= 150
Box= 600#1500
FrameNum= 30
SimTime= 100e-9
NewTimeGrid= 10000
IntResolution= 1e-9#s
FnlTimeStepNum= round(SimTime/IntResolution)

x= np.linspace(-Box, Box, GridNum)
y= np.linspace(-Box, Box, GridNum)
r= x[int(GridNum/2):GridNum]
X, Y= np.meshgrid(x, y)
R= np.sqrt(X**2+Y**2)
#%%
MassDen={}
eng={}
den={}

den_3d=np.zeros((GridNum, GridNum, int(FrameNum+1)))
en_3d=np.zeros((GridNum, GridNum, int(FrameNum+1)))

for n in range (0, int(FrameNum+1)):
    FileName= BaseFileName+str(n)+'.h5'
    path= BasePath+ FileName
    MassDen[n]= np.reshape(qReader(path) [0], (GridNum, GridNum))
    den[n]=MassDen[n]/(ELECMASS+PROTMASS*4)
    den_3d[:, :, n]= den[n]
    eng[n]= np.reshape(qReader(path) [1], (GridNum, GridNum))
    en_3d[:, :, n]= eng[n]
#%%
LineoutTime= np.zeros((int(GridNum/2), int(FrameNum+1)))
for n in range (0, int(FrameNum+1)):
    LineoutTime[:, n]= den[n][int(GridNum/2):int(GridNum), int(GridNum/2)]
#%%
time= np.linspace(0, SimTime, int(FrameNum+1))
f= interp2d(time, r, LineoutTime, kind='cubic')
Time= np.linspace(0, SimTime, NewTimeGrid)
LineOutTime= f(Time, r)
#%%
plt.figure(0)
plt.pcolormesh(Time, r, LineOutTime)
plt.xlabel('time frame')
plt.ylabel('r')
plt.title('Density Lineout - Time')

 #%%create a 1D density decay plot
Ntimestep= int(Time.shape[0])
#alpha= 5e-18 #m^3/s
alpha= 3e-16#5e-18 #m^3/s
PhotonNum= np.zeros((int(GridNum/2), int(NewTimeGrid-1)))
for f in range (0, int(GridNum/2)):
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
       
#%%
plt.pcolormesh(PhotonNum)
plt.xlabel('time frame')
plt.ylabel('r')
plt.title('Radial Photon Number - Time')
#%%Plot the recombination density change
plt.plot(Time, den, (Time+dt)[0:9999], nNew)
plt.xlabel('s')
plt.ylabel('m^-3')
#plt.yscale('log')

#plt.title('Alpha=1e-9cm^3/s_TimeStep~1000')

#%%Run from "Here"
#%%
'''
IntPhoNum= np.zeros(int(GridNum/2))
for f in range (0, int(GridNum/2)):
    frame= PhotonNum[f, :]
    IntPhoNum[f]= np.sum(frame[0:30])
#%%
plt.plot(r, IntPhoNum)
plt.xlabel('r ($\mu$m)')
plt.ylabel('Photon Number')
#%%
xvar= r
yvar= IntPhoNum
fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
p0= [10e20, 200, 0];
p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))#, epsfcn= 2.5e-34);
print(p1)
#%%
plt.plot(xvar, yvar, '.', label= 'plasma data')
plt.plot(xvar, fitfn(p1, xvar), '-', label='fitting')
plt.legend()

#%%
AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
#%%
plt.pcolormesh(x, y, AfterglowMap)
plt.axis('scaled')
plt.xlabel('$\mu$m')
plt.ylabel('$\mu$m')

#%%
AfterglowSide= np.sum(AfterglowMap, axis= 0)
plt.plot(AfterglowSide)
'''
#%%"Here"
#%%Afterglow at each time sec
AfterglowTime= np.zeros((int(GridNum), int(FnlTimeStepNum)))
for sec in range (0, int(FnlTimeStepNum)):
    IntPhoNum= np.zeros(int(GridNum/2))
    for f in range (0, int(GridNum/2)):
        frame= PhotonNum[f, :]
#        IntPhoNum[f]= np.sum(frame[0+sec*int(FnlTimeStepNum):int(FnlTimeStepNum)+sec*int(FnlTimeStepNum)])
        IntPhoNum[f]= np.sum(frame[0+sec*int(NewTimeGrid/FnlTimeStepNum):int(NewTimeGrid/FnlTimeStepNum)+sec*int(NewTimeGrid/FnlTimeStepNum)])

    xvar= r
    yvar= IntPhoNum
    fitfn= lambda p, xvar: p[0]*np.exp(-((xvar-p[2])**2/(2*p[1]**2))**1)

    errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
    p0= [10e20, 200, 0];
    p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    AfterglowMap= p1[0]*np.exp(-((R-p1[2])**2/(2*p1[1]**2))**1)
    AfterglowSide= np.sum(AfterglowMap, axis= 0)
    AfterglowTime[:, sec]= AfterglowSide

#%%
plt.figure(2)
plt.pcolormesh(AfterglowTime)
#%%Afterglow Integral over each time
AfterglowTime= np.zeros((int(GridNum), int(FnlTimeStepNum)))
for sec in range (0, int(FnlTimeStepNum)):
    IntPhoNum= np.zeros(int(GridNum/2))
    for f in range (0, int(GridNum/2)):
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
#%%
SAfterglow= AfterglowTime[50:100, :]
y_= y[50:100]
#%%
plt.figure(2)
Timeplt= np.linspace(0, SimTime*1e9, int(FnlTimeStepNum))
plt.pcolormesh(Timeplt, y_, SAfterglow)
plt.xlabel('Integration time (ns)')
plt.ylabel('Lineout ($\mu$m)')
plt.colorbar(label= '# of photons')
plt.title('w= 50 $\mu$m, n=1e16cm^-3')
#%%
#%%
Lineout= SAfterglow[:, 99]#/np.amax(SAfterglow[:, 99])
yI= np.linspace(np.amin(y_), np.amax(y_), 1000)
InitialPlasma= 1e17*np.exp(-(yI**2/(50)**2)**2.6)
#%%
AGdata= AfterglowTime[:, 99]
np.save('AGdata.npy', AGdata)
#%%
plt.title('w=50 Normalized Density/Intensity')
plt.plot(yI, InitialPlasma, label= 'Initial Plasma')
#plt.plot(y_, Lineout, label='Afterglow')
plt.xlabel('Lineout ($\mu$m)')
plt.legend()
#%%
Lineout1e16=Lineout
#Lineout5e16=Lineout
#Lineout1e17=Lineout
#Itl1e16= InitialPlasma
#Itl5e16= InitialPlasma
#Itl1e17= InitialPlasma
#%%
#Lineout2_6= Lineout
#Lineout1_8= Lineout
#Lineout1_0= Lineout
Itl2_6= InitialPlasma
#Itl1_8= InitialPlasma
#Itl1_0= InitialPlasma

#%%
plt.title('Ramp Study \n Normalized Density/Intensity')
#plt.plot(yI, Itl1e17, 'b--', label= 'Initial Plasma ^2.6')
plt.plot(y_, Lineout1e17, 'b-', label='Afterglow ^2.6')
#plt.plot(yI, Itl1_8, 'g--', label= 'Initial Plasma ^1.8')
#plt.plot(y_, Lineout1_8, 'g-', label='Afterglow ^1.8')
#plt.plot(yI, Itl1_0, 'y--', label= 'Initial Plasma ^1')
#plt.plot(y_, Lineout1_0, 'y-', label='Afterglow ^1')

plt.legend()

#%%
fig, ax1 = plt.subplots()
color1 = 'tab:red'
ax1.set_xlabel('r ($\mu$m)')
ax1.set_ylabel('Plasma Density (cm^-3)', color=color1)
ax1.plot(yI, Itl1e17, ':', label= 'Initial Plasma 1e17', color=color1)
ax1.plot(yI, Itl5e16, '--', label= 'Initial Plasma 5e16', color=color1)
ax1.plot(yI, Itl1e16, '-', label= 'Initial Plasma 1e16', color=color1)
ax1.legend(loc=2)
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()
color2 = 'tab:blue'
ax2.set_ylabel('# of photons', color=color2)  
ax2.plot(y_, Lineout1e17, ':', label='Afterglow 1e17', color=color2)
ax2.plot(y_, Lineout5e16, '--', label='Afterglow 5e16', color=color2)
ax2.plot(y_, Lineout1e16, '-', label='Afterglow 1e16', color=color2)
ax2.tick_params(axis='y', labelcolor=color2)
ax2.legend(loc=0)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

