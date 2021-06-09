#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 13:08:13 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
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
BasePath= 'AfterglowStudy/RampStudy/^2.6/'
BaseFileName= 'cap_2D_'
#%%
GridNum= 150
Box= 600#1500
FrameNum= 30
SimTime= 100
NewTimeGrid= int(GridNum/2)
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
plt.plot(LineoutTime[5])
#%%
plt.figure(2)
plt.pcolormesh(Time, r, LineOutTime/np.amax(LineOutTime))
plt.xlabel('time frame')
plt.ylabel('r')
plt.title('Density Lineout - Time')
plt.colorbar()
#%%Polynomial 
xvar, yvar= np.meshgrid(Time, r)
zvar= LineOutTime/np.amax(LineOutTime)
fitfn= lambda p, xvar, yvar: p[0]+p[1]*xvar+p[2]*yvar+p[3]*xvar**2+p[4]*yvar**2+p[5]*xvar*yvar+p[6]*xvar**2*yvar+p[7]*xvar*yvar**2 \
+p[8]*xvar**2*yvar**2+p[9]*xvar*yvar**3+p[10]*xvar**3*yvar

errfunc = lambda p: np.ravel(fitfn(p, xvar, yvar)) - np.ravel(zvar);
p0= [1, -0.1, -0.025, 0, 0, 0, 0, 0, 0, 0, 0];
p1, success = optimize.leastsq(errfunc, p0[:], epsfcn= 2.5e-20);
print(p1)
#%%Asymetry high order gaussian+ rotation 
xvar, yvar= np.meshgrid(Time, r)
zvar= LineOutTime/np.amax(LineOutTime)
fitfn= lambda p, xvar, yvar: (p[0]*np.exp(-(xvar*np.cos(-p[8])-yvar*np.sin(-p[8]))**2/(2*p[1]**2))+\
                              p[2]*np.exp(-((xvar*np.cos(-p[8])-yvar*np.sin(-p[8]))**3/(2*p[3])**3)))* \
                              (p[4]*np.exp(-(xvar*np.sin(-p[9])+yvar*np.cos(-p[9]))**2/(2*p[5]**2))+\
                               p[6]*np.exp(-((xvar*np.sin(-p[10])+yvar*np.cos(-p[10]))**3/(2*p[7])**3)))

errfunc = lambda p: np.ravel(fitfn(p, xvar, yvar)) - np.ravel(zvar);
p0= [1, 5, 1, 50, 1, 5, 1, 50, np.pi/4, np.pi/4, np.pi/4, np.pi/4];
p1, success = optimize.leastsq(errfunc, p0[:], epsfcn= 2.5e-20);
print(p1)
#%%Asymetry high order gaussian+ tanh rotation 
#Winner!
xvar, yvar= np.meshgrid(Time, r)
zvar= LineOutTime/np.amax(LineOutTime)
fitfn= lambda p, xvar, yvar:  (p[0]*np.exp(-(xvar*np.cos(-p[5])-yvar*np.sin(-p[5]))**2/(2*p[1]**2))+\
                              p[2]*np.exp(-((xvar*np.cos(-p[6])-yvar*np.sin(-p[6]))**3/(2*p[3]**2)**3)))* \
                              ((-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[4])+1)/np.amax(-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[4])+1))

errfunc = lambda p: np.ravel(fitfn(p, xvar, yvar)) - np.ravel(zvar);
p0= [1, 50, 1, 1, 50, 20, np.pi/4, np.pi/4];
p1, success = optimize.leastsq(errfunc, p0[:], epsfcn= 2.5e-20);
print(p1)
#%%tanh+ tanh rotation 
xvar, yvar= np.meshgrid(Time, r)
zvar= LineOutTime/np.amax(LineOutTime)
fitfn= lambda p, xvar, yvar:  ((-np.tanh((xvar*np.cos(-p[2])-yvar*np.sin(-p[2]))/p[0])+1)/np.amax(-np.tanh((xvar*np.cos(-p[2])-yvar*np.sin(-p[2]))/p[0])+1)+\
                              (-np.tanh((xvar*np.cos(-p[3])-yvar*np.sin(-p[3]))/p[1])+1)/np.amax(-np.tanh((xvar*np.cos(-p[3])-yvar*np.sin(-p[3]))/p[1])+1))*\
                              ((-np.tanh((xvar*np.sin(-p[6])+yvar*np.cos(-p[6]))/p[4])+1)/np.amax(-np.tanh((xvar*np.sin(-p[6])+yvar*np.cos(-p[6]))/p[4])+1)+\
                              (-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[5])+1)/np.amax(-np.tanh((xvar*np.sin(-p[7])+yvar*np.cos(-p[7]))/p[5])+1))

errfunc = lambda p: np.ravel(fitfn(p, xvar, yvar)) - np.ravel(zvar);
p0= [20, 20,  np.pi/4, np.pi/4, 20, 20,  np.pi/4, np.pi/4];
p1, success = optimize.leastsq(errfunc, p0[:], epsfcn= 2.5e-20);
print(p1)

#%%
#Axes3D.scatter
plt.figure(1)
plt.pcolormesh(Time, r, fitfn(p1, xvar, yvar)/np.amax(fitfn(p1, xvar, yvar)))
plt.colorbar()
plt.title('tanh+tanh')
#%%Error calculation
fit= fitfn(p1, xvar, yvar)/np.amax(fitfn(p1, xvar, yvar))
Error= sum(sum((zvar-fit)/zvar))
print(Error)