#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:31:07 2020

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;

#%%
lam0=800e-9 
c=3e8
k0= 2*np.pi/lam0
w0= c/lam0*2*np.pi
t= np.linspace(-500e-12, 500e-12, 5000)
#chirpR= 3e-15/100e-12
chirpR= -1e20
ft= 3.5e-6*t+2.65e-15
TW= 200e-12
Et= np.exp(-t**2/TW**2)*np.exp(-1j*chirpR*t**2)*np.exp(1j*w0*t)
#%%
plt.plot(t, abs(Et)**2)
#plt.plot(t, Et)
plt.xlabel('t (s)')
plt.ylabel('Normalized Intensity')

#%%
dz= 0.1
#k= 
Eprop= ifft(fft(Et)*np.exp(1j*k*dz))

#%%
dt= 200e-12
dx= dt*c
d= 1e-3/1200
dlam= 100e-9
order= 1

for thetai in np.linspace(0, np.pi/2, 100):
    
    for ang in np.linspace(np.pi/6, np.pi, 100):
        
        for L in np.linspace(0.1, 1, 40):
                
            thetad= np.arcsin(lam0/d-np.sin(thetai))
            dthetad= dlam*order/(d*np.cos(thetad))
            
            tdb= ang
            tdr= np.pi-ang-dthetad
            if tdr>np.pi/2:
                wb= (L/(np.tan(tdb)))
                wr= (L/(np.tan(np.pi-tdr)))
                w= wb-wr
                lb= L/np.sin(tdb)
                lr= L/np.sin(np.pi-tdr)
            
            
            else:
                wb= (L/(np.tan(tdb)))
                wr= (L/(np.tan(tdr)))
                w= wb+wr
                lb= L/np.sin(tdb)
                lr= L/np.sin(tdr)
            
            phi1= np.arcsin(lam0/d-np.sin(abs((np.pi/2-tdb))))
            phi2= np.arcsin(lam0/d-np.sin(abs(np.pi/2-tdr)))
            totb=lb+w*np.sin(phi1)
            totr=lr
            if abs(phi1-phi2)<0.005:
                if (abs(abs(totb-totr)*2-dx))<(0.05*dx):
                    
                    print(thetai, ang, L)
                
                
#%%
thetai= 0.47599888690754444
ang= 1.5020409320193622
L= 0.12307692307692308
d= 1e-3/1200
dlam= 100e-9
order= 1
thetad= np.arcsin(lam0/d-np.sin(thetai))
dthetad= dlam*order/(d*np.cos(thetad))

tdb= ang
tdr= np.pi-ang-dthetad
if tdr>np.pi/2:
    wb= (L/(np.tan(tdb)))
    wr= (L/(np.tan(np.pi-tdr)))
    w= wb-wr
    lb= L/np.sin(tdb)
    lr= L/np.sin(np.pi-tdr)


else:
    wb= (L/(np.tan(tdb))) 
    wr= (L/(np.tan(tdr)))
    w= wb+wr
    lb= L/np.sin(tdb)
    lr= L/np.sin(tdr)

phi1= np.arcsin(lam0/d-np.sin(abs(np.pi/2-tdb)))
phi2= np.arcsin(lam0/d-np.sin(abs(np.pi/2-tdr)))
totb=lb+w*np.sin(phi1)
totr=lr
print(phi1, phi2, phi1-phi2, totb, totr, (totb-totr)*2)