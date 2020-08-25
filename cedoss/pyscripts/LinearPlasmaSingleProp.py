#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 10:15:27 2020

Single particle tracking of the transfer matrix for a linear gradient plasma source
with the approximate expression.

@author: chris
"""
import numpy as np
import matplotlib.pyplot as plt
"""
#all units in cgs
c = 2.9979e10
e = 4.8032e-10
me = 9.1094e-28
pi = 3.14159

gamma = 19569.5
n = 3e16
sl = -9.91e16
rp = 43.55e-6 * 1e2

wp = np.sqrt(4*pi*n*e**2/me)
kb = wp/np.sqrt(2*gamma*c**2)
kr = (pi*e**2*sl*rp**2)/(gamma*c**2*me)

step = 0.1e-6 * 1e2
flat = 0.10 * 1e2
y0 = 3e-6 * 1e2
yp0 = 0

s_arr = np.linspace(0, flat, int(flat/step))
y = np.zeros(len(s_arr)); y[0]=y0
yp = np.zeros(len(s_arr)); yp[0]=yp0
y_uni = np.copy(y); yp_uni = np.copy(yp)

loop_arr = np.array(range(len(s_arr)-1))+1
for i in loop_arr:
    if i%100000==0:
        print(i/len(loop_arr)*100,"%")
    y[i]  = y[i-1]*np.cos(kb*step)       + yp[i-1]/kb*np.sin(kb*step) + 2*kr/kb**2*np.square(np.sin(kb*step/2))
    yp[i] = y[i-1]*(-kb)*np.sin(kb*step) + yp[i-1]*np.cos(kb*step)    + kr/kb*np.sin(kb*step)
    
    y_uni[i]  = y_uni[i-1]*np.cos(kb*step)       + yp_uni[i-1]/kb*np.sin(kb*step)
    yp_uni[i] = y_uni[i-1]*(-kb)*np.sin(kb*step) + yp_uni[i-1]*np.cos(kb*step)
"""
plt.figure(figsize=(12,6))
plt.title("Modified Betatron Motion due to Linear Plasma Gradient")
plt.ylabel("y "+r'$(\mu m)$')
plt.xlabel("s "+r'$(cm)$')
plt.plot(s_arr,y*1e4,c='red',ls='solid',label="Linear Gradient")
plt.plot(s_arr,y_uni*1e4,c='blue',ls='dashed',label="Uniform")
plt.axvline(2*pi/kb,c='black',ls='dotted',label=r'$2n\pi/k_\beta$')
plt.axvline(2*2*pi/kb,c='black',ls='dotted')
plt.grid(); plt.legend(); plt.show()

plt.figure(figsize=(12,6))
plt.title("Modified Betatron Motion due to Linear Plasma Gradient")
plt.ylabel("y' "+r'$(mrad)$')
plt.xlabel("s "+r'$(cm)$')
plt.plot(s_arr,yp*1e3,c='red',ls='solid',label="Linear Gradient")
plt.plot(s_arr,yp_uni*1e3,c='blue',ls='dashed',label="Uniform")
plt.axvline(2*pi/kb,c='black',ls='dotted',label=r'$2n\pi/k_\beta$')
plt.axvline(2*2*pi/kb,c='black',ls='dotted')
plt.grid(); plt.legend(); plt.show()