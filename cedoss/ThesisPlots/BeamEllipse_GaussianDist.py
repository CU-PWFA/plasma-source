#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:37:29 2023

Functions to plot Beam Ellipses for thesis Figures

This Version populates the beam ellispe with a Gaussian distribution

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

geo_emit = 1.0
alpha = -0.7
beta = 2.2

gamma = (1 + alpha**2)/beta

N = int(2e4)
x1r = np.random.uniform(0,1,N)
x2r = np.random.uniform(0,1,N)
Jx = -geo_emit * np.log(x1r)
phix = 2*np.pi*x2r
ux = np.sqrt(2*Jx)*np.cos(phix)
vx = -np.sqrt(2*Jx)*np.sin(phix)

ptcls_x = ux*np.sqrt(beta)
ptcls_xp = (vx-alpha*ux)/np.sqrt(beta)








x=np.linspace(-1,1,50000)*np.sqrt(beta*geo_emit)
func1 = (np.sqrt(beta*geo_emit-np.square(x))-alpha*x)/beta
func2 = (-np.sqrt(beta*geo_emit-np.square(x))-alpha*x)/beta

factor = 6
x_95=np.linspace(-1,1,50000)*np.sqrt(beta*geo_emit*factor)
func1_95 = (np.sqrt(beta*geo_emit*factor-np.square(x_95))-alpha*x_95)/beta
func2_95 = (-np.sqrt(beta*geo_emit*factor-np.square(x_95))-alpha*x_95)/beta

fig = plt.figure(figsize=(4,4))
ax = plt.Axes(fig,[0,0,1,1])
#ax.set_axis_off()
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['top'].set_position('zero')
ax.spines['right'].set_position('zero')

fig.add_axes(ax)

plt.scatter(ptcls_x,ptcls_xp,s=0.1,marker='.',c='deepskyblue')
plt.plot(x,func1,c='k')
plt.plot(x,func2,c='k')
plt.plot(x_95,func1_95,c='k',ls='dashed')
plt.plot(x_95,func2_95,c='k',ls='dashed')


extension = 3.9

plt.xlim([-np.sqrt(beta*geo_emit)*extension,np.sqrt(beta*geo_emit)*extension])
plt.ylim([-np.sqrt(gamma*geo_emit)*extension,np.sqrt(gamma*geo_emit)*extension])

plt.text(0.9*extension*np.sqrt(beta*geo_emit),0.1*np.sqrt(gamma*geo_emit),r'$x$')
plt.text(0.1*np.sqrt(beta*geo_emit),0.9*extension*np.sqrt(gamma*geo_emit),r"$x'$")

ax.set_yticklabels([])
ax.set_xticklabels([])


fun1 = func1[~np.isnan(func1)]
x_maxxp = x[np.argmax(fun1)]
y_maxxp = np.max(fun1)
x_maxx = x[-1]
y_maxx = fun1[-1]

plt.plot([0,x_maxx],[0,y_maxx],c='k',ls='dotted')
plt.plot([x_maxx,x_maxx],[0,0.9*y_maxxp],c='k',ls='dashed')
plt.text(-0.4*x_maxx,y_maxxp,r"$\sigma_{x'}$",color='k')
#plt.text(x_maxxp*0.6,y_maxxp*1.1,"Slope:"r'$-\gamma/\alpha$',color='k')

plt.plot([0,x_maxxp],[0,y_maxxp],c='k',ls='dotted')
plt.plot([0,0.9*x_maxx],[y_maxxp,y_maxxp],c='k',ls='dashed')
plt.text(0.9*x_maxx,-0.3*y_maxxp,r'$\sigma_x$',color='k')
#plt.text(x_maxx*1.1,y_maxx,"Slope:"r'$-\alpha/\beta$',color='k')

fun2 = func2[~np.isnan(func2)]
fun2_x = x[np.argmin(fun2)]*1.4
fun2_y = np.min(fun2)*1.2
plt.text(fun2_x,fun2_y,r'$\epsilon_{RMS}$')

fun3 = func2_95[~np.isnan(func2_95)]
fun3_x = x_95[np.argmin(fun3)]*1.3
fun3_y = np.min(fun3)*1.1
plt.text(fun3_x,fun3_y,r'$\epsilon_{95\%}$')

plt.show()