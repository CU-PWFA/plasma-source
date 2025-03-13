#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:13:19 2023

Functions to plot Beam Ellipses for thesis Figures

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

geo_emit = 1.0
alpha = -0.77
beta = 1.2

gamma = (1 + alpha**2)/beta

x=np.linspace(-1,1,50000)*np.sqrt(beta*geo_emit)
func1 = (np.sqrt(beta*geo_emit-np.square(x))-alpha*x)/beta
func2 = (-np.sqrt(beta*geo_emit-np.square(x))-alpha*x)/beta

fig = plt.figure(figsize=(4,4))
ax = plt.Axes(fig,[0,0,1,1])
#ax.set_axis_off()
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['top'].set_position('zero')
ax.spines['right'].set_position('zero')

fig.add_axes(ax)
plt.plot(x,func1,c='k')
plt.plot(x,func2,c='k')

extension = 1.9

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

plt.plot([0,x_maxx],[0,y_maxx],c='b',ls='dotted')
plt.plot([x_maxx,x_maxx],[0,0.95*y_maxxp],c='b',ls='dashed')
plt.text(-0.3*x_maxx,y_maxxp,r'$\sqrt{\gamma\epsilon}$',color='r')
plt.text(x_maxxp*0.6,y_maxxp*1.1,"Slope:"r'$-\gamma/\alpha$',color='r')

plt.plot([0,x_maxxp],[0,y_maxxp],c='r',ls='dotted')
plt.plot([0,0.95*x_maxx],[y_maxxp,y_maxxp],c='r',ls='dashed')
plt.text(0.9*x_maxx,-0.2*y_maxxp,r'$\sqrt{\beta\epsilon}$',color='b')
plt.text(x_maxx*1.1,y_maxx,"Slope:"r'$-\alpha/\beta$',color='b')

plt.show()