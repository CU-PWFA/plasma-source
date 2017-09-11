#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 11:19:02 2017

@author: chris
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

def func(x, y):
    return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

min_x = -1
min_y = -1
max_x = 2
max_y = 1
nx = 2**8
ny = 2**7

bgr = 0
numpts = 3000

grid_x, grid_y = np.mgrid[min_x:max_x:nx*1j, min_y:max_y:ny*1j]

points = np.random.rand(numpts, 2)

for i in range(len(points)):
    points[i] = (points[i]*[(max_x-min_x),(max_y-min_y)])+[min_x,min_y]
    
values = func(points[:,0], points[:,1])

grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear', fill_value = bgr)
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic', fill_value = bgr)

plt.subplot(221)
plt.imshow(func(grid_x, grid_y).T, extent=(min_x,max_x,min_y,max_y), origin='lower')
plt.plot(points[:,0], points[:,1], 'k.', ms=1)
plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(min_x,max_x,min_y,max_y), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(min_x,max_x,min_y,max_y), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(min_x,max_x,min_y,max_y), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()

