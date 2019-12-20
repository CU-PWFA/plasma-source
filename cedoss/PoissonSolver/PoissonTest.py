#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 10:45:03 2019

Going to  try a Poisson solver to calculate the electric field
from some arbitrary charge density.

A little sketchy with the boundary condiitons, but we'll see how it goes.

Base from https://www.codeproject.com/articles/1087025/using-python-to-solve-computational-physics-proble

This version uses a charge density field to do actual Poisson

@author: chris
"""

# Simple Numerical Laplace Equation Solution using Finite Difference Method
import numpy as np
import matplotlib.pyplot as plt

# Set maximum iteration
maxIter = 500

# Set Dimension and delta
xwindow = ywindow = 1e-3 * 1e2 #cm
delta = 0.5e-6 * 1e2 * 10 #cm

# Boundary condition
Ttop = 0
Tbottom = 0
Tleft = 0
Tright = 0

# Initial guess of interior grid

e = 4.8032e-10 #statcoul
pi = np.pi
n_cent = 5e16 #cm-3

# Set colour interpolation and colour map
colorinterpolation = 50
colourMap = plt.cm.jet #you can try: colourMap = plt.cm.coolwarm

# Set meshgrid
X, Y = np.meshgrid(np.arange(-1/2*xwindow, 1/2*ywindow, delta), np.arange(-1/2*xwindow, 1/2*ywindow,delta))
lenX = len(X); lenY = len(Y)

# Set array size and set the interior value with Tguess
T = np.empty((lenX, lenY)); nplasma = np.empty((lenX, lenY))

for i in range(len(X)):
    for j in range(len(Y)):
        T[i][j] = -pi*e*n_cent*(np.square(X[i][j])+np.square(Y[i][j]))
nplasma.fill(n_cent)

# Set Boundary condition
T[(lenY-1):, :] = Ttop
T[:1, :] = Tbottom
T[:, (lenX-1):] = Tright
T[:, :1] = Tleft

# Iteration (We assume that the iteration is convergence in maxIter = 500)
print("Iteration started")
for iteration in range(0, maxIter):
    if iteration%50==0: 
        print(iteration/maxIter*100,"%")
    for i in range(1, lenX-2):
        for j in range(1, lenY-2):
            T[i, j] = 0.25 * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] + delta**2*4*pi*nplasma[i][j])

print("Iteration finished")

# Configure the contour
plt.title("Contour of Temperature")
plt.contourf(X, Y, T, colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()

print("")