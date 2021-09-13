#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:24:29 2021

Calculating the electric field from a ring of negative charge with a linear gradient in y.
Testing the hypothesis that my field error is due to the asymmetric sheaths

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

e = 4.8032e-10
pi = np.pi

"""
n_cent = -2e16  #cm-3
radius = 80.00e-6 *1e2 #cm
thickness = 3.0e-6*1e2 #cm
slope  = 8e17 #cm-4
xwindow = ywindow = 170e-6 * 1e2  #cm
deltax = deltay = 2e-6 * 1e2  #cm
"""
"""
n_cent = -6.775e16  #cm-3
radius = 80.00e-6 *1e2 #cm
thickness = 10.0e-6*1e2 #cm
slope  = 1.68e18 #cm-4
xwindow = ywindow = 182e-6 * 1e2  #cm
deltax = deltay = 4e-6 * 1e2  #cm
"""
"""
n_cent = -7.0635e16  #cm-3
radius = 80.00e-6 *1e2 #cm
thickness = 10.0e-6*1e2 #cm
slope  = 2.3e16 #cm-4
xwindow = ywindow = 182e-6 * 1e2  #cm
deltax = deltay = 2e-6 * 1e2  #cm
"""
"""
#JanGrad
n_cent = -3.605e16#-3.5e16  #cm-3
radius = 100.00e-6 *1e2 #cm
thickness = 12.5e-6*1e2 #cm
slope  = 1.03e18#1.0e18 #cm-4
xwindow = ywindow = 230e-6 * 1e2  #cm
deltax = deltay = 4.6e-6 * 1e2  #cm
"""

#OctGrad
n_cent = -7.0e16#-3.5e16  #cm-3
radius = 80.00e-6 *1e2 #cm
thickness = 10e-6*1e2 #cm
slope  = 2.3e17-2e17#5.0e17#1.0e18 #cm-4
xwindow = ywindow = 182e-6 * 1e2  #cm
deltax = deltay = 1e-6 * 1e2  #cm

X = np.arange(-1/2*xwindow, 1/2*xwindow, deltax)
Y = np.arange(-1/2*ywindow, 1/2*ywindow, deltay)

lenX = len(X); lenY = len(Y)

dengrid = np.zeros((lenX, lenY))
for i in range(len(X)):
    for j in range(len(Y)):
        if (np.sqrt(np.square(X[i])+np.square(Y[j])) > radius) & (np.sqrt(np.square(X[i])+np.square(Y[j])) < radius+thickness):
            dengrid[i][j] = n_cent + Y[j]*slope
        else:
            dengrid[i][j] = 0
            
plt.title("Density")
plt.imshow(np.transpose(dengrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.clim(0,n_cent-slope*radius)#2.6e16,3.4e16)
plt.show()

#sys.exit()

centX = centY = total = 0
for i in range(len(X)):
    for j in range(len(Y)):
        centX = centX + X[i]*dengrid[i,j]
        centY = centY + Y[j]*dengrid[i,j]
        total = total + dengrid[i,j]
print("Center X: ",centX/total," cm")
print("Center Y: ",centY/total," cm")
ybar = centY/total

exgrid = np.zeros((lenX, lenY))
eygrid = np.zeros((lenX, lenY))
phigrid = np.zeros((lenX, lenY))
length = xwindow * 1e5
for m in range(len(X)):
    if m%10==0: 
        print(m/lenX*100,"%")
    for n in range(len(Y)):
        for i in range(len(X)):
            for j in range(len(Y)):
                if not (i==m and j==n):
                    dx = X[m] - X[i]
                    dy = Y[n] - Y[j]
                    d=np.sqrt(np.square(dx)+np.square(dy))
                    exgrid[m,n] = exgrid[m,n] + 2*e*deltax*deltay*dengrid[i,j]*dx/d**2
                    eygrid[m,n] = eygrid[m,n] + 2*e*deltax*deltay*dengrid[i,j]*dy/d**2
                    phigrid[m,n] = phigrid[m,n] + 2*e*deltax*deltay*dengrid[i,j]*np.log(length/d)

phigrid = phigrid - phigrid[int(lenX/2),int(lenY/2)]

plt.title("Ex")
plt.imshow(np.transpose(exgrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
#plt.clim(2.6e16,3.4e16)
plt.show()

plt.title("Ey")
plt.imshow(np.transpose(eygrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
#plt.clim(2.6e16,3.4e16)
plt.show()

plt.title("Phi")
plt.imshow(np.transpose(phigrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
#plt.clim(2.6e16,3.4e16)
plt.show()

Exvx_calc = exgrid[:,int(lenY/2)]

plt.title("E_x along x")
plt.xlabel("x [cm]")
plt.ylabel("E [cgs]")
plt.plot(X,Exvx_calc, label="Calculation")
plt.legend(); plt.grid(); plt.show()

Exvy_calc = eygrid[int(lenX/2),:]

plt.title("E_y along y")
plt.xlabel("y [cm]")
plt.ylabel("E [cgs]")
plt.plot(Y,Exvy_calc, label="Calculation")
plt.legend(); plt.grid(); plt.show()

plt.title("E_y along y")
#plt.xlabel("y [cm]")
#plt.ylabel("E [cgs]")
plt.xlabel("y [um]")
plt.ylabel("E [V/m]")
#plt.plot(Y,Exvy_calc, label="Calculation")
plt.plot(Y*1e4,Exvy_calc*2.998e4, label="Calculation")
#plt.ylim([Exvy_calc[int(lenY/2)]*1.5,0])
plt.ylim([-1.5e9,-1e9])
plt.legend(); plt.grid(); plt.show()

Eyvx_calc = eygrid[:,int(lenY/2)]

plt.title("E_y along x")
#plt.xlabel("x [cm]")
#plt.ylabel("E [cgs]")
plt.xlabel("x [um]")
plt.ylabel("E [V/m]")
plt.plot(Y*1e4,Eyvx_calc*2.998e4, label="Calculation")
plt.legend(); plt.grid(); plt.show()

print("Offset [SI] = ",Eyvx_calc[int(len(Eyvx_calc)/2)]*2.998e4)