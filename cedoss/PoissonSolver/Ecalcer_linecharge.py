#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 12:00:19 2019

Given a charge density distribution in 2D, simply sum up the electric field contributions
in all 3 dimensions for every point.

This version assumes E can be calculated from a sueprposition of line charges.

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

e = 4.8032e-10
pi = np.pi
"""
n_cent = 3.0e16  #cm-3
radius = 250e-6 *1e2 #cm
slope  = -9.91e16 #cm-4
xwindow = ywindow = 550e-6 * 1e2  #cm
deltax = deltay = 10e-6 * 1e2  #cm
"""
n_cent = 3.0e16  #cm-3
radius = 43.55e-6 *1e2 #cm
slope  = -9.91e16 #cm-4
xwindow = ywindow = 100e-6 * 1e2  #cm
deltax = deltay = 1e-6 * 1e2  #cm

X = np.arange(-1/2*xwindow, 1/2*xwindow, deltax)
Y = np.arange(-1/2*ywindow, 1/2*ywindow, deltay)

lenX = len(X); lenY = len(Y)

dengrid = np.zeros((lenX, lenY))
for i in range(len(X)):
    for j in range(len(Y)):
        if np.sqrt(np.square(X[i])+np.square(Y[j])) < radius:
            dengrid[i][j] = n_cent + Y[j]*slope
        else:
            dengrid[i][j] = 0
            
plt.title("Density")
plt.imshow(np.transpose(dengrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.clim(n_cent+slope*radius,n_cent-slope*radius)#2.6e16,3.4e16)
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
Exvx_theory = 2*pi*e*n_cent*X

plt.title("E_x along x")
plt.xlabel("x [cm]")
plt.ylabel("E [cgs]")
plt.plot(X,Exvx_calc, label="Calculation")
plt.plot(X,Exvx_theory,label="Theory")
plt.ylim([min(Exvx_calc)*1.2,max(Exvx_calc)*1.2])
plt.legend(); plt.grid(); plt.show()

Exvy_calc = eygrid[int(lenX/2),:]
Exvy_theory = 2*pi*e*n_cent*(Y-2*ybar)
Exvy_theory4 = 2*pi*e*n_cent*(Y-2*ybar) + 1.5*pi*e*slope*Y**2

plt.title("E_y along y")
plt.xlabel("y [cm]")
plt.ylabel("E [cgs]")
plt.plot(Y,Exvy_calc, label="Calculation")
plt.plot(Y,Exvy_theory,label="Theory")
plt.plot(Y,Exvy_theory4,label="Theory4")
plt.ylim([min(Exvy_calc)*1.2,max(Exvy_calc)*1.2])
plt.legend(); plt.grid(); plt.show()

plt.plot(Y[5:35],(Exvy_theory4[5:35] - Exvy_calc[5:35]))
plt.show()

phivx_calc = phigrid[:,int(lenY/2)]
phivx_theory = -pi*e*n_cent*np.square(X)

plt.title("Electric potential along x")
plt.xlabel("x [cm]")
plt.ylabel("V [cgs]")
plt.plot(X,phivx_calc,label="Calculation")
plt.plot(X,phivx_theory,label="Theory")
plt.legend(); plt.grid(); plt.show()

phivy_calc = phigrid[int(lenX/2),:]
phivy_theory = -pi*e*n_cent*np.square(Y-2*ybar)
phivy_theory4 = -pi*e*n_cent*np.square(Y-2*ybar)-1/2*pi*slope*Y**3*e

plt.title("Electric potential along y")
plt.xlabel("y [cm]")
plt.ylabel("V [cgs]")
plt.plot(Y,phivy_calc,label="Calculation")
plt.plot(Y,phivy_theory+max(phivy_calc),label="Theory")
plt.plot(Y,phivy_theory4+max(phivy_calc),label="Theory4")
plt.legend(); plt.grid(); plt.show()

phigrid_theory = np.zeros((lenX, lenY))
for i in range(lenX):
    for j in range(lenY):
        if np.sqrt(np.square(X[i])+np.square(Y[j])) < radius:
            phigrid_theory[i,j] = -pi*e*n_cent*(np.square(X[i])+np.square(Y[j]-2*ybar))-1/2*pi*slope*e*Y[j]**3 - 1/2*pi*slope*e*X[i]**2*Y[j]
        else:
            phigrid_theory[i,j] = phigrid[i,j]
        
plt.title("Phi_theory")
plt.imshow(np.transpose(phigrid_theory), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
#plt.clim(2.6e16,3.4e16)
plt.show()

plt.title("Phi_theory - Phi")
plt.imshow(np.transpose(phigrid_theory-phigrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
#plt.clim(2.6e16,3.4e16)
plt.show()

eygrid_theory = np.zeros((lenX, lenY))
exgrid_theory = np.zeros((lenX, lenY))
for i in range(lenX):
    for j in range(lenY):
        if np.sqrt(np.square(X[i])+np.square(Y[j])) < radius:
            exgrid_theory[i,j] = 2*pi*e*n_cent*X[i] + pi*e*slope*X[i]*Y[j]
            eygrid_theory[i,j] = 2*pi*e*n_cent*(Y[j]-2*ybar)+3/2*pi*slope*e*Y[j]**2 + 1/2*pi*e*slope*X[i]**2
        else:
            exgrid_theory[i,j] = exgrid[i,j]
            eygrid_theory[i,j] = eygrid[i,j]
        
plt.title("Ey_theory")
plt.imshow(np.transpose(eygrid_theory), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.show()

plt.title("Ey_theory - Ey")
plt.imshow(np.transpose(eygrid_theory-eygrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.show()

plt.title("Ex_theory")
plt.imshow(np.transpose(exgrid_theory), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.show()

plt.title("Ex_theory - Ex")
plt.imshow(np.transpose(exgrid_theory-exgrid), extent=(X[0],X[-1],Y[0],Y[-1]), origin = 'lower')
plt.colorbar()
plt.show()