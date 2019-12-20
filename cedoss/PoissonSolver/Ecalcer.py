#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:56:04 2019

Given a charge density distribution in 2D, simply sum up the electric field contributions
in all 3 dimensions for every point.

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

e = 4.8032e-10
pi = np.pi

n_cent = 3.0e16 #cm-3
radius = 400e-6 *1e2 #cm
slope  = -9.91e16 #cm-4

xwindow = ywindow = 1e-3 * 1e2 * 1 #cm
deltax = deltay = 0.5e-6 * 1e2 * 50  #cm

zwindow = xwindow * 10
deltaz = deltax * 10

testfac = 0

X = np.arange(-1/2*xwindow, 1/2*xwindow, deltax)
Y = np.arange(-1/2*ywindow, 1/2*ywindow, deltay)
Z = np.arange(-1/2*zwindow, 1/2*zwindow, deltaz)

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
plt.clim(2.6e16,3.4e16)
plt.show()

exgrid = np.zeros((lenX, lenY))
eygrid = np.zeros((lenX, lenY))
phigrid = np.zeros((lenX, lenY))
for m in range(len(X)):
    if m%10==0: 
        print(m/lenX*100,"%")
    for n in range(len(Y)):
        for i in range(len(X)):
            for j in range(len(Y)):
                if not (i==m and j==n):
                    for k in range(len(Z)):
                        dx = X[m] - X[i]
                        dy = Y[n] - Y[j]
                        d=np.sqrt(np.square(dx)+np.square(dy)+np.square(Z[k]))
                        exgrid[m,n] = exgrid[m,n] + e*deltax*deltay*deltaz*dengrid[i,j]*dx/d**3
                        eygrid[m,n] = eygrid[m,n] + e*deltax*deltay*deltaz*dengrid[i,j]*dy/d**3
                        phigrid[m,n] = phigrid[m,n] + e*deltax*deltay*deltaz*dengrid[i,j]/d

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
Exvx_theory_10x = Exvx_theory*testfac

plt.title("E_x along x")
plt.xlabel("x [cm]")
plt.ylabel("E [cgs]")
plt.plot(X,Exvx_calc, label="Calculation")
plt.plot(X,Exvx_theory,label="Theory")
if testfac != 0:
    plt.plot(X,Exvx_theory_10x,label="Theory x " + str(testfac))
plt.legend(); plt.grid(); plt.show()

Exvy_calc = eygrid[int(lenX/2),:]
Exvy_theory = 2*pi*e*n_cent*X
Exvy_theory_10x = Exvy_theory*testfac

plt.title("E_y along y")
plt.xlabel("y [cm]")
plt.ylabel("E [cgs]")
plt.plot(Y,Exvy_calc, label="Calculation")
plt.plot(Y,Exvy_theory,label="Theory")
if testfac != 0:
    plt.plot(Y,Exvy_theory_10x,label="Theory x " + str(testfac))
plt.legend(); plt.grid(); plt.show()

phivx_calc = phigrid[:,int(lenY/2)]
phivx_theory = -pi*e*n_cent*np.square(X)
phivx_theory_10x = phivx_theory*testfac

plt.title("Electric potential along x")
plt.xlabel("x [cm]")
plt.ylabel("V [cgs]")
plt.plot(X,phivx_calc,label="Calculation")
plt.plot(X,phivx_theory,label="Theory")
if testfac != 0:
    plt.plot(X,phivx_theory_10x,label="Theory x " + str(testfac))
plt.legend(); plt.grid(); plt.show()

phivy_calc = phigrid[int(lenX/2),:]
phivy_theory = -pi*e*n_cent*np.square(Y)
phivy_theory_10x = phivy_theory*testfac

plt.title("Electric potential along y")
plt.xlabel("y [cm]")
plt.ylabel("V [cgs]")
plt.plot(Y,phivy_calc,label="Calculation")
plt.plot(Y,phivy_theory,label="Theory")
if testfac != 0:
    plt.plot(Y,phivy_theory_10x,label="Theory x " + str(testfac))
plt.legend(); plt.grid(); plt.show()