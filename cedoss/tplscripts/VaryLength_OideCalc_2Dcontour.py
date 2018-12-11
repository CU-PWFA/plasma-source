#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:53:54 2018

Taking into account asymmetric beam profile, caclulates parameters
for maximum luminosity over a 2D contour of initial beta & density, and
each case is optomized over thickness.

@author: chris
"""

import sys
sys.path.insert(0, "../")
from modules import TPLFocalLength as Foc
from modules import OideCalc as Oide
from modules import CalcEmitGrowth as W2
import numpy as np
import matplotlib.pyplot as plt
from beamprop_v2 import BeamPropFuncs as PProp

func =1 

#ILC
emity = 35e-9 *100 #cm-rad
emitx = 10e-6 *100 #cm-rad
gam = Foc.gam_def*50
sigmaE = 0.00124
delta = np.sqrt(3)*sigmaE
ymode = 1 #1 for nm, 0 for um
lumi=1

emit = emitx

d_set = 0
len_max = 400 #um

beta_i_arr = np.linspace(0.05, 1, 30) #cm
n0_arr = np.linspace(1e20, 1e21, 30) #cm^-3
sigmin_image = np.zeros((len(beta_i_arr),len(n0_arr)))
optlen_image = np.zeros((len(beta_i_arr),len(n0_arr)))
sigx_image = np.zeros((len(beta_i_arr),len(n0_arr)))
sigy_image = np.zeros((len(beta_i_arr),len(n0_arr)))

for x in range(len(beta_i_arr)):
    if x%5 == 0: print(x/len(beta_i_arr)*100,"%");
    for y in range(len(n0_arr)):
        
        beta_i = beta_i_arr[x]
        n0 = n0_arr[y]
        
        len_arr = np.linspace(50,len_max,101)/1e4 #cm
        sig_product_arr = np.zeros(len(len_arr))
        sigx_arr = np.zeros(len(len_arr))
        sigy_arr = np.zeros(len(len_arr))
        
        for i in range(len(len_arr)):
            L = len_arr[i]
            K = Foc.Calc_K(n0, gam)
            focal = Foc.Calc_Focus_KLength(K, L)
            
            KLls_set = [K, L, Oide.Get_ls_thick(K,L,beta_i,d_set)]
            beta_f = Foc.Calc_ThickBetaStar_DeltaOff_UnNormalized(K,L,beta_i,d_set)
            
            F_val = Oide.F_Oide(KLls_set)
            
            if lumi==1:
                sigx_arr[i] = Oide.Calc_SigOide(F_val, emitx, gam, beta_f)
                sigy_arr[i] = Oide.Calc_SigOide(F_val, emity, gam, beta_f)
                sig_product_arr[i]=sigx_arr[i]*sigy_arr[i]
        sigmin_image[x][y] = np.min(sig_product_arr)
        optlen_image[x][y] = len_arr[np.argmin(sig_product_arr)]
        sigx_image[x][y] = sigx_arr[np.argmin(sig_product_arr)]
        sigy_image[x][y] = sigy_arr[np.argmin(sig_product_arr)]

if func == 0:
    minloc = PProp.PlotContour(sigmin_image, beta_i_arr, n0_arr, r'$\beta_i$ [cm]', r'$n_0\mathrm{\ cm^{-3}}$')

#def PlotContour(contour, x_arr, y_arr, x_label, y_label, simple = False):
if func == 1:
        # find location of min(B)
    x_arr = beta_i_arr
    y_arr = n0_arr
    x_label =r'$\beta_i$ [cm]'
    y_label =r'$n_0\mathrm{\ cm^{-3}}$'
    
    i_Bmin_x = np.argmin(np.min(sigmin_image,1))
    i_Bmin_y = np.argmin(np.min(sigmin_image,0))
    Bmin_x = x_arr[i_Bmin_x]
    Bmin_y = y_arr[i_Bmin_y]
    Bmin   = np.min(np.min(sigmin_image))
    
    print('Optimal beta [cm]  : ',Bmin_x)
    print('Optimal den [cm^-3]: ',Bmin_y)
    print('Sigx [nm]  : ',sigx_image[i_Bmin_x][i_Bmin_y]*1e7)
    print('Sigy [nm]  : ',sigy_image[i_Bmin_x][i_Bmin_y]*1e7)
    print('Length [um]: ',optlen_image[i_Bmin_x][i_Bmin_y]*1e4)
    
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,sigmin_image*1e14,100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$(\sigma_x\times\sigma_y)_{min}\mathrm{\ nm^2}$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.show()
    
    ##Now for optimal length    
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,optlen_image*1e4,100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$TPL_{L}\mathrm{\ \mu m}$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.show()

    ##Sigx
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,sigx_image*1e7,100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$\sigma_x\mathrm{\ nm}$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.show()
    
    ##Sigy
    X = np.tile(x_arr.reshape(-1,1),(1,len(y_arr)))
    Y = np.tile(y_arr.T,(len(x_arr),1))
    fig, axes = plt.subplots(1,1, sharey=True)
    plt.contourf(X,Y,sigy_image*1e7,100,\
                cmap=plt.get_cmap('Vega20c'),\
                linewidth=2.0)
    cbar = plt.colorbar()
    plt.scatter(Bmin_x,Bmin_y,color='k')
    cbar.ax.set_ylabel(r'$\sigma_y\mathrm{\ nm}$')
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.show()