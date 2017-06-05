#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:05:07 2017

@author: litos
"""

import numpy as np
from calc_rms import calc_rms
from frac_filter import frac_filter

def calc_twiss(parts,frac=1.0):
    
    npart = len(parts)
    x  = np.zeros(npart)
    xp = np.zeros(npart)
    y  = np.zeros(npart)
    yp = np.zeros(npart)
    z  = np.zeros(npart)
    gb = np.zeros(npart)
    for i in range(0,npart):
        [x[i],xp[i],y[i],yp[i],z[i],gb[i]] = parts[i][:]
        
    # calculate raw values
    rms_x   = calc_rms(x)
    rms_xp  = calc_rms(xp)
    rms_xxp = calc_rms(x*xp)
    rms_y   = calc_rms(y)
    rms_yp  = calc_rms(yp)
    rms_yyp = calc_rms(y*yp)
    rms_z   = calc_rms(z)
    rms_gb  = calc_rms(gb)

    avg_gb  = np.mean(gb)

    avg_x2  = np.mean(x**2)
    avg_xp2 = np.mean(xp**2)
    avg_xxp = np.mean(x*xp)
  
    rms_x_eps   = avg_gb*np.sqrt(avg_x2*avg_xp2-avg_xxp**2)
    rms_x_beta  = avg_gb*(rms_x**2)/rms_x_eps
    rms_x_gamma = avg_gb*(rms_xp**2)/rms_x_eps
    rms_x_alpha = -avg_gb*rms_xxp/rms_x_eps
    
    avg_y2  = np.mean(y**2)
    avg_yp2 = np.mean(yp**2)
    avg_yyp = np.mean(y*yp)
    
    rms_y_eps   = avg_gb*np.sqrt(avg_y2*avg_yp2-avg_yyp**2)
    rms_y_beta  = avg_gb*(rms_y**2)/rms_y_eps
    rms_y_gamma = avg_gb*(rms_yp**2)/rms_y_eps
    rms_y_alpha = -avg_gb*rms_yyp/rms_y_eps
    
    # filter particles based on emittance
    # and recalculate values
    if frac<1.0 :

        xpn = rms_x_alpha*x+rms_x_beta*xp
        x_eps2 = x**2 + xpn**2
        ix = frac_filter(x_eps2,frac)
        
        ypn = rms_y_alpha*y+rms_y_beta*yp
        y_eps2 = y**2 + ypn**2
        iy = frac_filter(y_eps2,frac)
        
        ixy = np.intersect1d(ix,iy)
        
        x  = x[ixy]
        xp = xp[ixy]
        y  = y[ixy]
        yp = yp[ixy]
        z  = z[ixy]
        gb = gb[ixy]
    
        rms_x   = calc_rms(x)
        rms_xp  = calc_rms(xp)
        rms_xxp = calc_rms(x*xp)
        rms_y   = calc_rms(y)
        rms_yp  = calc_rms(yp)
        rms_yyp = calc_rms(y*yp)
        rms_z   = calc_rms(z)
        rms_gb  = calc_rms(gb)
    
        avg_gb  = np.mean(gb)
    
        avg_x2  = np.mean(x**2)
        avg_xp2 = np.mean(xp**2)
        avg_xxp = np.mean(x*xp)
      
        rms_x_eps   = avg_gb*np.sqrt(avg_x2*avg_xp2-avg_xxp**2)
        rms_x_beta  = avg_gb*(rms_x**2)/rms_x_eps
        rms_x_gamma = avg_gb*(rms_xp**2)/rms_x_eps
        rms_x_alpha = -avg_gb*rms_xxp/rms_x_eps
        
        avg_y2  = np.mean(y**2)
        avg_yp2 = np.mean(yp**2)
        avg_yyp = np.mean(y*yp)
        
        rms_y_eps   = avg_gb*np.sqrt(avg_y2*avg_yp2-avg_yyp**2)
        rms_y_beta  = avg_gb*(rms_y**2)/rms_y_eps
        rms_y_gamma = avg_gb*(rms_yp**2)/rms_y_eps
        rms_y_alpha = -avg_gb*rms_yyp/rms_y_eps
    
    # output
    rms_parts = [rms_x,rms_xp,rms_y,rms_yp,rms_z,rms_gb]
    rms_X = [rms_x_eps,rms_x_beta,rms_x_alpha,rms_x_gamma]
    rms_Y = [rms_y_eps,rms_y_beta,rms_y_alpha,rms_y_gamma]
    
    return [rms_parts,rms_X,rms_Y]