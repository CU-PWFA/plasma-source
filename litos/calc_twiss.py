#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:05:07 2017

@author: litos
"""

import numpy as np
from calc_rms import calc_rms

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
    
    rms_x   = calc_rms(x,frac)
    rms_xp  = calc_rms(xp,frac)
    rms_xxp = calc_rms(x*xp,frac)
    rms_y   = calc_rms(y,frac)
    rms_yp  = calc_rms(yp,frac)
    rms_yyp = calc_rms(y*yp,frac)
    rms_z   = calc_rms(z,frac)
    rms_gb  = calc_rms(gb,frac)

    avg_gb  = np.mean(gb)

    avg_x2  = np.mean(x**2)
    avg_xp2 = np.mean(xp**2)
    avg_xxp = np.mean(x*xp)
    
    rms_x_eps   = avg_gb*np.sqrt(avg_x2*avg_xp2-avg_xxp**2)
    rms_x_beta  = (rms_x**2)/(rms_x_eps*avg_gb)
    rms_x_gamma = (rms_xp**2)/(rms_x_eps*avg_gb)
    rms_x_alpha = -rms_xxp/(rms_x_eps*avg_gb)
    
    avg_y2  = np.mean(y**2)
    avg_yp2 = np.mean(yp**2)
    avg_yyp = np.mean(y*yp)
    
    rms_y_eps   = avg_gb*np.sqrt(avg_x2*avg_xp2-avg_xxp**2)
    rms_y_beta  = (rms_y**2)/(rms_y_eps*avg_gb)
    rms_y_gamma = (rms_yp**2)/(rms_y_eps*avg_gb)
    rms_y_alpha = -rms_yyp/(rms_y_eps*avg_gb)
    
    rms_parts = [rms_x,rms_xp,rms_y,rms_yp,rms_z,rms_gb]
    rms_X = [rms_x_eps,rms_x_beta,rms_x_alpha,rms_x_gamma]
    rms_Y = [rms_y_eps,rms_y_beta,rms_y_alpha,rms_y_gamma]
    
    return [rms_parts,rms_X,rms_Y]