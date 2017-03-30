#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:57:12 2017

@author: litos
"""

import numpy as np
import matplotlib.pyplot as plt


def draw_ellipse(eps,T,T_mat=0,plot_on=True,npts=101):
    """ Draws beam phase space ellipse.
    
        Parameters
        ----------
        eps
            beam emittance in mm-mrad
        T
            Twiss parameter array:
                T[0] = beta (m)
                T[1] = alpha
                T[2] = gamma (1/m)
        T_mat
            Matched Twiss parameter array:
                T[0] = beta_mat (m)
                T[1] = alpha_mat
                T[2] = gamma_mat (1/m)
        npts
            number of points to draw

        Returns
        -------
        x
            x points drawn
        xp
            x' points drawn
    """
    [beta,alpha,gamma] = T
    sig   = np.sqrt(eps*beta)
    phi   = np.linspace(0,np.pi,npts)
    r_max = sig
    x     = r_max*np.cos(phi)

    # regular phase space ellipse
    # note use of abs to avoide edge cases that might go negative
    xp_top = (x>=-sig)*(x<=sig)*(1./beta)*\
    (-alpha*x+np.sqrt(np.abs(eps*beta-x**2)))
    xp_bot = (x>=-sig)*(x<=sig)*(1./beta)*\
    (-alpha*x-np.sqrt(np.abs(eps*beta-x**2)))
    xp     = np.append(xp_top,xp_bot)
    
    # normalized coordinates
    if T_mat == 0:
        T_mat = T            
    [beta_mat,alpha_mat,gamma_mat] = T_mat
    xpn_top = alpha_mat*x+beta_mat*xp_top
    xpn_bot = alpha_mat*x+beta_mat*xp_bot
    xpn     = np.append(xpn_top,xpn_bot)
    
    ## plot results
    if (plot_on):
        plt.plot(x,xp_top,'b.',x,xp_bot,'b.')
        plt.xlabel(r'$x$ (m)')
        plt.ylabel(r'$x^{\prime}$')
        plt.show()

    x = np.append(x,x)

    return [x,xp,xpn]
