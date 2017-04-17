#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:06:38 2017

@author: litos
"""

import numpy as np

def plasma_ramp(np0,shape,z,args,updown='up'):
    """ Creates plasma density profile as a func. of z.
    
        Parameters
        ----------
        np0
            cm^-3; flat-top plasma density
        shape
            string; ramp shape
        z
            m; longitudinal coordinates
        args
            shape-dependent arguments
            L: m; full length of ramp
            sig: m; characteristic falloff length of ramp
            P: exponential parameter (not used for all shapes)
        updown
            string; whether to create up or down ramp

        Returns
        -------
        npl
            cm^-3; plasma density as a func. of z
    """
    
    if shape.lower() == 'gauss': # normal gauss. func.
        [L,sig] = args[0:2]
        if updown.lower() == 'up':
            npl = np0*((z<L)*np.exp(-((z-L)**2)/(2*(sig**2))) +\
                       (z>=L)*1)
        else:
            npl = np0*((z>L)*np.exp(-((z-L)**2)/(2*(sig**2))) +\
                       (z<=L)*1)
    elif shape.lower() == 'gen_gauss': # generalized gauss. func.
        [L,sig,P] = args[0:3]
        npl = np0*((z<L)*np.exp(-((z-L)**P)/(2*(sig**P))) +\
                   (z>=L)*1)
    elif shape.lower() == 'trap': # trapezoidal func.
        [L,sig] = args[0:2]
        npl = np0*((z<=L-2*sig)*0 +\
                   (z>L-2*sig)*(z<L)*(z-(L-2*sig))/(2*sig) +\
                   (z>=L)*1)
    elif shape.lower() == 'sigmoid': # sigmoidal func.
        [L,sig] = args[0:2]
        if updown.lower() == 'up':
            #        npl = np0*(1/(1+np.exp(-(z-(L-8*sig))/sig)))
            npl = np0*(1/(1+np.exp(-(z-(L-2*sig))/(sig/4))))
        else:
            npl = np0*(1/(1+np.exp(+(z-(L+2*sig))/(sig/4))))
    elif shape.lower() == 'genlog': # generalized logistical fun.
        [L,sig,P] = args[0:3]
        npl = np0*(1/((1+np.exp(-(z-(L-8*sig))/sig))**P))
    elif shape.lower() == 'gomp': # gompertz func.
        [L,sig] = args[0:2]
        npl = np0*(1-np.exp(-np.exp((z-L)/sig)))
    elif shape.lower() == 'xu3': # Xu PRL 2016 ramp shapes
        [L,sig] =  args[0:2]
        npl = np0*((z<L)*(1/(1-2*(z-L)/sig)) +\
                   (z>=L)*1)
    elif shape.lower() == 'xu4': # Xu PRL 2016 ramp shapes
        [L,sig] =  args[0:2]
        npl = np0*((z<L)*(1/((1-(z-L)/sig)**2)) +\
                   (z>=L)*1)
    elif shape.lower() == 'xu5': # Xu PRL 2016 ramp shapes
        [L,sig] =  args[0:2]
        npl = np0*((z<L)*(1/((1-(z-L)/(2*sig))**4)) +\
                   (z>=L)*1)
    else:
        print('bad plasma ramp shape' ) 

    return npl