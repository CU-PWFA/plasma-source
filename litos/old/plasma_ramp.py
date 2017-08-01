#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:06:38 2017

@author: litos
"""

import numpy as np

def plasma_ramp(npl0,shape,s,args,updown='up',dgds0=0):
    """ Creates plasma density profile as a func. of s.
    
        Parameters
        ----------
        npl0
            cm^-3; flat-top plasma density
        shape
            string; ramp shape
        s
            m; longitudinal coordinates
        args
            shape-dependent arguments
            L: m; full length of ramp
            hw: m; half-width, half-max of ramp
            P: exponential parameter (not used for all shapes)
        updown
            string; whether to create up or down ramp
        dgds0
            m^-1; energy gain per unit length at flat-top density

        Returns
        -------
        npl
            cm^-3; plasma density as a func. of s
    """
    
    if shape.lower() == 'gauss': # normal gauss. func.
        [L,hw] = args[:2]
        sig = hw/(np.sqrt(2*np.log(2)))
        if updown.lower() == 'up':
            npl = npl0*((s<L)*np.exp(-((s-L)**2)/(2*(sig**2))) +\
                       (s>=L)*1)
        else:
            npl = npl0*((s<L)*np.exp(-(s**2)/(2*(sig**2))) +\
                       (s>=L)*0)
    elif shape.lower() == 'sigmoid': # sigmoidal func.
        [L,hw] = args[:2]
        sig = hw
        if updown.lower() == 'up':
            #        npl = npl0*(1/(1+np.exp(-(s-(L-8*sig))/sig)))
#            npl = npl0*(1/(1+np.exp(-(s-(L-2*sig))/(sig/4))))
            npl = npl0*(1+np.exp(np.log(1e-3)*(s-(L-sig))/sig))
        else:
#            npl = npl0*(1/(1+np.exp(+(s-(L+2*sig))/(sig/4))))
            npl = npl0*(1+np.exp(np.log(1e-3)*((L-sig)-s)/sig))
    elif shape.lower() == 'gen_gauss': # generalized gauss. func.
        [L,hw,P] = args[:3]
        sig = hw/(np.sqrt(2*np.log(2)))
        npl = npl0*((s<L)*np.exp(-((s-L)**P)/(2*(sig**P))) +\
                   (s>=L)*1)
    elif shape.lower() == 'trap': # trapezoidal func.
        [L,hw] = args[:2]
        npl = npl0*((s<=L-2*hw)*0 +\
                   (s>L-2*hw)*(s<L)*(s-(L-2*hw))/(2*hw) +\
                   (s>=L)*1)
    elif shape.lower() == 'genlog': # generalized logistical fun.
        [L,hw,P] = args[:3]
        sig = hw
        npl = npl0*(1/((1+np.exp(-(s-(L-8*sig))/sig))**P))
    elif shape.lower() == 'gomp': # gompertz func.
        [L,hw] = args[:2]
        sig = hw
        if updown.lower() == 'up':
            npl = npl0*((s<L)*np.exp(1-0.99(s-L)/sig -\
                             np.exp(-0.99*(s-L)/sig)) +\
                       (s>=L)*1)
        else:
            npl = npl0*((s>L)*np.exp(1+0.99(s-L)/sig -\
                             np.exp(+0.99*(s-L)/sig)) +\
                       (s<=L)*1)
    elif shape.lower() == 'ngomp': # negative-gompertz func.
        [L,hw] = args[:2]
        sig = hw
        if updown.lower() == 'up':
            npl = npl0*(np.exp(1+1.46(s-L)/sig-np.exp(+1.46*(s-L)/sig)))
        else:
            npl = npl0*(np.exp(1-1.46(s-L)/sig-np.exp(-1.46*(s-L)/sig)))       
    elif shape.lower() == 'xu3': # Xu PRL 2016 ramp shapes
        [L,hw] =  args[:2]
        sig = 2*hw
        if updown.lower() == 'up':
            npl = npl0*((s<L)*(1/(1-(s-L)/sig)) +\
                       (s>=L)*1)
        else:
            npl = npl0*((s>L)*(1/(1+(s-L)/sig)) +\
                       (s<=L)*1)
    elif shape.lower() == 'xu4': # Xu PRL 2016 ramp shapes
        [L,sig] =  args[:2]
        sig = hw/(np.sqrt(2)-1)
        if updown.lower() == 'up':
            npl = npl0*((s<L)*(1/((1-(s-L)/sig)**2)) +\
                       (s>=L)*1)
        else:
            npl = npl0*((s>L)*(1/((1+(s-L)/sig)**2)) +\
                       (s<=L)*1)
    elif shape.lower() == 'xu5': # Xu PRL 2016 ramp shapes
        [L,hw] =  args[:2]
        sig = hw/(2**(1/4)-1)
        if updown.lower() == 'up':
            npl = npl0*((s<L)*(1/((1-(s-L)/(2*sig))**4)) +\
                       (s>=L)*1)
        else:
            npl = npl0*((s>L)*(1/((1+(s-L)/(2*sig))**4)) +\
                       (s<=L)*1)
    else:
        print('bad plasma ramp shape' ) 

    dgds = dgds0*np.sqrt(npl/npl0)*(2*np.sqrt(npl/npl0)-1) # wake strength ~sqrt(np), phase ~sqrt(np)

    return [npl, dgds]