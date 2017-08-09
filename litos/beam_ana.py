#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:20:28 2017

@author: mike
"""

import numpy as np
from collections import defaultdict
import mvee as mvee
import scipy.spatial as spatial

def calc_Bmag(Tb,Tl):
    """ Calculates Bmag.
    
        Bmag measures matching of a beam to a lattice.
    
        Parameters
        ----------
        Tb
            Twiss parameter array for beam:
                Tb[0] = beta (m)
                Tb[1] = alpha
                Tb[2] = gamma (1/m)
        Tl
            Twiss parameters for lattice:
                Tl[0] = beta (m)
                Tl[1] = alpha
                Tl[2] = gamma (1/m)

        Returns
        -------
        Bmag
            matching parameter >= 1
                =1 -> matched beam
                >1 -> mismatched beam
    """
    [beta_b,alpha_b,gamma_b] = Tb
    [beta_l,alpha_l,gamma_l] = Tl
    
    Bmag = (1/2)*(beta_b*gamma_l-2*alpha_b*alpha_l+gamma_b*beta_l)

    return Bmag

def calc_M(Tb,Tl):
    """ Calculates mismatch parameter M.
    
        Bmag measures matching of a beam to a lattice.
    
        Parameters
        ----------
        Tb
            Twiss parameter array for beam:
                Tb[0] = beta (m)
                Tb[1] = alpha
                Tb[2] = gamma (1/m)
        Tl
            Twiss parameters for lattice:
                Tl[0] = beta (m)
                Tl[1] = alpha
                Tl[2] = gamma (1/m)

        Returns
        -------
        Bmag
            matching parameter >= 1
                =1 -> matched beam
                >1 -> mismatched beam
    """
    Bmag = calc_Bmag(Tb,Tl)    
    M = Bmag + np.sqrt(Bmag**2 - 1)

    return M

def real2norm_coords(x,xp,beta,alpha):
    u = x/np.sqrt(beta)
    v = xp*np.sqrt(beta)+x*(alpha/beta)
    return [u,v]

def norm2real_coords(u,v,beta,alpha):
    x  = u*np.sqrt(beta)
    xp = (v-(alpha/beta)*x)/np.sqrt(beta)
    return [x,xp]

def calc_95ellipse(ebeam,step,tol=0.01):
    lips = defaultdict(dict)

    [u,v] = real2norm_coords(ebeam[step]["x"],ebeam[step]["xp"],\
                             ebeam[step]["beta"],ebeam[step]["alpha"])
    
    P = np.vstack((u,v)).T

    ET = mvee.EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(P,tol)

    lips["center"] = center
    lips["radii"]  = radii
    lips["rot"]    = rotation
    
    return lips
    