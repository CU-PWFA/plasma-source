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

def calc_95ellipse(u,v,hires=True,tol=0.01):
    lips = defaultdict(dict)

    # define points
    P = np.vstack((u,v)).T
    npoint = len(P)

    # use only points of convex hull for speed
    hull = spatial.ConvexHull(P)
    P = hull.points
    P = np.unique(P,axis=0)

    # find minimum volume enclosing ellipsoid for 100% of points
    ET = mvee.EllipsoidTool()
    (cent, rad, rot) = ET.getMinVolEllipse(P,tol)

    # define ellipse foci
    focmag = np.sqrt((rad[0]/2)**2+(rad[1]/2)**2)
    foc1 = np.dot([0,focmag],rot)+cent
    foc2 = -foc1
    
    # calculate string length of points wrt ellipse foci
    s = np.sqrt((u-foc1[0])**2+(v-foc1[1])**2) +\
        np.sqrt((u-foc2[0])**2+(v-foc2[1])**2)

    # reduce ellipse to only include 95% of particles
    if (hires):
        # cut a few points at a time until 5% are outside of ellipse
        i_s = np.argsort(s)
        ifrac = 1
        while len(i_s)>=0.95*npoint:
            # subtract 0.05*npart*(1/2)**ifrac particles
            sub = round(0.05*npoint*(0.5**ifrac))
            if sub==0:
                break
            i_s = i_s[:-sub]
            u = u[i_s]
            v = v[i_s]
    
            # define points
            P = np.vstack((u,v)).T
        
            # use only points of convex hull for speed
            hull = spatial.ConvexHull(P)
            P = hull.points
            P = np.unique(P,axis=0)
        
            # find minimum volume enclosing ellipsoid for 100% of points
            ET = mvee.EllipsoidTool()
            (cent, rad, rot) = ET.getMinVolEllipse(P,tol)
    
            # define ellipse foci
            focmag = np.sqrt((rad[0]/2)**2+(rad[1]/2)**2)
            foc1 = np.dot([0,focmag],rot)
            foc2 = -foc1
            
            # calculate string length of points wrt ellipse foci
            s = np.sqrt((u-foc1[0])**2+(v-foc1[0])**2) +\
                np.sqrt((u-foc2[0])**2+(v-foc2[0])**2)
            i_s = np.argsort(s)
    
            ifrac += 1
            continue
    else:
        # cut the outer 5% of points in one go
        i_s = np.argsort(s)
        i_s = i_s[:-round(0.05*len(i_s))]

    u = u[i_s]
    v = v[i_s]

    # define points
    P = np.vstack((u,v)).T

    # use only points of convex hull for speed
    hull = spatial.ConvexHull(P)
    P = hull.points
    P = np.unique(P,axis=0)

    # find minimum volume enclosing ellipsoid for 100% of points
    ET = mvee.EllipsoidTool()
    (cent, rad, rot) = ET.getMinVolEllipse(P,tol)
    
    # define ellipse foci
    focmag = np.sqrt((rad[0]/2)**2+(rad[1]/2)**2)
    foc1 = np.dot([0,focmag],rot)
    foc2 = -foc1

    lips["center"] = cent
    lips["radii"]  = rad
    lips["rot"]    = rot
    lips["foci"]   = [foc1,foc2]
    lips["ecc"]    = np.sqrt(1-(max(rad)/min(rad))**2)
    lips["area"]   = np.pi*rad[0]*rad[1]

    return lips

#def calc_95ellipseB(u,v,beta,alpha,hires=True,tol=0.01):
#    lips = defaultdict(dict)
#
#    # define points
#    P = np.vstack((u,v)).T
#    npoint = len(P)
#
#    # use only points of convex hull for speed
#    hull = spatial.ConvexHull(P)
#    P = hull.points
#    P = np.unique(P,axis=0)
#
#    # angle of major axis
#    gamma = (1+alpha**2)/beta
#    theta = 0.5*np.arctan2(2*alpha,gamma-beta)
#    
#    # rotate points
#    rot = [[ np.cos(theta), np.sin(theta)],\
#           [-np.sin(theta), np.cos(theta)]]
#    P = np.dot(rot,P)
#
#    # first guess for foci: half of max
#    focd = max(abs(P[:,0]))
#    foc1 = [-focd/2,0]
#    foc2 = [+focd/2,0]
#    
#    # calculate string length of points wrt ellipse foci
#    s = np.sqrt((u-foc1[0])**2+(v-foc1[1])**2) +\
#        np.sqrt((u-foc2[0])**2+(v-foc2[1])**2)
#    
#    
#    
#    args = [gb0,eps0,beta0,alpha0,gamma0,\
#        np0,shape,z,Lp0,\
#        targ_beta,targ_alpha,targ_gamma]
#    
#    x0 = [focd,ecc]
#    ## perform fit
#    fitres = minimize(max_s,x0,(args,),\
#                      bounds=bnds,\
#                      method='L-BFGS-B',\
#                      options={'disp':False})
##                  constraints=cons,\
#    
#    
#    
#    
#
#    lips["center"] = cent
#    lips["radii"]  = rad
#    lips["rot"]    = rot
#    lips["foci"]   = [foc1,foc2]
#    lips["ecc"]    = np.sqrt(1-(max(rad)/min(rad))**2)
#    lips["area"]   = np.pi*rad[0]*rad[1]
#
#    return lips
#
#def find_max_s(x,args):
#    [focd,ecc] = x
#    [u,v] = args
#    foc1 = [-focd/2,0]
#    foc2 = [+focd/2,0]
#    s = np.sqrt((u-foc1[0])**2+(v-foc1[1])**2) +\
#        np.sqrt((u-foc2[0])**2+(v-foc2[1])**2)
#    return max(s)