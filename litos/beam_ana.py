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
import matplotlib.pyplot as plt

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

def calc_frac_ellipse(u,v,frac=0.95,hires=False,tol=0.01):
    lips = defaultdict(dict)

    # define points
    npoint = len(u)
    P = np.vstack((u,v)).T

    # use only points of convex hull for speed
    # NOTE: *could* be useful, but need to figure out
    # how to make it select fewer points!
#    hull = spatial.ConvexHull(P,qhull_options='QbB')
#    P = hull.points
#    P = np.unique(P,axis=0)
#    print(len(hull.points))
#    print(len(P))

    # find minimum volume enclosing ellipsoid for 100% of points
    ET = mvee.EllipsoidTool()
    (cent, rad, rot) = ET.getMinVolEllipse(P,tol)

    # center points
    Pc = np.zeros([len(P),2])
    Pc[:,0] = P[:,0] - cent[0]
    Pc[:,1] = P[:,1] - cent[1]

    # rotate points
    R = np.rot90(rot)
    Pr = np.zeros([len(Pc),2])
    for i in range(len(Pc)):
        Pr[i,:] = np.dot(R,Pc[i,:].T)

#    # plot to check
#    figA, axA1 = plt.subplots(1,1,sharey=True)
#    axA1.scatter(u,v,c='r',s=0.5)
#    axA1.scatter(Pr[:,0],Pr[:,1],c='b',s=0.5)
#    axA1.set_xlim([-1.1*max(abs(u)),+1.1*max(abs(u))])
#    axA1.set_ylim([-1.1*max(abs(v)),+1.1*max(abs(v))])
#    plt.show()

    # normalize points
    Pn = np.zeros([len(Pr),2])
    Pn[:,0] = Pr[:,0]/max(abs(Pr[:,0]))
    Pn[:,1] = Pr[:,1]/max(abs(Pr[:,1]))
    
    # calculate radius^2 of all points
    Cn = [np.mean(Pn[:,0]),np.mean(Pn[:,1])]
    r2 = (Pn[:,0]-Cn[0])**2 + (Pn[:,1]-Cn[1])**2
    Pt = Pn
    
    # reduce ellipse to only include fraction of particles
    if (hires):
        # cut a few points at a time until (1-frac) are outside of ellipse
        # NOTE: not obvious that this method is measurably better
        # than the other method
        i_r2 = np.argsort(r2)
        Pt = Pt[i_r2,:]
        ipow = 1
        while len(i_r2)>=frac*npoint:
            # subtract (1-frac)*npart*(1/2)**ipow particles
            sub = round((1-frac)*npoint*(0.5**ipow))
            if sub==0:
                break
            Pt = Pt[:-sub,:]
            # calculate radius^2 of all points
            Ct = [np.mean(Pt[:,0]),np.mean(Pt[:,1])]
            r2 = (Pt[:,0]-Ct[0])**2 + (Pt[:,1]-Ct[1])**2
            i_r2 = np.argsort(r2)
            Pt = Pt[i_r2,:]
            ipow += 1
            continue
    else:
        # cut the outer fraction of points in one go
        i_r2 = np.argsort(r2)
        i_r2 = i_r2[:-round((1-frac)*len(i_r2))]
        Pt = Pn[i_r2,:]
    
#    # plot to check
#    figB, axB1 = plt.subplots(1,1,sharey=True)
##    axB1.scatter(u,v,c='r',s=0.5)
#    axB1.scatter(Pn[:,0],Pn[:,1],c='b',s=0.5)
#    axB1.scatter(Pt[:,0],Pt[:,1],c='r',s=0.5)
##    axB1.set_xlim([-1.1*max(abs(u)),+1.1*max(abs(u))])
##    axB1.set_ylim([-1.1*max(abs(v)),+1.1*max(abs(v))])
#    plt.show()
    
    # un-normalize points
    Pt[:,0] = Pt[:,0]*max(abs(Pr[:,0]))
    Pt[:,1] = Pt[:,1]*max(abs(Pr[:,1]))
    # rotate back to normal orientation
    for i in range(len(Pt)):
        Pt[i,:] = np.dot(np.linalg.inv(R),Pt[i,:])
    # un-center points
    Pt[:,0] = Pt[:,0] + cent[0]
    Pt[:,1] = Pt[:,1] + cent[1]
    
#    # plot
#    axA1.scatter(Pt[:,0],Pt[:,1],c='g',s=0.5)
#    plt.show()
    
#    # plot
#    figC, axC1 = plt.subplots(1,1,sharey=True)
#    axC1.scatter(u,v,c='r',s=1.0)
#    axC1.scatter(Pt[:,0],Pt[:,1],c='g',s=0.1)
#    axC1.set_xlim([-1.1*max(abs(u)),+1.1*max(abs(u))])
#    axC1.set_ylim([-1.1*max(abs(v)),+1.1*max(abs(v))])
#    plt.show()

    # use only points of convex hull for speed
    # NOTE: *could* be useful, but need to figure out
    # how to make it select fewer points!
#    hull  = spatial.ConvexHull(Pt)
#    Pt = hull.points
#    Pt = np.unique(Pt,axis=0)

    # find minimum volume enclosing ellipsoid for remaining points
    ET = mvee.EllipsoidTool()
    (cent, rad, rot) = ET.getMinVolEllipse(Pt,tol)

    # define ellipse foci
    R    = np.rot90(rot)
    focd = 2*np.sqrt((rad[0]/2)**2+(rad[1]/2)**2)
    foc1 = np.dot(np.linalg.inv(R),[focd/2,0])
    foc2 = -foc1

    lips["center"] = cent
    lips["radii"]  = rad
    lips["rot"]    = rot
    lips["focd"]   = focd
    lips["foci"]   = [foc1,foc2]
    lips["ecc"]    = np.sqrt(1-(min(rad)/max(rad))**2)
    lips["area"]   = np.pi*rad[0]*rad[1]
    lips["Pt"]     = Pt

    return lips
