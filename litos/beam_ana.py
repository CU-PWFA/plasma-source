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
import mike_math as mm

def calc_ebeam_cent(ebeam,step=0,frac=1.00):
    ebeam_cent = defaultdict(dict)
    
    x  = ebeam[step]["x"] # m
    xp = ebeam[step]["xp"] # rad

    cent_x = mm.calc_mean(x,frac)
    cent_xp = mm.calc_mean(xp,frac)
    cent_xxp = mm.calc_mean(x*xp,frac)
    
    gb = ebeam[step]["gb"]
    avg_gb = np.mean(gb)
    
    ebeam_cent["x"]       = cent_x
    ebeam_cent["xp"]      = cent_xp
    ebeam_cent["xxp"]     = cent_xxp
    
    return ebeam_cent
  
def calc_ebeam_rms(ebeam,step=0,frac=1.00):
    ebeam_rms = defaultdict(dict)
    
    x  = ebeam[step]["x"] # m
    xp = ebeam[step]["xp"] # rad

    rms_x = mm.calc_rms(x,frac)
    rms_xp = mm.calc_rms(xp,frac)
    rms_xxp = mm.calc_rms(x*xp,frac)
    
    gb = ebeam[step]["gb"]
    avg_gb = np.mean(gb)
    
    avg_x2  = np.mean((x-np.mean(x))**2)
    avg_xp2 = np.mean((xp-np.mean(xp))**2)
    avg_xxp = np.mean((x-np.mean(x))*(xp-np.mean(xp)))

    rms_x_eps   = avg_gb*np.sqrt(avg_x2*avg_xp2-avg_xxp**2)
#    rms_x_eps   = avg_gb*np.sqrt((rms_x**2)*(rms_xp**2)-avg_xxp**2)
    rms_x_beta  = avg_gb*(avg_x2)/rms_x_eps
    rms_x_gamma = avg_gb*(avg_xp2)/rms_x_eps
    rms_x_alpha = -avg_gb*avg_xxp/rms_x_eps
    rms_x_phase = np.arctan2(2*rms_x_alpha,rms_x_gamma-rms_x_beta)/2
    
    ebeam_rms["x"]       = rms_x
    ebeam_rms["xp"]      = rms_xp
    ebeam_rms["xxp"]     = rms_xxp
    ebeam_rms["x_eps"]   = rms_x_eps
    ebeam_rms["x_beta"]  = rms_x_beta
    ebeam_rms["x_alpha"] = rms_x_alpha
    ebeam_rms["x_gamma"] = rms_x_gamma
    ebeam_rms["x_phase"] = rms_x_phase
    
    return ebeam_rms
    
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
    v = xp*np.sqrt(beta)+x*alpha/np.sqrt(beta)
    return [u,v]

def norm2real_coords(u,v,beta,alpha):
    x  = u*np.sqrt(beta)
    xp = (v-alpha*u)/np.sqrt(beta)
    return [x,xp]

def norm2act_coords(u,v,beta,alpha):
    J = np.sqrt(u**2 + v**2)/2
    phi = np.arctan2(-v,u)
    return [J,phi]

def real2act_coords(x,xp,beta,alpha):
    [u,v]   = real2norm_coords(x,xp,beta,alpha)
    [J,phi] = norm2act_coords(u,v,beta,alpha)
    return [J,phi]

def calc_frac_ellipse(u,v,frac=0.95,hires=False,tol=0.001):
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


    if frac<1.00:
        # center points
        Pc = np.zeros([len(P),2])
        Pc[:,0] = P[:,0] - cent[0]
        Pc[:,1] = P[:,1] - cent[1]
    
        # rotate points
        R = np.rot90(rot)
        Pr = np.zeros([len(Pc),2])
        for i in range(len(Pc)):
            Pr[i,:] = np.dot(R,Pc[i,:].T)
    
        # normalize points
        Pn = np.zeros([len(Pr),2])
        Pn[:,0] = Pr[:,0]/max(abs(Pr[:,0]))
        Pn[:,1] = Pr[:,1]/max(abs(Pr[:,1]))
        
        # calculate radius^2 of all points
#        Cn = [np.mean(Pn[:,0]),np.mean(Pn[:,1])]
        Cn = ([(np.median(Pn[:,0])+np.mean(Pn[:,0]))/2,\
               (np.median(Pn[:,1])+np.mean(Pn[:,1]))/2])
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
                # subtract 10% of target cut number at a pass
#                sub = round((1-frac)*npoint*0.10)
                # subtract 1 particle at a time
#                sub = 1
                if sub==0:
                    break
                Pt = Pt[:-sub,:]
                # calculate radius^2 of all points
#                Ct = [np.mean(Pt[:,0]),np.mean(Pt[:,1])]
                Ct = ([(np.median(Pn[:,0])+np.mean(Pn[:,0]))/2,\
                       (np.median(Pn[:,1])+np.mean(Pn[:,1]))/2])
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

        # un-normalize points
        Pt[:,0] = Pt[:,0]*max(abs(Pr[:,0]))
        Pt[:,1] = Pt[:,1]*max(abs(Pr[:,1]))
        # rotate back to normal orientation
        for i in range(len(Pt)):
            Pt[i,:] = np.dot(np.linalg.inv(R),Pt[i,:])
        # un-center points
        Pt[:,0] = Pt[:,0] + cent[0]
        Pt[:,1] = Pt[:,1] + cent[1]

        # use only points of convex hull for speed
        # NOTE: *could* be useful, but need to figure out
        # how to make it select fewer points!
    #    hull  = spatial.ConvexHull(Pt)
    #    Pt = hull.points
    #    Pt = np.unique(Pt,axis=0)
    
        # find minimum volume enclosing ellipsoid for remaining points
        ET = mvee.EllipsoidTool()
        (cent, rad, rot) = ET.getMinVolEllipse(Pt,tol)

    else:
        Pt = P

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
