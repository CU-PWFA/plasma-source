#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 14:17:52 2017

@author: mike
"""

import numpy as np
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing
import time as time

def make_ebeam(s,twiss,parts):
    ebeam = defaultdict(dict)
    ebeam[0]["s"]     = s
    ebeam[0]["twiss"] = twiss
    ebeam[0]["parts"] = parts
    return ebeam

def append_ebeam_step(ebeam,step,s,twiss,parts):
    ebeam[step]["s"]     = s
    ebeam[step]["twiss"] = twiss
    ebeam[step]["parts"] = parts
    return ebeam

def make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz):
    """generate Twiss parameter dictionary object"""
    twiss = defaultdict(dict)
    twiss["beta"]  = beta
    twiss["alpha"] = alpha
    twiss["gamma"] = gamma
    twiss["eps"]   = eps
    twiss["gbC"]   = gbC
    twiss["dgb"]   = dgb
    twiss["dz"]    = dz
    return twiss

def make_parts(twiss,npart,dist):
    """generate macro particle dictionary object"""
    parts = defaultdict(dict)

    """Old, non-parallel way:
#    beta  = twiss["beta"]
#    alpha = twiss["alpha"]
#    gamma = twiss["gamma"]
#    eps   = twiss["eps"]
#    gbC   = twiss["gbC"]
#    dgb   = twiss["dgb"]
#    dz    = twiss["dz"]
#         
#    x  = np.zeros(npart)
#    xp = np.zeros(npart)
#    y  = np.zeros(npart)
#    yp = np.zeros(npart)
#    z  = np.zeros(npart)
#    gb = np.zeros(npart)
#    
#    for i in range(0,npart):
#        rndx = np.random.uniform(0,1)
#        rndy = np.random.uniform(0,1)
#        if dist.lower() == 'gauss': # gaussian distribution
#            rx   = np.sqrt(2)*np.sqrt(eps/gbC)*\
#                    np.sqrt(-np.log(rndx))
#            ry   = np.sqrt(2)*np.sqrt(eps/gbC)*\
#                    np.sqrt(-np.log(rndy))
#        elif dist.lower() == 'uniform': # uniform distribution
#            rx   = np.sqrt(rndx)*np.sqrt(eps/gbC)
#            ry   = np.sqrt(rndy)*np.sqrt(eps/gbC)
#        else :
#            print('no particle distribution given')
#            return
#                    
#        phix = np.random.uniform(0,2*np.pi)
#        ux = rx*np.cos(phix)
#        vx = rx*np.sin(phix)
#        x[i]  = ux*np.sqrt(beta)
#        xp[i] = (vx-(alpha/beta)*x[i])/np.sqrt(beta)
#        
#        phiy = np.random.uniform(0,2*np.pi)
#        uy = ry*np.cos(phiy)
#        vy = ry*np.sin(phiy)
#        y[i]  = uy*np.sqrt(beta)
#        yp[i] = (vy-(alpha/beta)*y[i])/np.sqrt(beta)
#                
#        z[i]  = dz*np.random.uniform(-1,1)
#        gb[i] = gbC*(1+dgb*np.random.uniform(-1,1))
#
#    parts["x"]  = x
#    parts["xp"] = xp
#    parts["y"]  = y
#    parts["yp"] = yp
#    parts["z"]  = z
#    parts["gb"] = gb
#    parts["npart"] = npart
#    parts["dist"]  = dist
    """
            
    # propagate individual particles in parallel
    num_cores = multiprocessing.cpu_count()
    phase6D = Parallel(n_jobs=num_cores)\
        (delayed(make_part)(twiss,dist,i)\
         for i in range(0,npart))
    phase6D = np.reshape(phase6D,[npart,6])

    # return transported beam params
    parts["x"]  = phase6D[:,0]
    parts["xp"] = phase6D[:,1]
    parts["y"]  = phase6D[:,2]
    parts["yp"] = phase6D[:,3]
    parts["z"]  = phase6D[:,4]
    parts["gb"] = phase6D[:,5]
    parts["npart"] = npart
    parts["dist"]  = dist
    return parts

def make_part(twiss,dist,i=0):
    """generate single macro particle phase space parameters"""
    beta  = twiss["beta"]
    alpha = twiss["alpha"]
    gamma = twiss["gamma"]
    eps   = twiss["eps"]
    gbC   = twiss["gbC"]
    dgb   = twiss["dgb"]
    dz    = twiss["dz"]
    
    iseed = int((time.time()%1e2)*1e7 + i)
    np.random.seed(iseed)    
    rndx = np.random.uniform(0,1)
    rndy = np.random.uniform(0,1)
    if dist.lower() == 'gauss': # gaussian distribution
        rx   = np.sqrt(2)*np.sqrt(eps/gbC)*\
                np.sqrt(-np.log(rndx))
        ry   = np.sqrt(2)*np.sqrt(eps/gbC)*\
                np.sqrt(-np.log(rndy))
    elif dist.lower() == 'uniform': # uniform distribution
        rx   = np.sqrt(rndx)*np.sqrt(eps/gbC)
        ry   = np.sqrt(rndy)*np.sqrt(eps/gbC)
    else :
        print('no particle distribution given')
        return
                
    phix = np.random.uniform(0,2*np.pi)
    ux = rx*np.cos(phix)
    vx = rx*np.sin(phix)
    x  = ux*np.sqrt(beta)
    xp = (vx-(alpha/beta)*x)/np.sqrt(beta)
    
    phiy = np.random.uniform(0,2*np.pi)
    uy = ry*np.cos(phiy)
    vy = ry*np.sin(phiy)
    y  = uy*np.sqrt(beta)
    yp = (vy-(alpha/beta)*y)/np.sqrt(beta)
            
    z  = dz*np.random.uniform(-1,1)
    gb = gbC*(1+dgb*np.random.uniform(-1,1))
    
    phase6D = [x,xp,y,yp,z,gb]
    return phase6D