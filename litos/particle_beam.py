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
    ebeam = append_ebeam_step(ebeam,0,s,twiss,parts)
    return ebeam

def append_ebeam_step(ebeam,step,s,twiss,parts):
    ebeam[step]["s"]     = s
    ebeam[step]["beta"]  = twiss["beta"]
    ebeam[step]["alpha"] = twiss["alpha"]
    ebeam[step]["gamma"] = twiss["gamma"]
    ebeam[step]["eps"]   = twiss["eps"]
    ebeam[step]["gbC"]   = twiss["gbC"]
    ebeam[step]["dgb"]   = twiss["dgb"]
    ebeam[step]["dz"]    = twiss["dz"]
    ebeam[step]["x"]     = parts["x"]
    ebeam[step]["xp"]    = parts["xp"]
    ebeam[step]["y"]     = parts["y"]
    ebeam[step]["yp"]    = parts["yp"]
    ebeam[step]["z"]     = parts["z"]
    ebeam[step]["gb"]    = parts["gb"]
    ebeam[step]["npart"] = parts["npart"]
    ebeam[step]["dist"]  = parts["dist"]
    return ebeam

def get_twiss(ebeam,step=0):
    """generate Twiss parameter dictionary object from ebeam"""
    twiss = defaultdict(dict)
    twiss["beta"]  = ebeam[step]["beta"]
    twiss["alpha"] = ebeam[step]["alpha"]
    twiss["gamma"] = ebeam[step]["gamma"]
    twiss["eps"]   = ebeam[step]["eps"]
    twiss["gbC"]   = ebeam[step]["gbC"]
    twiss["dgb"]   = ebeam[step]["dgb"]
    twiss["dz"]    = ebeam[step]["dz"]
    return twiss

def get_parts(ebeam,step=0):
    """generate macro particle dictionary object from ebeam"""
    parts = defaultdict(dict)
    parts["x"]     = ebeam[step]["x"]
    parts["xp"]    = ebeam[step]["xp"]
    parts["y"]     = ebeam[step]["y"]
    parts["yp"]    = ebeam[step]["yp"]
    parts["z"]     = ebeam[step]["z"]
    parts["gb"]    = ebeam[step]["gb"]
    parts["npart"] = ebeam[step]["npart"]
    parts["dist"]  = ebeam[step]["dist"]
    return parts

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

    # propagate individual particles in parallel
    num_cores = multiprocessing.cpu_count()
    phase6D = Parallel(n_jobs=num_cores)\
        (delayed(make_part)(twiss,dist,i)\
         for i in range(0,npart))
    phase6D = np.reshape(phase6D,[npart,6])

    # return transported beam params
    parts["x"]     = phase6D[:,0]
    parts["xp"]    = phase6D[:,1]
    parts["y"]     = phase6D[:,2]
    parts["yp"]    = phase6D[:,3]
    parts["z"]     = phase6D[:,4]
    parts["gb"]    = phase6D[:,5]
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