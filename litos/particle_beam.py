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

def make_ebeam(s,i_twiss,i_parts):
    """create ebeam dictionary object"""
    ebeam = defaultdict(dict)
    append_ebeam_step(ebeam,0,s,i_twiss,i_parts)
    return ebeam

def append_ebeam_step(ebeam,step,s,i_twiss,i_parts):
    """append single step to ebeam dictionary object"""
    ebeam[step]["s"]     = s
    ebeam[step]["beta"]  = i_twiss["beta"]
    ebeam[step]["alpha"] = i_twiss["alpha"]
    ebeam[step]["gamma"] = i_twiss["gamma"]
    ebeam[step]["eps"]   = i_twiss["eps"]
    ebeam[step]["gbC"]   = i_twiss["gbC"]
    ebeam[step]["dgb"]   = i_twiss["dgb"]
    ebeam[step]["dz"]    = i_twiss["dz"]
    ebeam[step]["x"]     = i_parts["x"]
    ebeam[step]["xp"]    = i_parts["xp"]
    ebeam[step]["y"]     = i_parts["y"]
    ebeam[step]["yp"]    = i_parts["yp"]
    ebeam[step]["z"]     = i_parts["z"]
    ebeam[step]["gb"]    = i_parts["gb"]
    ebeam[step]["npart"] = i_parts["npart"]
    ebeam[step]["dist"]  = i_parts["dist"]
    return

def get_twiss(ebeam,step=0):
    """generate Twiss parameter dictionary object from ebeam"""
    twiss = defaultdict(dict)
    twiss[step]["beta"]  = ebeam[step]["beta"]
    twiss[step]["alpha"] = ebeam[step]["alpha"]
    twiss[step]["gamma"] = ebeam[step]["gamma"]
    twiss[step]["eps"]   = ebeam[step]["eps"]
    twiss[step]["gbC"]   = ebeam[step]["gbC"]
    twiss[step]["dgb"]   = ebeam[step]["dgb"]
    twiss[step]["dz"]    = ebeam[step]["dz"]
    return twiss

def get_parts(ebeam,step=0):
    """generate macro particle dictionary object from ebeam"""
    parts = defaultdict(dict)
    parts[step]["x"]     = ebeam[step]["x"]
    parts[step]["xp"]    = ebeam[step]["xp"]
    parts[step]["y"]     = ebeam[step]["y"]
    parts[step]["yp"]    = ebeam[step]["yp"]
    parts[step]["z"]     = ebeam[step]["z"]
    parts[step]["gb"]    = ebeam[step]["gb"]
    parts[step]["npart"] = ebeam[step]["npart"]
    parts[step]["dist"]  = ebeam[step]["dist"]
    return parts

def make_twiss(beta,alpha,gamma,eps,gbC,dgb,dz,step=0):
    """generate Twiss parameter dictionary object"""
    twiss = defaultdict(dict)
    twiss[step]["beta"]  = beta
    twiss[step]["alpha"] = alpha
    twiss[step]["gamma"] = gamma
    twiss[step]["eps"]   = eps
    twiss[step]["gbC"]   = gbC
    twiss[step]["dgb"]   = dgb
    twiss[step]["dz"]    = dz
    return twiss

def make_parts(i_twiss,npart,dist,step=0):
    """generate macro particle dictionary object"""
    parts = defaultdict(dict)

    # propagate individual particles in parallel
    num_cores = multiprocessing.cpu_count()
    phase6D = Parallel(n_jobs=num_cores)\
        (delayed(make_part)(i_twiss,dist,i)\
         for i in range(0,npart))
    phase6D = np.reshape(phase6D,[npart,6])

    # return transported beam params
    parts[step]["x"]     = phase6D[:,0]
    parts[step]["xp"]    = phase6D[:,1]
    parts[step]["y"]     = phase6D[:,2]
    parts[step]["yp"]    = phase6D[:,3]
    parts[step]["z"]     = phase6D[:,4]
    parts[step]["gb"]    = phase6D[:,5]
    parts[step]["npart"] = npart
    parts[step]["dist"]  = dist
    return parts

def make_part(i_twiss,dist,i=0):
    """generate single macro particle phase space parameters"""
    beta  = i_twiss["beta"]
    alpha = i_twiss["alpha"]
    gamma = i_twiss["gamma"]
    eps   = i_twiss["eps"]
    gbC   = i_twiss["gbC"]
    dgb   = i_twiss["dgb"]
    dz    = i_twiss["dz"]
    
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