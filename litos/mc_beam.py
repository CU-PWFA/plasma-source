#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:30:34 2017

@author: litos
"""

import numpy as np

def mc_beam(Tbeam,npart=100,dE=0):
    
    [gb0,eps,beta,alpha,gamma] = Tbeam
    
    u  = np.zeros(npart)
    v  = np.zeros(npart)
    x  = np.zeros(npart)
    xp = np.zeros(npart)
    gb = np.zeros(npart)
        
    for i in range(0,npart):
        s   = np.random.uniform(0,1)
        
#        r   = np.sqrt(s)*np.sqrt(eps/gb0) # square distribution
        r   = np.sqrt(2)*np.sqrt(eps/gb0)*\
                np.sqrt(-np.log(s)) # gaussian distribution

        phi = np.random.uniform(0,2*np.pi)
        
        u[i] = r*np.cos(phi)
        v[i] = r*np.sin(phi)
    
        x[i]  = u[i]*np.sqrt(beta)
        xp[i] = (v[i]-(alpha/beta)*x[i])/np.sqrt(beta)
                
        gb[i] = gb0*(1+dE*np.random.uniform(-1,1))

    mcb = [x,xp,u,v,gb]
    
    return mcb