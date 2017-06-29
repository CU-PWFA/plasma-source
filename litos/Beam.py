#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 11:34:00 2017

@author: mike
"""
import numpy as np

class Beam(object):
    """Electron Beam"""
    def __init__(self,gb0=0,eps0=0,beta0=0,alpha0=0,gamma0=0,dz0=0,dgb0=0,npart=0,dist='gauss'):
        """Initialize Beam"""
        self.gbC   = gb0
        self.eps   = eps0
        self.beta  = beta0
        self.alpha = alpha0
        self.gamma = gamma0
        self.dz    = dz0
        self.dgb   = dgb0
        self.npart = npart
        self.dist  = dist
        if self.npart>0:
            self.gen_parts(self.npart,self.dist)
        else:
            self.parts = []
            
    def gen_parts(self,npart,dist):
        """generate macro-particle distribution"""
        x  = np.zeros(npart)
        xp = np.zeros(npart)
        y  = np.zeros(npart)
        yp = np.zeros(npart)
        z  = np.zeros(npart)
        gb = np.zeros(npart)
            
        for i in range(0,npart):
            rndx = np.random.uniform(0,1)
            rndy = np.random.uniform(0,1)
            if dist.lower() == 'gauss': # gaussian distribution
                rx   = np.sqrt(2)*np.sqrt(self.eps/self.gbC)*\
                        np.sqrt(-np.log(rndx))
                ry   = np.sqrt(2)*np.sqrt(self.eps/self.gbC)*\
                        np.sqrt(-np.log(rndy))
            elif dist.lower() == 'uniform': # uniform distribution
                rx   = np.sqrt(rndx)*np.sqrt(self.eps/self.gbC)
                ry   = np.sqrt(rndy)*np.sqrt(self.eps/self.gbC)
            else :
                print('no particle distribution given')
                return
                        
            phix = np.random.uniform(0,2*np.pi)
            ux = rx*np.cos(phix)
            vx = rx*np.sin(phix)
            x[i]  = ux*np.sqrt(self.beta)
            xp[i] = (vx-(self.alpha/self.beta)*x[i])/np.sqrt(self.beta)
            
            phiy = np.random.uniform(0,2*np.pi)
            uy = ry*np.cos(phiy)
            vy = ry*np.sin(phiy)
            y[i]  = uy*np.sqrt(self.beta)
            yp[i] = (vy-(self.alpha/self.beta)*y[i])/np.sqrt(self.beta)
                    
            z[i]  = self.dz*np.random.uniform(-1,1)
            gb[i] = self.gbC*(1+self.dgb*np.random.uniform(-1,1))
    
            self.parts = np.zeros((npart,6))
            for i in range(0,npart):
                self.parts[i][:] = [x[i],xp[i],y[i],yp[i],z[i],gb[i]]

    def propagate(self,s,np,dgds):
        """propagate beam"""
        
        

        
# define initial beam
gb0    = 20000 # relativistic lorentz factor
eps0   = 5e-6  # m-rad, normalized emittance

beta0  = 0.10 # m
alpha0 = 0.00
gamma0 = (1.0+alpha0**2)/beta0 # 1/m
dE0    = 0.01
npart  = 3     

        
beam = Beam(gb0,eps0,beta0,alpha0,gamma0,0,dE0,npart,'gauss')
print(beam.beta)
print(beam.parts[1,:])