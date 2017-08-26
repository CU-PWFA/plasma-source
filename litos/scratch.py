#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 14:49:29 2017

@author: mike
"""

import numpy as np
import beam_ana as ba
import matplotlib.pyplot as plt
import mvee as mvee

if __name__ == "__main__":
    
    #step = 2244
    #step = 0
    #step = 293 # vac beta min
    step = 709 # plas flat top start
    
    # get normalized coords
    [u,v] = ba.real2norm_coords(ebeam[step]["x"],ebeam[step]["xp"],\
                                ebeam[step]["beta"],ebeam[step]["alpha"])
    
    #u = ebeam[step]["x"]
    #v = ebeam[step]["xp"]
    
    # find the ellipse
    lips = ba.calc_frac_ellipse(u,v,frac=0.95,hires=False)
    rad  = lips["radii"]
    cent = lips["center"]
    rot  = lips["rot"]
    focd = lips["focd"]
    [foc1,foc2] = lips["foci"]
    ecc = lips["ecc"]
    area = lips["area"]
    Pt   = lips["Pt"]
    
    # print the fractional emittance
    emit = area*np.mean(ebeam[step]["gb"])/np.pi
    print('emittance: ',emit)
    
    
    # plot data points
    figX, axX1 = plt.subplots(1,1,sharey=True)
    axX1.scatter(u,v,c='r',s=0.5)
    axX1.scatter(Pt[:,0],Pt[:,1],c='g',s=0.5)
    axX1.set_xlim([-1.1*max(abs(u)),+1.1*max(abs(u))])
    axX1.set_ylim([-1.1*max(abs(v)),+1.1*max(abs(v))])
    
    # plot foci
    foc1 = np.dot([0,focd/2],rot)
    foc2 = -foc1
    axX1.scatter(foc1[0]+cent[0],foc1[1]+cent[1],c='b',s=10.0)
    axX1.scatter(foc2[0]+cent[0],foc2[1]+cent[1],c='b',s=10.0)
    
    # plot axes
    theta = np.arctan2(foc1[1],foc1[0])
    a = focd/np.sqrt(2-ecc**2)
    b = a*np.sqrt(1-ecc**2)
    majax_pt = [a*np.cos(theta),a*np.sin(theta)]
    minax_pt = [b*np.cos(theta+np.pi/2),b*np.sin(theta+np.pi/2)]
    
    majax = np.zeros([100,2])
    majax[:,0] = np.linspace(-majax_pt[0],+majax_pt[0],100)
    majax[:,1] = np.linspace(-majax_pt[1],+majax_pt[1],100)
    
    minax = np.zeros([100,2])
    minax[:,0] = np.linspace(-minax_pt[0],+minax_pt[0],100)
    minax[:,1] = np.linspace(-minax_pt[1],+minax_pt[1],100)
    
    axX1.plot(majax[:,0]+cent[0],majax[:,1]+cent[1],c='b',lw=0.5)
    axX1.plot(minax[:,0]+cent[0],minax[:,1]+cent[1],c='b',lw=0.5)
    
    # plot ellipse
    phi = np.linspace(0,2*np.pi,100)
    elx = rad[0]*np.cos(phi)
    ely = rad[1]*np.sin(phi)
    for i in range(len(elx)):
        [elx[i],ely[i]] = np.dot([elx[i],ely[i]], rot)
    axX1.plot(elx+cent[0],ely+cent[1],c='b',lw=1.0)
    
    plt.show()
