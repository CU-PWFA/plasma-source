#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 9 11:26:00 2017

Taken from public GitHub repository of:
@author: Michael Imelfort
Version: June 6, 2013

Modified by:
@author: mike
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg
from random import random
import scipy.spatial as spatial

class EllipsoidTool:
    """Some stuff for playing with ellipsoids"""
    def __init__(self): pass
    
    def getMinVolEllipse(self, P=None, tolerance=0.01):
        """ Find the minimum volume ellipsoid which holds all the points
        
        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!
        
        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]
        
        Returns:
        (center, radii, rotation)
        
        """
        
        hull = spatial.ConvexHull(P,qhull_options='QbBJ')
        P = hull.points
        P = np.unique(P,axis=0)
        
        (N, d) = np.shape(P)
        d = float(d)
    
        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)]) 
        QT = Q.T
        
        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse 
        center = np.dot(P.T, u)
    
        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) - 
                       np.array([[a * b for b in center] for a in center])
                       ) / d
                       
        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)
        
        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4./3.*np.pi*radii[0]*radii[1]*radii[2]

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False, cageColor='b', cageAlpha=0.2):
        """Plot an ellipsoid"""
        make_ax = ax == None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)
        
        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i,j],y[i,j]] = np.dot([x[i,j],y[i,j]], rotation) + center
    
        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0],0.0],
                             [0.0,radii[1]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)
    
    
            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                ax.plot(X3, Y3, color=cageColor)
    
        # plot ellipsoid
        ax.plot(x, y, color=cageColor, alpha=cageAlpha)
        
        if make_ax:
            plt.show()
#            plt.close(fig)
#            del fig
       


# example use 
if __name__ == "__main__":
    # make random points
    P = np.reshape([random()*100 for i in range(2000)],(1000,2))

#    hull = spatial.ConvexHull(P)
#    P = hull.points
#    P = np.unique(P,axis=0)

    # find the ellipsoid
    ET = EllipsoidTool()
    (center, radii, rotation) = ET.getMinVolEllipse(P, .01)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot points
    ax.scatter(P[:,0], P[:,1], color='g', marker='*', s=100)

    # plot ellipsoid
    ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=True)
    
    plt.show()
#    plt.close(fig)
#    del fig