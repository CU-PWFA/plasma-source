#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 20:33:52 2021

@author: valentinalee
"""

import numpy as np
from scipy import optimize
from sympy import symbols, Eq, solve
from scipy.interpolate import griddata
import os

def Analysis(Path, DataNumber):
    '''
    Path: str, path to where 02_MultishotsAvgResults dir and 03_PeaksFoundByClickyTool dir are at
    ----------
    return:
    PV_avg: float, average peak-valley value
    Separation: First Peaks Separation at the middle of the signal
    '''
    def y_valley(x):
        y_valley= ((-a2/np.sqrt(a2**2+1))-(a1/np.sqrt(a1**2+1)))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))*x+ \
            ((-b2/np.sqrt(a2**2+1))-(b1/(np.sqrt(a1**2+1))))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))
        return y_valley
    
    def ResultProcess(x_valley, data, BG):
        '''
        x_valley: int, the x coord. of the valley line
        data: 2D array, the 2D array of data (or signal) that is being processed
        --------
        return: 
        ResultList: 2D List, ([x_valley, Intensity of Peak1, Intensity of Peak2, Intensity of Valley, Distance of Peak1 and Peak2])
        '''
#        global a_bis_v, x_grid, y_grid, a1, b1, a2, b2, x, y, griddata
        y_current= y_valley(x_valley)
        b_bis_v_current= y_current- a_bis_v*x_valley
        
        eq0= Eq(a_bis_v*x- y+ b_bis_v_current)
        eq1= Eq(a1*x- y+ b1)
        eq2= Eq(a2*x- y+ b2)
        x1y1= solve((eq0,eq1), (x, y))
        x2y2= solve((eq0,eq2), (x, y))
        
        points= np.array((x_grid.flatten(), y_grid.flatten())).T
        values= data.flatten()
        BGvalues= BG.flatten()
        z0= griddata(points, values, (x_valley, y_current))
        z1= griddata(points, values, (x1y1[x], x1y1[y]))
        z2= griddata(points, values, (x2y2[x], x2y2[y]))
        z3= griddata(points, BGvalues, (x_valley, y_current))
        dist= np.sqrt((float(x1y1[x])-float(x2y2[x]))**2+ (float(x1y1[y])-float(x2y2[y]))**2)
        return [x_valley, float(z1), float(z2), float(z0), dist, float(z3)]
    
    
    Path= os.path.expanduser(Path)
    WhateverData= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_Signal.npy')
#    WhateverData= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_Data.npy')
    BG= np.load(Path+'/02_MultishotsAvgResults'+'/01_npy'+'/'+DataNumber+'_BG.npy')
    Coords= np.load(Path+'/03_PeaksFoundByClickyTool'+'/'+DataNumber+'_PeaksInfo.npy')

    X= np.linspace(0, 499, 500)
    Y= np.linspace(0, 499, 500)
    x_grid, y_grid= np.meshgrid(X, Y)
    
    
    
    Peak1X= Coords[:, 0][Coords[:, 0]!=0]
    Peak1Y= Coords[:, 1][Coords[:, 1]!=0]
    Peak2X= Coords[:, 2][Coords[:, 2]!=0]
    Peak2Y= Coords[:, 3][Coords[:, 3]!=0]
    
    xvar= Peak1X
    yvar= Peak1Y
    fitfn= lambda p, xvar: p[0]*xvar+p[1]
    
    errfunc= lambda p, xvar, yvar: fitfn(p, xvar)- yvar
    p0= [-1, 300]
    p1, success= optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    a1= p1[0]
    b1= p1[1]
    
    xvar= Peak2X
    yvar= Peak2Y
    fitfn= lambda p, xvar: p[0]*xvar+p[1]
    
    errfunc= lambda p, xvar, yvar: fitfn(p, xvar)- yvar
    p0= [-1, 0]
    p1, success= optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
    a2= p1[0]
    b2= p1[1]
    
    a_bis= ((-a2/np.sqrt(a2**2+1))-(a1/np.sqrt(a1**2+1)))/((-1/np.sqrt(a2**2+1))-(1/np.sqrt(a1**2+1)))

    a_bis_v= -1/a_bis
    x, y= symbols('x y')
    ResultList= [ResultProcess(x, WhateverData, BG) for x in range(150, 350)]
#TODO Change a way to save the result so you can plot it in different way before you finalize how you are going to show the result
#    PV= [(ResultList[x][1]+ResultList[x][2])/2-ResultList[x][3] for x in range(0, np.shape(ResultList)[0])]
    np.save(Path+'/04_Results/'+DataNumber+'_ResultList', np.array(ResultList))
    np.save(Path+'/04_Results/'+DataNumber+'_ab', np.array([a1, b1, a2, b2]))

    PV= [((ResultList[x][1]+ResultList[x][2])/2-ResultList[x][3])/ResultList[x][5] for x in range(0, np.shape(ResultList)[0])]
    PV_avg= sum(PV)/ len(PV)
    Angle= np.arctan((a1-a2)/(1+a1*a2))
#    Separation= ResultList[int(np.shape(ResultList)[0]/2)][4]

    return PV_avg, Angle
#    return ResultList

