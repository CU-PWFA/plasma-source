#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:29:57 2017

@author: chris
"""

import sys

import numpy as np
import GaussianBeam
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, "../python")

from ionization import ionization
from ionization import adk

def BuildDomain(size,step):
    return np.arange(-size/2,size/2,step)

def IntensityFromSpotSizes(wy,wz,xrange,yrange,zrange,step,I0,w0):
    I=np.empty([len(xrange),len(yrange),len(zrange)])
    i_counter=np.arange(len(xrange))
    j_counter=np.arange(len(yrange))
    k_counter=np.arange(len(zrange))
    for i in i_counter:
        for j in j_counter:
            for k in k_counter:
                I[i][j][k]=GaussianBeam.GaussianBeamIntensity_SpotArray_2D(I0,wy[i],wz[i],w0,yrange[j],zrange[k])
    return I

def IonFracFromIntensity(I,chi,delt_t):
    H=np.empty([len(I[:,0,0]),len(I[0,:,0]),len(I[0,0,:])])
    i_counter=np.arange(len(I[:,0,0]))
    j_counter=np.arange(len(I[0,:,0]))
    for i in i_counter:
        for j in j_counter:
            E = ionization.field_from_intensity(I[i][j]/1e14)
            H[i][j] = adk.gaussian_frac(chi,E,delt_t/1e-15,1)
    return H

#Give the data as a 2D plane of data[i][j] where we plot for i and vary j
def VarianceCut(data,axis,offset,number,spacing,unit,plot=['Title','x','f(x)','Legend'],both=False):
    center = round(len(data[0,:])/2)+offset
    
    plt.figure(figsize=(10,6))
    cmap = plt.get_cmap('jet_r')
    for i in np.arange(number):
        data_cut = data[:,center+i*spacing]
        col = cmap(float(i)/number)
        plt.plot(axis,data_cut,c=col,label=str(i*spacing*unit))  
        if both:
            data_cut = data[:,center-i*spacing]
            plt.plot(axis,data_cut,c=col)
    plt.title(plot[0])
    plt.xlabel(plot[1])
    plt.ylabel(plot[2])
    plt.legend(title=plot[3])
    plt.grid()
    plt.show()
    
def RobertRoll(data):
    return np.rollaxis(data,0,3)
    
def ImageCut(data,x,y,z,x_off=0,y_off=0,z_off=0,zoom=1,units='',label='',label_units='',color=0):
    xrange_z=[i*zoom for i in x]
    yrange_z=[i*zoom for i in y]
    zrange_z=[i*zoom for i in z]

    data_yz=np.transpose(data[round(len(x)/2)+x_off,:,:])
    data_xy=np.transpose(data[:,round(len(y)/2)+y_off,:])
    data_xz=np.transpose(data[:,:,round(len(z)/2)+z_off])
    
    if color == 1:
        plt.set_cmap('plasma')
    else:
        plt.set_cmap('viridis')
    
    gridSize = (2,5)
    plt.figure(figsize=(16,9))
    gridspec.GridSpec(gridSize[0], gridSize[1])
    
    plt.subplot2grid(gridSize, (0,0), rowspan=2)
    plt.imshow(data_yz, interpolation="none", origin="lower",
                 extent=[yrange_z[0],yrange_z[-1],zrange_z[0],zrange_z[-1]],
                 aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.ylabel('y '+units+' - Jet')
    plt.xlabel('x '+units+' - Beam')
    plt.title(label+'; x='+str(xrange_z[round(len(x)/2)+x_off])+' '+units)

    plt.subplot2grid(gridSize, (0,1), colspan=4)
    plt.imshow(data_xy, interpolation="none", origin="lower",
               extent=[xrange_z[0],xrange_z[-1],zrange_z[0],zrange_z[-1]],
               aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.xlabel('z '+units+' - Laser')
    plt.ylabel('y '+units+' - Jet')
    plt.title(label+'; y='+str(yrange_z[round(len(y)/2)+y_off])+' '+units)

    plt.subplot2grid(gridSize, (1,1), colspan=4)
    plt.imshow(data_xz, interpolation="none", origin="lower",
               extent=[xrange_z[0],xrange_z[-1],yrange_z[0],yrange_z[-1],],
               aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.xlabel('z '+units+' - Laser')
    plt.ylabel('x '+units+' - Beam')
    plt.title(label+'; z='+str(zrange_z[round(len(z)/2)+z_off])+' '+units)
    
    plt.tight_layout()
    plt.show()