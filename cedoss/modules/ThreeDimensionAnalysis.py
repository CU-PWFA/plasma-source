#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:29:57 2017
Contains functions for building 3D Intensity and Ionization fraction profiles,
and functions to plot cuts along the various axes and 2D planes

@author: chris
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, "../../python")

from ionization import ionization
from ionization import adk

#Modifies np.arange so that the window is centered at zero
def BuildDomain(size,step):
    return np.arange(-size/2,size/2,step)

#Takes a 3D intensity and returns a 3D ionization fraction
#!!I want to put this in a better place eventually!!
#  I - 3D intensity in W/cm^2
#  chi - Ionization energy of species
#  delt_t - Temporal pulse length of beam in s
def IonFracFromIntensity(I,chi,delt_t):
    H=np.empty([len(I[:,0,0]),len(I[0,:,0]),len(I[0,0,:])])
    i_counter=np.arange(len(I[:,0,0]))
    j_counter=np.arange(len(I[0,:,0]))
    for i in i_counter:
        for j in j_counter:
            E = ionization.field_from_intensity(I[i][j]/1e14)
            H[i][j] = adk.gaussian_frac(chi,E,delt_t/1e-15,1)
    return H

#If using data from roberts code, we need to roll the axis to fit my convention
# Robert: [beam,jet,laser]   Doss: [laser,beam,jet]
#  data - 3D array in [beam,jet,laser]
def RobertRoll(data):
    return np.rollaxis(data,2)

#Takes a 3D array, and returns a list containing the x,y,z offset of array
# indices from the middle to the maximum value, as well as the maximum value
#  data - 3D array
#  returns - [x_off,y_off,z_off,max(data)]
def GetMaximumOffset(data):
    ret = [0,0,0,0]
    for i in range(len(data[:,0,0])):
        for j in range(len(data[0,:,0])):
            for k in range(len(data[0,0,:])):
                if data[i,j,k] > ret[3]:
                    ret=[i,j,k,data[i,j,k]]
    return [ret[0]-round(len(data[:,0,0])/2),
            ret[1]-round(len(data[0,:,0])/2), 
            ret[2]-round(len(data[0,0,:])/2), ret[3]]

#Take a 2D array of data, plots cuts in one direction while varying the other
#  data - 2D array of [i][j], where i is what is plotted and j is varied
#  axis - the range of which [i] is plotted against
#  offset - number of array indices off the center of [j] to offset origin by
#  number - number of lines to plot iteratively further from origin
#  spacing - number of array indices that separate sequential lines
#  unit - factor relating separation distance per array index
#  plot - labels for the plot, axes, and legend
#  both - set to True to interate in both directions from the origin
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

#Takes a 3D array, and plots the 2D planes cut along the 3 axes
#  data - 3D array in [laser,beam,jet]
#  x,y,z - axes in x,y,z
#  x,y,z_off - number of array indices to offset the axes by
#  zoom - factor to multiply the axes by
#  units - String for the units of the axes
#  label - String for the data that is being plotted
#  label_units - String for the units of the data being plotted
#  color - set to 1 for 'plasma', otherwise defaults to 'viridis'
def ImageCut(data,x,y,z,x_off=0,y_off=0,z_off=0,zoom=1,units='',label='',label_units='',color=0):
    xrange_z=[i*zoom for i in x]
    yrange_z=[i*zoom for i in y]
    zrange_z=[i*zoom for i in z]

    data_yz=np.transpose(data[round(len(x)/2)+x_off,:,:])
    data_xz=np.transpose(data[:,round(len(y)/2)+y_off,:])
    data_xy=np.transpose(data[:,:,round(len(z)/2)+z_off])
    
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
    plt.ylabel('z '+units+' - Jet')
    plt.xlabel('y '+units+' - Beam')
    plt.title(label+'; x='+str(xrange_z[round(len(x)/2)+x_off])+' '+units)

    plt.subplot2grid(gridSize, (0,1), colspan=4)
    plt.imshow(data_xz, interpolation="none", origin="lower",
               extent=[xrange_z[0],xrange_z[-1],zrange_z[0],zrange_z[-1]],
               aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.xlabel('x '+units+' - Laser')
    plt.ylabel('z '+units+' - Jet')
    plt.title(label+'; y='+str(yrange_z[round(len(y)/2)+y_off])+' '+units)

    plt.subplot2grid(gridSize, (1,1), colspan=4)
    plt.imshow(data_xy, interpolation="none", origin="lower",
               extent=[xrange_z[0],xrange_z[-1],yrange_z[0],yrange_z[-1],],
               aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.xlabel('x '+units+' - Laser')
    plt.ylabel('y '+units+' - Beam')
    plt.title(label+'; z='+str(zrange_z[round(len(z)/2)+z_off])+' '+units)
    
    plt.tight_layout()
    plt.show()