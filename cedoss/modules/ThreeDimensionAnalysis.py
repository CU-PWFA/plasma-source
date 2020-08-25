#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 13:29:57 2017
Contains functions for building 3D Intensity and Ionization fraction profiles,
and functions to plot cuts along the various axes and 2D planes.

Second half contains functions to fit various profiles to 1D tanh, 2D tanh,
or 1D Gaussian profiles; and plot their comparisons.

@author: chris
"""

import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import optimize

sys.path.insert(0, "../../python")

from ionization import ionization
from ionization import adk

#Modifies np.arange so that the window is centered at zero
def BuildDomain(size,step):
    return np.arange(-size/2,size/2,step)

#Takes a 3D intensity and returns a 3D ionization fraction
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

#Takes a 1D intensity and returns a 3D ionization fraction
#  I - 1D intensity in W/cm^2
#  chi - Ionization energy of species
#  delt_t - Temporal pulse length of beam in s
def IonFracFromIntensity_1D(I,chi,delt_t):
    H=np.zeros(len(I))
    E = ionization.field_from_intensity(I/1e14)
    H = adk.gaussian_frac(chi,E,delt_t/1e-15,1)
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

def VarianceCut_Prod(data,axis,offset,number,spacing,unit,plot=['Title','x','f(x)','Legend'],both=False):
    center = round(len(data[0,:])/2)+offset
    
    plt.figure(figsize=(10,6))
    plt.rcParams.update({'font.size': 16})
    cmap = plt.get_cmap('jet_r')
    for i in np.arange(number):
        data_cut = data[:,center+i*spacing]
        col = cmap(float(i)/number)
        if i == 0:
            plt.plot(axis,data_cut,c=col, label='0.00')
        else:
            plt.plot(axis,data_cut,c=col, label=r'$+$'+'%s' % float('%.3g' % (i*spacing*unit)))
        if both & (i>0):
            data_cut = data[:,center-i*spacing]
            plt.plot(axis,data_cut,'--',c=col, label=r'$-$'+'%s' % float('%.3g' % (i*spacing*unit)))
    #plt.title(plot[0])
    plt.xlabel(plot[1])
    plt.ylabel(plot[2])
    plt.legend(title=plot[3])
    #plt.grid()
    plt.tight_layout()
    #plt.savefig('/home/chris/Desktop/fig5.eps',format='eps',bbox_inches='tight',dpi=150)
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
    
def ImageCut_xy_Production(data,x,y,z,x_off=0,y_off=0,z_off=0,zoom=1,units='',label='',label_units='',color=0):
    xrange_z=[i*zoom for i in x]
    yrange_z=[i*zoom for i in y]
    zrange_z=[i*zoom for i in z]

    #data_yz=np.transpose(data[round(len(x)/2)+x_off,:,:])
    data_xz=np.transpose(data[round(len(x)*1/4):round(len(x)*3/4),round(len(y)/2)+y_off,:])*10
    #data_xy=np.transpose(data[:,:,round(len(z)/2)+z_off])
    
    if color == 1:
        plt.set_cmap('plasma')
    else:
        plt.set_cmap('viridis')
        
    
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}

    plt.rc('font', **font)
    plt.figure(figsize=(66.6,2))
    plt.imshow(data_xz, interpolation="none", origin="lower",
               extent=[xrange_z[0]/2,xrange_z[-1]/2,zrange_z[0],zrange_z[-1]],
               aspect='equal')
    CB=plt.colorbar(orientation = 'horizontal', pad =  .4)
    CB.set_label(label+' '+label_units)
    #plt.xlabel('x '+units+' - Laser')
    plt.ylabel('z '+units+' - Jet')
    plt.title(label+'; y='+str(yrange_z[round(len(y)/2)+y_off])+' '+units)
    #plt.savefig('demo.png', transparent=True,bbox_inches='tight')
    plt.savefig('/home/chris/Desktop/fig_strip.png',transparent=True,bbox_inches='tight')
    plt.show()

#Given the plasma density and 3D axes, return the total amount of charge in C
#   data = 3D array of plasma density (units 1e17 cm-3)
#   x,y,z = axes for each direction (units um)
def TotalCharge(data,x,y,z):
    sumden = 0
    e = 1.6022e-19
    volume = (x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])/(1e18)
    for i in range(len(data[:,0,0])):
        for j in range(len(data[0,:,0])):
            for k in range(len(data[0,0,:])):
                sumden = sumden + (data[i,j,k] * 1e6 * volume)
    return sumden * e * 1e17
    
"""
Data Fit Functions
"""    

#Generic function for a double tanh profile.  Flat in the center with ramps
# on either side which go to zero.  Centered around x = 0
#  p - parameters: [a (~top length), b (~ramp length), n_0 (density at center)]
#  x - array of distances
def DoubleTanh(p, x):
    return ((.5 + .5*np.tanh((x+p[0])/p[1])) * 
            (.5 - .5*np.tanh((x-p[0])/p[1]))) * p[2]
    
#Same as above, but a and b cannot be negative.  When a and b are negative
# there is a chance the fit can be a strange, Gaussian-like peak
def DoubleTanhAbs(p, x):
    return ((.5 + .5*np.tanh((x+abs(p[0]))/abs(p[1]))) *
            (.5 - .5*np.tanh((x-abs(p[0]))/abs(p[1])))) * p[2]
    
#DoubleTanh with a vertical offset
#  p - parameters: [a (~top length), b (~ramp length), n_0 (density at center)]
#  x - array of distances
def DoubleTanhOffset(p, x):
    return (((.5 + .5*np.tanh((x+p[0])/p[1])) * 
            (.5 - .5*np.tanh((x-p[0])/p[1])))+p[3]/p[2]) * p[2]
    
#Doubletanh with a flattop slant (for axial gas jet profiles)
#  p - parameters: [a (~top length), b (~ramp length), n_0 (density at center)]
#  x - array of distances
def DoubleTanhSlant(p, x):
    return ((.5 + .5*np.tanh((x+p[0])/p[1])) * 
            (.5 - .5*np.tanh((x-p[0])/p[1]))) * p[2] * (p[3]*x+1)
    
#Finds the effective distance for an ellipse, y^2 * (scl*z)^2
#  y,z - array of distances in the elliptical plane
#  scl - scaling factor between the narrow and wide waists
def EllipDist(y, z, scl):
    return np.sqrt(np.square(y) + np.square(scl * z))
    
#Generic function for Gaussian
#  p - parameters: [n_0 (density at center), sigma (distribution), x_0 (offset)]
#  x - array of distances
def Gaussian(p, x):
    return  p[0]*np.exp(-.5*np.square(x-p[2])/np.square(p[1]))

#Generic function for Gaussian
#  p - parameters: [n_0 (density at center), sigma (distribution), x_0 (offset)]
#  x - array of distances
def GaussianOffset(p, x):
    return  p[0]*(np.exp(-.5*np.square(x-p[2])/np.square(p[1]))+p[3]/p[0])

#Generic function for a Super Gaussian
#  p - parameters: [n_0 (density at center), sigma (distribution), x_0 (offset), p (power)]
#  x - array of distances
def SuperGaussian(p, x):
    return  p[0]*np.exp((-.5*np.square(x-p[2])/np.square(p[1])**p[3]))

#Generic function for Lorentzian
#  p - parameters: [A (density at center*factors), gamma (distribution), x_0 (offset)]
#  x - array of distances
def Lorentz(p, x):
    return  p[0]/(2*np.pi)*(p[1]/(np.square(x-p[2])+np.square(0.5*p[1])))

#A doubletanh multiplied by a Gaussian.  Disregards the n_0 for the Gaussian function.
#  p_tanh - parameters: [a (~top length), b (~ramp length), n_0 (density at center)]
#  p_gauss - parameters: [n_0 (density at center), sigma (distribution), x_0 (offset)]
#  x - array of distances
def DoubleTanh_Gaussian(p_tanh, p_gauss, x):
    p_gauss[0] = 1
    return DoubleTanh(p_tanh, x) * Gaussian(p_gauss, x)

#A doubletanh multiplied by a Lorentzian.  Disregards the n_0 for the DoubleTanh function.
#  p_tanh - parameters: [a (~top length), b (~ramp length), n_0 (density at center)]
#  p_lorentz - parameters: [A (density at center*factors), gamma (distribution), x_0 (offset)]
#  x - array of distances
def DoubleTanh_Lorentzian(p_tanh, p_lorentz, x):
    p_tanh[2] = 1
    return DoubleTanh(p_tanh, x) * Lorentz(p_lorentz, x)

#Basic exponential fit to f(x)=ae^(bx)+c, where p=[a,b,c]
def Exponential(p, x):
    return p[0]*np.exp(p[1]*x)+p[2]

#Lorentz, but with p[3] as constant offset in value
def LorentzOffset(p, x):
    return p[0]/(2*np.pi)*(p[1]/(np.square(x-p[2])+np.square(0.5*p[1])))+p[3]
    
#Attempts to fit a given 1D data set to a DoubleTanh Profile.
#  data - 1D data set for which to fit
#  axis - the independent variable which corresponds to data (ie: distance)
#  guess - inital guess of parameters [a, b, n_0]  (See DoubleTanh)
#  datlabel - what to call inital data in the plot's legend
#Returns the parameter array which fits the data
def FitDataDoubleTanh(data, axis, guess = [0.,0.,0.], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, DoubleTanh, guess, datlabel)
    print("a = " + str(p1[0]))
    print("b = " + str(p1[1]))
    print("n_0 = " + str(p1[2]))
    print()
    return p1

def FitDataDoubleTanhAbs(data, axis, guess = [0.,0.,0.], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, DoubleTanhAbs, guess, datlabel)
    p1[0] = abs(p1[0])
    p1[1] = abs(p1[1])
    print("a = " + str(p1[0]))
    print("b = " + str(p1[1]))
    print("n_0 = " + str(p1[2]))
    print()
    return p1

#Attempts to fit a given 1D data set to a Gaussian Profile.
#  data - 1D data set for which to fit
#  axis - the independent variable which corresponds to data (ie: distance)
#  guess - inital guess of parameters [n_0, sigma, x_0]  (See Gaussian)
#  datlabel - what to call inital data in the plot's legend
#Returns the parameter array which fits the data
def FitDataGaussian(data, axis, guess = [0.,0.,0.], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, Gaussian, guess, datlabel)
    print("n_0 = " + str(p1[0]))
    print("sig = " + str(p1[1]))
    print("x_0 = " + str(p1[2]))
    print()
    return p1

#Attempts to fit a given 1D data set to a Super Gaussian Profile.
#  data - 1D data set for which to fit
#  axis - the independent variable which corresponds to data (ie: distance)
#  guess - inital guess of parameters [n_0, sigma, x_0]  (See Gaussian)
#  datlabel - what to call inital data in the plot's legend
#Returns the parameter array which fits the data
def FitDataSuperGaussian(data, axis, guess = [0.,0.,0.,1], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, SuperGaussian, guess, datlabel)
    print("n_0 = " + str(p1[0]))
    print("sig = " + str(p1[1]))
    print("x_0 = " + str(p1[2]))
    print("pow = " + str(p1[3]))
    print()
    return p1

#Attempts to fit a given 1D data set to a Lorentzian Profile.
#  data - 1D data set for which to fit
#  axis - the independent variable which corresponds to data (ie: distance)
#  guess - inital guess of parameters [n_0, gamma, x_0]  (See Lorentz)
#  datlabel - what to call inital data in the plot's legend
#Returns the parameter array which fits the data
def FitDataLorentz(data, axis, guess = [0.,0.,0.], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, Lorentz, guess, datlabel)
    print("n0p = " + str(p1[0]))
    print("gam = " + str(p1[1]))
    print("x_0 = " + str(p1[2]))
    print()
    return p1

def FitDataExponential(data, axis, guess = [0.,0.,0.], datlabel = 'Simulation'):
    p1 = FitDataSomething(data, axis, Exponential, guess, datlabel)
    print("a = " + str(p1[0]))
    print("b = " + str(p1[1]))
    print("c = " + str(p1[2]))
    print()
    return p1

#Generic function to fit data to a function and plot results.  See above functions
def FitDataSomething(data, axis, function, guess = [0.,0.,0.], datlabel = 'Simulation',log=False,supress=False):
    errfunc = lambda p, x, y: function(p, x) - y
    p0 = guess
    p1, success = optimize.leastsq(errfunc,p0[:], args=(axis, data))
    if supress==False:
        if log==True:
            plt.semilogy(axis, data, label=datlabel)
        else:
            plt.plot(axis, data, label=datlabel)
        plt.plot(axis, function(p1,axis), label="Fitted "+ function.__name__ +" profile")
        plt.title("Comparison of data with "+ function.__name__ +" profile")
        plt.xlabel("Distance from axis (microns)")
        plt.ylabel("Density (10^17 cm^-3)")
        plt.legend()
        plt.grid()
        plt.show()
    return p1

#Plots 2D images of a 2D data from simulation and an analytical comparison given
# by the fits of n(y) and n(z).  Rather than fitting the full 2D profile, this 
# function assumes the fit of n(y) with an elliptical correction r = y + s*z,
# where the scaling factor s is the ratio between the first DoubleTanh parameters
# of n(y) and n(z).  Effectively, s is the ratio of flat top lengths
#  data - 2D data to be compared against.  Given in data[y][z]
#  y,zrange - axes of the y and z directions as arrays
#  py,pz - parameter arrays from the output of FitDataDoubleTanh (See DoubleTanh)
#  units - string of the units of the axes
#  label - string of what is plotted in the 2D plane
#  label_units - string of the units of label, what is plotted in the 2D plane
#Returns the 2D array of the difference between data and analytical
def Plot2DimDataTanh(data, yrange, zrange, py, pz, units='',label='',label_units='', slant=0):
    
    p1 = [py[0],py[1],py[2],py[0]/pz[0]]
        
    yaxis=np.reshape(yrange, (len(yrange), 1))
    zaxis=np.reshape(zrange, (1, len(zrange)))

    print("a = " + str(p1[0]))
    print("b = " + str(p1[1]))
    print("n_0 = " + str(p1[2]))
    print("scl = " + str(p1[3]))
    
    approx = DoubleTanh(p1,EllipDist(yaxis,zaxis,p1[3]))*(slant*zaxis+1)
    
    plt.set_cmap('plasma')
    gridSize = (2,3)
    plt.figure(figsize=(16,9))
    gridspec.GridSpec(gridSize[0], gridSize[1])
    
    plt.subplot2grid(gridSize, (0,0), rowspan=2)
    plt.imshow(np.transpose(data), interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],zaxis[0,0],zaxis[0,-1]],
                 aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.ylabel('z '+units+' - Jet')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Simulated ' + label)
    
    plt.subplot2grid(gridSize, (0,1), rowspan=2)
    plt.imshow(np.transpose(approx), interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],zaxis[0,0],zaxis[0,-1]],
                 aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.ylabel('z '+units+' - Jet')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Approximate ' + label)
    
    difference = np.transpose(data - approx)
    
    plt.subplot2grid(gridSize, (0,2), rowspan=2)
    plt.imshow(difference, interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],zaxis[0,0],zaxis[0,-1]],
                 aspect='auto')
    plt.set_cmap('viridis')
    CB=plt.colorbar()
    CB.set_label('Density Difference '+label_units)
    plt.ylabel('z '+units+' - Jet')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Density Difference (sim - fit)')
    
    plt.tight_layout()
    plt.show()
    return difference

#Plots 2D images of a 2D data from simulation and an analytical comparison given
# by the fits of n(y) and n(x).  Rather than fitting the full 2D profile, this 
# function assumes the 2D fit is n(y)*n(x)/n_0, or a product of a double tanh
# and a Gaussian.  
#  data - 2D data to be compared against.  Given in data[y][x]
#  y,xrange - axes of the y and x directions as arrays
#  py,px - parameter arrays from the output of FitDataDoubleTanh (See DoubleTanh)
#               and the output of FitDataGaussian (See Gaussian), respectively
#  units - string of the units of the axes
#  label - string of what is plotted in the 2D plane
#  label_units - string of the units of label, what is plotted in the 2D plane
#Returns the 2D array of the difference between data and analytical
def Plot2DimTanhxGaussian(data, yrange, xrange, py, px, units='',label='',label_units=''):
    print("a = " + str(py[0]))
    print("b = " + str(py[1]))
    print("n_0 = " + str(py[2]))
    print("sig = " + str(px[1]))
    print("x_0 = " + str(px[2]))
    
    yaxis=np.reshape(yrange, (len(yrange), 1))
    xaxis=np.reshape(xrange, (1, len(xrange)))
    
    ny = DoubleTanh(py,yaxis)
    nx = Gaussian(px,xaxis)
    approx = ny * nx / px[0]
    
    plt.set_cmap('plasma')
    gridSize = (2,3)
    plt.figure(figsize=(16,9))
    gridspec.GridSpec(gridSize[0], gridSize[1])
    
    plt.subplot2grid(gridSize, (0,0), rowspan=2)
    plt.imshow(np.transpose(data), interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],xaxis[0,0],xaxis[0,-1]],
                 aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.ylabel('z '+units+' - Laser')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Simulated ' + label)
    
    plt.subplot2grid(gridSize, (0,1), rowspan=2)
    plt.imshow(np.transpose(approx), interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],xaxis[0,0],xaxis[0,-1]],
                 aspect='auto')
    CB=plt.colorbar()
    CB.set_label(label+' '+label_units)
    plt.ylabel('x '+units+' - Laser')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Approximate ' + label)
    
    difference = np.transpose(data - approx)
    
    plt.subplot2grid(gridSize, (0,2), rowspan=2)
    plt.imshow(difference, interpolation="none", origin="lower",
                 extent=[yaxis[0,0],yaxis[-1,0],xaxis[0,0],xaxis[0,-1]],
                 aspect='auto')
    plt.set_cmap('viridis')
    CB=plt.colorbar()
    CB.set_label('Density Difference '+label_units)
    plt.ylabel('x '+units+' - Laser')
    plt.xlabel('y '+units+' - Beam')
    plt.title('Density Difference (sim - fit)')
    
    plt.tight_layout()
    plt.show()
    return difference

#Plots 2D contour of inital intensity before a refraction propagation
def PlotInitialIntensity(Ii,x,y):
    plt.figure(figsize=(5,5))
    plt.set_cmap('viridis')
    plt.imshow(np.abs(Ii),extent=[y[0],y[-1],x[0],x[-1]],aspect=y[-1]/x[-1])
    CB=plt.colorbar()
    CB.set_label('Intensity (10^14 W/cm^2)')
    plt.xlabel('Wide Axis: y (um)')
    plt.ylabel('Narrow Axis : x (um)')
    plt.title('Initial Intensity at edge')
    plt.show()

#I wanted this to work, but it just isnt at the moment...
#Currently bypassing optimization by assuming scl = p1[3] = a_y / a_z
# This works pretty well, so I am reluctant to fix 2D optimization
"""
def Fit2DimTanh(data,yrange,zrange, p0 = [0.,0.,0.,0.]):
    effdist = lambda y, z, scl: np.sqrt(np.square(y) + np.square(scl * z))
    
    fitfunc = lambda p, y, z: ((.5 + .5 * np.tanh((effdist(y, z, p[3]) + p[0])/p[1])) * 
                            (.5 - .5 * np.tanh((effdist(y, z, p[3]) - p[0])/p[1]))) * p[2]   
    errfunc = lambda p, y, z, val: fitfunc(p, y, z) - val
    p1, success = optimize.leastsq(errfunc, p0[:], args=(yrange, zrange, data))
    print("a = " + str(p1[0]))
    print("b = " + str(p1[1]))
    print("den = " + str(p1[2]))
    print("scl = " + str(p1[3]))
    return p1
"""