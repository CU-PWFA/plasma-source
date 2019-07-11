#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 18:16:14 2019

@author: robert
"""

import numpy as np
import os
import re
import matplotlib.pyplot as plt

__path__ = None
__spec__ = None
datasets = {}
N = None
l = None

def setup(path, spec):
    """ Find all the files and setup variables for spectrometer analysis.
    
    Parameters
    ----------
    path : string
        The path containing the dataset folder.
    spec : string
        The id of the spectrometer.
    """
    global __path__ 
    global __spec__
    global datasets
    global l
    global N
    __path__ = path

    __spec__ = spec
    directories = os.listdir(__path__)
    dataset = None
    # We do all this in case two spectrometers with different id's are connected
    for i in range(len(directories)):
        dataset = directories[i]
        rex = re.compile(spec+'_'+dataset+'_(?P<shot>.*).npy')
        file_names = os.listdir(__path__+directories[i])
        shots = []
        for i in range (len(file_names)):
            match = rex.search(file_names[i])
            if match:
                shot = int(match.group('shot'))
                shots.append(shot)
        if len(shots) != 0:
            datasets[dataset] = {'shots' : shots}
    # Load one dataset to get the wavelengths
    if dataset is not None and len(shots) > 0:
        spectrum = np.load(path+dataset+'/michaelito_'+dataset+'_0000.npy').item()
        l = spectrum['lambda']
        N = len(l)
    else:
        print('No non-empty datasets found')
        return
    print('Spectrometer range %0.2fnm to %0.2fnm' %(l[0], l[-1]))
    print('Number of pixels', N)
    print('Datasets found:', datasets.keys())
    return l, N

def integrate_dataset(dataset):
    """ Add all the shots in a dataset together.
    
    Parameters
    ----------
    dataset : int or string
        The dataset number to add everything together for.
    """
    dataset = str(dataset)
    if dataset not in datasets:
        print("Dataset number was not detected in the files.")
        return None
    file = __path__ + str(dataset) + '/' + __spec__ +'_' + str(dataset) + '_'
    I = np.zeros(N, dtype='double')
    for i in datasets[dataset]['shots']:
        filename = file + '%04d.npy' % i
        I += np.load(filename).item()['I']
    return I

def plot_spectrum(*args, xlim=None, ylim=None, lines=None, title=None):
    """ Plot the passed spectrums with the passed lines. 
    
    Parameters
    ----------
    *args : arrays of floats
        The spectrum to plot.
    xlim, optional : tuple or array
        Start and end values for the x axis.
    ylim, optional : tuple or array
        Start and end for the y axis.
    lines, optional : array of floats
        Spectral lines to mark on the plot.
    title, optional : string
        Title of the plot.
    """
    if xlim is None:
        xlim = (l[0], l[-1])
    plt.figure(figsize=(8, 2), dpi=150)
    ax = plt.subplot()
    peak = 0
    colors = [plt.cm.brg(i) for i in np.linspace(0, 1, len(args))]
    ax.set_prop_cycle('color', colors)
    for arg in args:
        plt.plot(l, arg, linewidth=0.2)
        amax = np.amax(arg)
        if amax > peak:
            peak = amax
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('Counts')
    plt.xlim(xlim)
    if ylim is None:
        ylim = (-1000, 1.05*peak)
    plt.ylim(ylim)
    if lines is not None:
        for i in range(len(lines)):
            line = lines[i]
            plt.plot([line, line], [ylim[0], ylim[1]], 'k--', linewidth=0.2)
    if title is not None:
        plt.title(title)
    plt.show()
    
def plot_line(intensity, density, lam, xlim=None, ylim=None, lines=None):
    """ Plot a given line for several datasets, normalized. 
    
    Parameters
    ----------
    intensity : 2D array
        First dimension is datasets at different densities, second is spectrums.
    density : array of floats
        Density for each dataset.
    lam : float
        Wavelength (nm) of the line to look at.
    xlim, optional : tuple or array
        Start and end values for the x axis.
    ylim, optional : tuple or array
        Start and end for the y axis.
    lines, optional : array of floats
        Spectral lines to mark on the plot.
    """
    M = np.shape(intensity)[0]
    if xlim is None:
        xlim = (lam-10, lam+10)
    if ylim is None:
        ylim = (0, 10000)
    
    plt.figure(figsize=(8, 3), dpi=150)
    ax = plt.subplot()
    colors = [plt.cm.brg(i) for i in np.linspace(0, 1, M+4)]
    ax.set_prop_cycle('color', colors)
    linewidth = 0.4
    for i in range(M):
        plt.plot(l, intensity[i], linewidth=linewidth, label='%0.2E' % density[i])
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('Counts')
    plt.xlim(xlim)
    plt.ylim(ylim)
    if lines is not None:
        for i in range(len(lines)):
            line = lines[i]
            plt.plot([line, line], [ylim[0], ylim[1]], 'k--', linewidth=0.2)
    plt.legend(title='Densities ($cm^{-3}$)', loc='upper left', prop={'size':8})
    plt.show()
    
    # Find the index of the line
    ind = np.argmin(abs(l - lam))
    ylim = (0, 1.2)
    plt.figure(figsize=(8, 3), dpi=150)
    ax = plt.subplot()
    ax.set_prop_cycle('color', colors)
    for i in range(M):
        plt.plot(l, intensity[i]/intensity[i][ind], linewidth=linewidth, label='%0.2E' % density[i])
    plt.xlabel('$\lambda$ (nm)')
    plt.ylabel('Counts')
    plt.xlim(xlim)
    plt.ylim(ylim)
    if lines is not None:
        for i in range(len(lines)):
            line = lines[i]
            plt.plot([line, line], [ylim[0], ylim[1]], 'k--', linewidth=0.2)
    plt.legend(title='Densities ($cm^{-3}$)', loc='upper left', prop={'size':8})
    plt.show()