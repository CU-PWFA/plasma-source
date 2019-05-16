#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:26:34 2019

@author: robert
"""

import sys
import numpy as np
from gdsCAD import *
from scipy.interpolate import interp1d


def main(fileName, layers, min_feature=2):
    """ Convert a numpy array of radius and phase values to a mask.
    
    Wraps and discretizes the phase in the input file and generates a gds file
    with the appropriate mask writing for phase information. Radius should be
    in um and phase in rad.
    
    Note: one day this should be switched to gdspy which is actively maintained
    and works with python 3. It also doesn't require modifying the source code
    in order to get it to work.
    
    Parameters
    ----------
    fileName : string
        Name of the file to convert, will also be the name of the output file.
    layers : int
        The number of discrete phase steps, the GDS will have Nlayers=layers/2.
    min_feature : float
        The minimum feature size in um, small features will be converted.
    """
    print("Converting the array in "+fileName+" to a GDS file")
    data = np.load(fileName)
    print("Successfully loaded "+fileName)
    try:
        r = data[0]
        phi = data[1]
    except:
        print("The numpy file doesn't have the correct format, it should be an array with:")
        print("[r, phi]")
    print("Converting phase into layered rings")
    phif = interp1d(r, phi, kind='cubic')
    print("Creating nm resolution arrays")
    N = int(r[-1]*1000)
    # Radius array with exactly one point per nm (nm is base scale of gds)
    r_nm = np.linspace(0, N/1000, N, endpoint=False)
    phi_nm = phif(r_nm)
    phi_nm = phi_nm % (2*np.pi)
    phi_d = np.floor(phi_nm*layers/(2*np.pi))
    
    # Parse the phase to find the start and end of each ring
    rings = {0 : {'r':0, 'l':0}}
    j = 0
    phi_cur = 0
    for i in range(len(r_nm)):
        if phi_d[i] != phi_cur:
            j += 1
            phi_cur = phi_d[i]
            r_tran = np.round(r_nm[i], 3)
            rings[j] = {'r':r_tran,'l':phi_cur}
            rings[j-1]['rf'] = r_tran
    # Have to set the outside of the last ring and make sure it is always large enough
    rings[j]['rf'] = r[-1]+min_feature
    print("Finished parsing ring boundaries")
    
    print("Converting any undersized features")
    del rings[0]
    delta = np.zeros(len(rings)-1)
    for i in range(1, len(rings)):
        delta[i-1] = rings[i+1]['r']-rings[i]['r']
    
    # Find each region where the phase is changing too fast
    small = delta < min_feature
    region = {}
    j = 0
    for i in range(1, len(small)):
        # Beginning of a region
        if small[i] == True and small[i-1] == False:
            region[j] = {'start' : i+1}
        # End of a region
        if small[i] == False and small[i-1] == True:
            region[j]['end'] = i+1
            j += 1
        if i == len(small)-1 and j in region:
            region[j]['end'] = i+2
    
    for i in range(len(region)):
        start_ind = region[i]['start']
        end_ind = region[i]['end']
        start = rings[start_ind]['r']
        end = rings[end_ind]['r']
        # Remove the old rings
        for j in range(start_ind, end_ind):
            del rings[j]
        dr = end - start
        n = int(np.floor(dr / min_feature))
        dr = dr / n
        # Add in new rings
        for j in range(start_ind, start_ind+n):
            r_start = rings[j-1]['r']
            r = np.round(r_start+dr, 3)
            rf = np.round(r+dr, 3)
            sel = np.logical_and(r_nm >= r_start, r_nm < r)
            l = np.floor(np.average(phi_nm[sel])*layers/(2*np.pi))
            rings[j] = {'r':r, 'l':l, 'rf':rf}
    
    # Begin making the gds file
    print("Building GDS file")
    top = core.Cell("TOP")
    # You will never have 2**9 layers
    for i in range(10):
        if layers <= 2**i:
            Nlayers = i
    
    
    for key, item in rings.items():
        depth = layers-item['l']-1
        radius = item['rf']
        inner_radius = item['r']
        for k in range(Nlayers):
            etch = 2**(Nlayers-k-1)
            if (depth/etch) >= 1:
                ring = shapes.Disk((0, 0), radius, inner_radius=inner_radius, layer=etch, number_of_points=500)
                top.add(ring)
                depth = depth % etch
    
    # Add the top-cell to a layout and save
    layout = core.Layout("LAYOUT")
    layout.add(top)
    print("Beginning file saving")
    layout.save(fileName[:-4]+".gds")
    print("Finished with file saving")
    
if __name__ == "__main__":
    args = sys.argv
    if len(args) == 4:
        fileName = args[1] # .npy file with the array in it
        layers = int(args[2]) # Number of layers to discretize the phase into
        min_feature = float(args[3])
        main(fileName, layers, min_feature)
    else:
        print("Usage: python numpyArrayToGDS.py <fileName> <layers> <min_feature[um]>")