#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 13:22:22 2019

@author: robert
"""

import sys
import numpy as np
from gdsCAD import *


def main(fileName, cellSize):
    """ Convert a numpy array to a GDS file with a grid of squares.
    
    All identical elements of an array are placed on their own layer. This code
    assumes all array elements are multiples of some base element. It is 
    meant to be used on an array of phase for each pixel. Higher phase
    therefore corresponds to a shallower etch. 
    
    Note: one day this should be switched to gdspy which is actively maintained
    and works with python 3. It also doesn't require modifying the source code
    in order to get it to work.
    
    Parameters
    ----------
    fileName : string
        name of the file to convert, will also be the name of the output file.
    cellSize : float
        The size of each square in um, precision down to nm is included in gds.
    """
    print("Converting the array in "+fileName+" to a GDS file")
    data = np.load(fileName)
    print("Successfully loaded "+fileName)
    try:
        height, width = np.shape(data)
    except:
        print("The numpy file doesn't have the correct format, it must be a single array")
        return
    print("Height: %d" % height)
    print("Width: %d" % width)
    elements = np.unique(data)
    base = elements[1]
    mask = elements[-1]-data
    mask = mask/base
    Nlayers = len(elements)/2
    print("Layers: %d" % Nlayers)
    print("Layers should be half the number of discrete depths")

    # Create the grid and dictionary for storing the unit cells
    unitCells = {}
    grid = core.Cell("GRID")
    for i in range(Nlayers):
        unitCells[i+1] = core.Cell("Layer%d" % (i+1))
        square = shapes.Rectangle((0.0, 0.0), (1.0, 1.0), layer=(int)(i+1))
        unitCells[i+1].add(square)
    print("Unit cells created")
    
    print("Begining layer construction")
    for i in range(height):
        for j in range(width):
            depth = mask[i, j]
            for k in range(Nlayers):
                etch = Nlayers-k
                if (depth/etch) >= 1:
                    cell = core.CellReference(unitCells[etch], origin=(j, height-i-1))
                    grid.add(cell)
                    depth = depth % etch
    print("Finished with layer construction")

    scaledGrid = core.CellReference(grid, origin=(0, 0), magnification=(float)(cellSize))
    
    top = core.Cell("TOP")
    top.add(scaledGrid)
    
    # Add the top-cell to a layout and save
    layout = core.Layout("LAYOUT")
    layout.add(top)
    print("Beginning file saving")
    layout.save(fileName[:-4]+".gds")
    print("Finished with file saving")

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 3:
        fileName = args[1] # .npy file with the array in it
        cellSize = args[2] # Size of a single pixel in the array (um), XX.XXX
        main(fileName, cellSize)
    else:
        print("Usage: python numpyArrayToGDS.py <fileName> <cellSize[um]>")
        