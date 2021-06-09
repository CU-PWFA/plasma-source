#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:00:10 2020

@author: valentinalee
"""

#Run this in the directry where all the n_GP_PW files sit
#%%
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from os import path

#%%
ELECMASS = 9.10938356e-31     # kG
PROTMASS = 1.672621898e-27    # kG
n0=1e23 
kb=1.38064852e-23  # J/K 
#%% Call data
def qReader(FileName):
    hf=h5py.File(FileName,'r')
    Flu= hf['fluids']
    q=np.array(Flu.get('q'))
    Den=q[:, 0]
    Eng=q[:, 4]
    return Den, Eng

def q_meshReader(FileName):
    hf= h5py.File(FileName,'r')
    Flu= hf['fluids']['domain']['face']
    cells= np.array(Flu.get('cells'))
    vertices= np.array(Flu.get('vertices'))
    return cells, vertices

#%%
SimExt= np.zeros(4)
for density in range (1, 101, 1):
    
    for power in range (1, 11):
        
        for pw in range (30, 160, 5):

            CheckRun= 'ScanResults/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'.npy'
            if path.exists(CheckRun) == False:

                Dir= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'/'
                File= 'cap_2D.pre'
                NewPreFile= Dir+File
                if path.exists(Dir)==True:
                    reading_file = open(NewPreFile, "r")
                    for line in reading_file:
                        stripped_line = line.strip()
                        if (stripped_line.find('lowerBounds =') ==0):
                            Box= float(stripped_line[16:24])
                        if (stripped_line.find('numZones =') ==0):
                            GridSize= int(stripped_line[16:19])
                        if (stripped_line.find('$ SIM_TIME =') ==0):
                            SimTime=float(stripped_line[13:19])
                        if (stripped_line.find('$ N_FRAME') ==0):
                            Frame=int(stripped_line[12:15])
                    print(Dir)
                    print(Box, GridSize, SimTime, Frame)
                    Box= Box*1e6
                    SimExt[0]=Box #in micron
                    SimExt[1]=GridSize
                    SimExt[2]=SimTime #in s
                    SimExt[3]=Frame
                    np.save('ScanResults/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_SimExt', SimExt)
                    BaseFileName= 'cap_2D_'
                    
                    grid_x= np.linspace(-Box*1e-6, Box*1e-6, GridSize)
                    grid_y= np.linspace(-Box*1e-6, Box*1e-6, GridSize)
                    grid_X, grid_Y= np.meshgrid(grid_x, grid_y)
                    
                    den_3d=np.zeros((GridSize, GridSize, Frame))
    
                    BasePath= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'/'
                    Gridpath= BasePath+ 'cap_2DGrid.h5'
                    cells, vertices= q_meshReader(Gridpath)
        
                    for n in range (0, Frame):
                        FileName= BaseFileName+str(n)+'.h5'
                        Path= BasePath+ FileName
                        print(Path)
                        den= qReader(Path) [0]/(ELECMASS+PROTMASS*2*2)/(100)**3
                        grid=np.zeros((len(den), 2))
        
                        for cell in range (0, len(den)):
                            corner1x= vertices[cells[cell, 0], 0]
                            corner1y= vertices[cells[cell, 0], 1]
                            corner2x= vertices[cells[cell, 1], 0]
                            corner2y= vertices[cells[cell, 1], 1]
                            corner3x= vertices[cells[cell, 2], 0]
                            corner3y= vertices[cells[cell, 2], 1]
                            corner4x= vertices[cells[cell, 3], 0]
                            corner4y= vertices[cells[cell, 3], 1]
                
                            grid[cell, 0]= (corner1x+corner2x+corner3x+corner4x)/4
                            grid[cell, 1]= (corner1y+corner2y+corner3y+corner4y)/4
        
                        den_3d[:, :, n] = griddata(grid, den, (grid_X, grid_Y), method='linear')     
                    LineOutArray= np.zeros((GridSize, Frame))
                    for n in range (0, Frame):
                        LineOutArray[:, n]= den_3d[:, :, n][:, int(GridSize/2)]
                    np.save('ScanResults/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw), LineOutArray)
                    time= np.linspace(0, SimTime*1e9, Frame)
                    r= np.linspace(-Box, Box, GridSize)
                    plt.pcolormesh(time, r, LineOutArray)
                    plt.colorbar()
                    plt.ylabel('LineOut $\mu$m')
                    plt.xlabel('Time (ns)')
                    plt.title('n'+str(density)+'_GP'+str(power)+'_PW'+str(pw))
                    plt.savefig('ExpansionPlots/n'+str(density)+'_GP'+str(power)+'_PW'+str(pw))
                    plt.close()
    
    
    
