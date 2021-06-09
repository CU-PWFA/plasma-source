#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:51:59 2020

@author: valentinalee
"""
#%%
from os import path

for density in range (1, 11, 1):
    
    for power in range (1, 11):
        
        for pw in range (30, 160, 10):
            Dir= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'/'
            File= 'cap_2D.pre'
            NewPreFile= Dir+File
            if path.exists(Dir)==True:
                reading_file = open(NewPreFile, "r")
                for line in reading_file:
                    stripped_line = line.strip()
                    if (stripped_line.find('lowerBounds =') ==0):
                        Box= int(stripped_line[16:24])
                    if (stripped_line.find('numZones =') ==0):
                        GridSize= int(stripped_line[16:19])
                    if (stripped_line.find('$ SIM_TIME =') ==0):
                        SimTime=float(stripped_line[13:18])
                    if (stripped_line.find('$ N_FRAME') ==0):
                        Frame=int(stripped_line[12:14])
                print(Dir)
                print(Box)
