#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:58:20 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#%%
c= 3e8
#%%
df= pd.read_excel(r'/media/valentinalee/Transcend/_CUPWFAWork/Data/FrogData072621/FrogData072621.ods')
record= df.to_numpy()
record= sorted(record, key=lambda x: x[1])
#%%
Count= []
Delay= []
for data_number in range(len(record)):
    if data_number < len(record)-1:
        a= []
        with open('/media/valentinalee/Transcend/_CUPWFAWork/Data/FrogData072621/michaelito1_'+record[int(data_number)][0]+'.txt',"r") as r:
            a=r.readlines()
        a= [x.replace("\n","") for x in a]
        a= a[a.index(">>>>>Begin Spectral Data<<<<<")+1:]
        
        lam= []
        count= []
        for line in a:
            lam.append(float(line.split()[0]))
            count.append(float(line.split()[1]))
        
        lam= np.array(lam)
        count= np.array(count)
        
        Count.append(count)
        Delay.append(record[int(data_number)][1])

    
Count= np.array(Count)
Delay= np.array(Delay)
#%%
plt.pcolormesh(lam[625: 875], Delay*1e-3/c*1e12, Count[:, 625: 875])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Delay (ps)')