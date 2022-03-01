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
file_name= 'Gratings2000.ods'
df= pd.read_excel(r'/media/valentinalee/Transcend/MainFrogData081121/'+file_name)
record= df.to_numpy()
record= sorted(record, key=lambda x: x[1])
#%%
Count= []
Delay= []
date= '081121'
for data_number in range(len(record)):
    if data_number < len(record)-1:
        data= np.load('/media/valentinalee/Transcend/MainFrogData'+date+'/'+\
                      str(int(record[int(data_number)][0]))+'/michaelito_'+\
                      str(int(record[int(data_number)][0]))+'_0001.npy', allow_pickle=True)
        
        lam= data.item()['lambda']
        count= data.item()['I']
                
        Count.append(count)
        Delay.append(record[int(data_number)][1])

    
Count= np.array(Count)
Delay= np.array(Delay)
#%%
plt.pcolormesh(lam[625: 875], Delay*2*1e-3/c*1e12, Count[:, 625: 875])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Delay (ps)')
plt.title('Grating Location= 2000')
plt.colorbar()