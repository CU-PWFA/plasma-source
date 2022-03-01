#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  8 15:00:59 2022

@author: valentinalee
"""

import numpy as np
import pandas as pd
import os
from os import path
import sys
import glob
import matplotlib.pyplot as plt

#%%
#load csv file
CSVPath= os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/05_CSVFile/')
df_A= pd.read_csv(CSVPath+ 'PressureScan.csv')
df_A= df_A.dropna()
df_A['Data Number']= df_A['Data Number'].astype('int').astype('str')

#define result path
ResultPath= os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/04_Results')

FinalScanResult= []
#laod ResultList and ab
#calculate PV avg and angle
#append final result
for name in glob.glob(ResultPath+ '/??????????_ab.npy'):
    DataNumber= name[-17: -7]
    print(DataNumber)
    ResultList= np.load(ResultPath+ '/'+ DataNumber+ '_ResultList.npy')
    PV= ((ResultList[:, 1]+ResultList[:, 2])/2-ResultList[:, 3])/ResultList[:, 5]
    PV_avg= sum(PV)/ len(PV)
    ab= np.load(ResultPath+ '/'+ DataNumber+ '_ab.npy')
    Angle= np.arctan((ab[0]-ab[2])/(1+ab[0]*ab[2]))
    
    OnlyN= DataNumber[6:10].lstrip('0')
    idx= int(df_A[df_A['Data Number']== OnlyN].index.values)
    Pressure= df_A['Gas Pressure Mix'][idx]
    
    CurrentResult= (Pressure, PV_avg, Angle)
    FinalScanResult.append(CurrentResult)
    
FinalScanResultArray= np.array(FinalScanResult)
#%%

#np.save('NormalizedByBG-angle', ResultArray)

#%%
plt.figure()
plt.title('Signal')
plt.plot(FinalScanResultArray[:, 0], FinalScanResultArray[:, 1], 'o')
plt.xscale('log')
plt.ylabel('(P-V)/BG')
plt.xlabel('Pressure (mbar)')
#%%
plt.figure()
plt.title('Signal')
plt.plot(FinalScanResultArray[:, 0], FinalScanResultArray[:, 2], 'o')
plt.xscale('log')
plt.ylabel('Angle(rad)')
plt.xlabel('Pressure (mbar)')
