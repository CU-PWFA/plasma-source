#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 10:44:47 2021

@author: valentinalee
"""

#%%
import pandas as pd
import os
from os import path
import sys

sys.path.insert(0, os.path.expanduser('~/Desktop/Guerrero DAQ local copy/plasma-source-master/lee/03_DataAnalysis/04_Functions/'))
from MultishotAvg import MultishotAvg

CSVPath= os.path.expanduser('~/Desktop/Guerrero DAQ local copy/plasma-source-master/lee/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/05_CSVFile/')
df_A= pd.read_csv(CSVPath+ 'PressureScan.csv')

df_A= df_A.dropna()
df_A['Data Number']= df_A['Data Number'].astype('int').astype('str')

BasePath= '~/Desktop/Guerrero DAQ local copy/DAQ_Data/mnt/md0/DAQ/IMAGE/year_2021/month_10/day_29/' #DataDir

BaseResultDir= '~/Desktop/Guerrero DAQ local copy/plasma-source-master/lee/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/'
for idx in range(0, df_A.shape[0]-1):

    if df_A['Gas Pressure Mix'][idx]== df_A['Gas Pressure Mix'][idx+1]:
        if df_A['Background or not'].str[0][idx]=='D':
            yy= df_A['Date '].str[6:8][idx]
            dd= df_A['Date '].str[3:5][idx]
            mm= df_A['Date '].str[0:2][idx]
            DataNum= df_A['Data Number'].str[:][idx]
            DataDir= BasePath+yy+mm+dd+DataNum.zfill(4)+'/'
            BGDataNum= df_A['Data Number'].str[:][idx+1]
            BGDir= BasePath+yy+mm+dd+BGDataNum.zfill(4)+'/'
        else:
            yy= df_A['Date '].str[6:8][idx]
            dd= df_A['Date '].str[3:5][idx]
            mm= df_A['Date '].str[0:2][idx]
            DataNum= df_A['Data Number'].str[:][idx+1]
            DataDir= BasePath+yy+mm+dd+DataNum.zfill(4)+'/'
            BGDataNum= df_A['Data Number'].str[:][idx]
            BGDir= BasePath+yy+mm+dd+BGDataNum.zfill(4)+'/'
        print(DataDir, BGDir)
        CheckRun= os.path.expanduser(BaseResultDir+'01_npy/')+ yy+ mm+ dd+ DataNum.zfill(4)+ '_Avg.npy'
        if path.exists(CheckRun) == False:
            MultishotAvg(DataDir, BGDir, BaseResultDir, shotN= 50)
#%%
#BaseResultDir= '~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/'
#MultishotAvg(DataDir, BGDir, BaseResultDir+'01_npy/', BaseResultDir+'02_Plots/', shotN= 5)
