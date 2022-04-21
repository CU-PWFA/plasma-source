#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 20:33:52 2021

@author: valentinalee
"""

import numpy as np
import pandas as pd
import os
from os import path
import sys
import glob
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.realpath('../../03_DataAnalysis/04_Functions/'))
from AnalysisFn import Analysis
import multiprocessing
import functools
#%%
#TODO run all files
#use a for loop do calculate it for every files you have and also extract the pressure value for csv file


def do_computation(Path, df_A, args):
    i, name = args
    print(f'{multiprocessing.current_process()}: starting computation {i} ({name})')
    DataNumber= name[0:10]
    #run through all the files and find results and save them in a array
    PV_avg, Separation= Analysis(Path, DataNumber)
    #also read df to get the pressure
    OnlyN= DataNumber[6:10].lstrip('0')
    idx= int(df_A[df_A['Data Number']== OnlyN].index.values)
    Pressure= df_A['Gas Pressure Mix'][idx]
    CurrentResult= (Pressure, PV_avg, Separation)
    #CurrentResult = 0.1
    print(f'{multiprocessing.current_process()}: done with computation {i} ({name}): {CurrentResult}')
    return i, CurrentResult

if __name__ == '__main__':

    #load csv file
    print('\033[31mBLAH\033[m')
    CSVPath= os.path.realpath('../../03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/05_CSVFile/')
    df_A= pd.read_csv(CSVPath+ '/PressureScan.csv')
    df_A= df_A.dropna()
    df_A['Data Number']= df_A['Data Number'].astype('int').astype('str')

    #define the directory of clicky files
    Path= os.path.realpath('../../03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan')
    ClickyPath= os.path.realpath('../../03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/03_PeaksFoundByClickyTool')
    ResultPath= os.path.realpath('../../03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/04_Results')
    CheckList= []
    os.chdir(os.path.realpath(ClickyPath))


    data_files = list(glob.glob('??????????_PeaksInfo.npy'))
    print(data_files)
    results = [None] * len(data_files)
    processes = 8
    if processes > len(data_files):
        processes = len(data_files)
    partial = functools.partial(do_computation, Path, df_A)
    with multiprocessing.Pool(processes) as pool:
        for i, result in pool.imap_unordered(partial, enumerate(data_files)):
            results[i] = result

    ResultArray= np.array(results)
    np.save(ResultPath+ '/NormalizedByBG-angle', ResultArray)


'''

Result= []
#%%
for name in glob.glob('??????????_PeaksInfo.npy'):
    DataNumber= name[0:10]
    if (DataNumber not in CheckList):
        print(DataNumber)
#run through all the files and find results and save them in a array
        PV_avg, Separation= Analysis(Path, DataNumber)
        CheckList.append(DataNumber)
#also read df to get the pressure
        OnlyN= DataNumber[6:10].lstrip('0')
        idx= int(df_A[df_A['Data Number']== OnlyN].index.values)
        Pressure= df_A['Gas Pressure Mix'][idx]
        CurrentResult= (Pressure, PV_avg, Separation)
        Result.append(CurrentResult)

#%%
ResultArray= np.array(Result)
np.save(ResultPath+ '/NormalizedByBG-angle', ResultArray)
'''
#%%
'''
ResultArray= np.load('NormalizedByBG-angle.npy')
#%%
plt.title('Signal')
plt.plot(ResultArray[:, 0], ResultArray[:, 1], 'o')
plt.xscale('log')
plt.ylabel('(P-V)/BG')
plt.xlabel('Pressure (mbar)')
#%%
plt.title('Signal')
plt.plot(ResultArray[:, 0], ResultArray[:, 2], 'o')
plt.xscale('log')
plt.ylabel('Angle(rad)')
plt.xlabel('Pressure (mbar)')
'''
