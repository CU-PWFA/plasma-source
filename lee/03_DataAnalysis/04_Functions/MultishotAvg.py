#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 10:11:50 2021

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from photutils.centroids import centroid_com
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
from scipy import optimize
import os
import glob
from scipy.ndimage.interpolation import shift


def MultishotAvg(DataDir, BGDir, ResultDir, shotN= 50, filterSigma= 2):
    '''
    DataDir: Directory of the data data
    BGDir: Directory of the background data
    ResultDir: Directory where the average signal, data, background (npy) result will be saved in /01_npy and 
                the average signal plot will be saved in /02_Plots
    shotN: shot number, default 50 shots
    filterSignam: the std of the Gaussian filter applied on the average signal 
    '''
    
    os.chdir(os.path.expanduser(BGDir))
    for name in glob.glob('17529184_??????????_????.tiff'):
        DataNumber= name[9:19]    
    print(DataNumber)
    
    for shot in range (0, shotN):
        strN= str(shot)
        FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
        print(FileName)
        if shot==0:
            BG= np.array(Image.open(FileName))
            TotalBG= np.zeros((BG.shape[0], BG.shape[1]))
            BG_CX0, BG_CY0= centroid_com(BG)
        
        BG= np.array(Image.open(FileName))    
        BG_cx, BG_cy= centroid_com(BG)
    #move BG so that centroid goes to center
        BG= shift(BG, ((BG_CY0-BG_cy), (BG_CX0-BG_cx)))
        TotalBG= TotalBG+ BG
    
    AvgBG= TotalBG/shotN
    
    def shiftBG(s1, s2):
        ShiftedBG= shift(AvgBG, (s1, s2))
        return ShiftedBG
    
    os.chdir(os.path.expanduser(DataDir))
    for name in glob.glob('17529184_??????????_????.tiff'):
        DataNumber= name[9:19]    
    print(DataNumber)

    for shot in range (0, shotN):
        if shot== 0:
            strN= str(shot)
            FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
            
        #load the tiff image and change it to array
            data= np.array(Image.open(FileName))
            errorFn= lambda p1: np.ravel(shiftBG(*p1)- data)
            p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
            data= shift(data, (-p[0], -p[1]))
            TotalData= np.zeros((data.shape[0], data.shape[1]))
            signal= data-AvgBG
            signalFiltered= gaussian_filter(signal, sigma=6)

            sumSignal= np.sum(-signalFiltered, axis= 1)[200:400]
            peaks, properties = find_peaks(sumSignal, prominence= \
                                (np.amax(np.sum(-signalFiltered, axis= 1))-np.amin(np.sum(-signalFiltered, axis= 1)))*0.015, \
                                width= (10, 100))
            
            if len(peaks)>1:
                PeakWidth= properties["widths"]
                sortedArray= sorted(range(len(PeakWidth)), key=lambda k: PeakWidth[k])
                Peak= int(peaks[sortedArray[0]]+200)

                #                PeakInt= [np.sum(-signalFiltered, axis= 1)[int(200+peaks[n])] for n in range (0, len(peaks))]
#                sortedArray= sorted(range(len(PeakInt)), key=lambda k: PeakInt[k])
#                Peak= int(peaks[sortedArray[-1]]+200)
            elif len(peaks)==1:
                Peak0= peaks[0]+200
            else:
                print('No peak found in '+FileName)
                continue
            
        strN= str(shot)
        FileName= '17529184_'+DataNumber+'_00'+strN.zfill(2)+'.tiff'
        print(FileName)
        
    #load the tiff image and change it to array
        data= np.array(Image.open(FileName))
        errorFn= lambda p1: np.ravel(shiftBG(*p1)- data)
        p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
        print('Moved Coords:', p)
        data= shift(data, (-p[0], -p[1]))
        signal= data-AvgBG
        signalFiltered= gaussian_filter(signal, sigma=6)
        sumSignal= np.sum(-signalFiltered, axis= 1)[200:400]
        peaks, properties = find_peaks(sumSignal, prominence= \
                                       (np.amax(np.sum(-signalFiltered, axis= 1))-np.amin(np.sum(-signalFiltered, axis= 1)))*0.015, \
                                       width= (10, 100))

        if 'Peak0' in locals():
            if len(peaks)>1:
                print(peaks)
                DisToSeed= [abs((p+200)-Peak0) for p in peaks]
                print(DisToSeed)
                sortedArray= sorted(range(len(DisToSeed)), key=lambda k: DisToSeed[k])
                Peak= int(peaks[sortedArray[0]]+200)
            elif len(peaks)==1:
                Peak= peaks[0]+200
            else:
                print('No peaks found in '+FileName)
                continue
            movedData= shift(data, (Peak0-Peak, 0))
            TotalData= TotalData+movedData            

        else:
            if len(peaks)>1:
                PeakWidth= properties["widths"]
                sortedArray= sorted(range(len(PeakWidth)), key=lambda k: PeakWidth[k])
                Peak= int(peaks[sortedArray[0]]+200)
    
    #            PeakInt= [np.sum(-signalFiltered, axis= 1)[int(200+peaks[n])] for n in range (0, len(peaks))]
    #            sortedArray= sorted(range(len(PeakInt)), key=lambda k: PeakInt[k])
    #            Peak= int(peaks[sortedArray[-1]]+200)
            elif len(peaks)==1:
                Peak= peaks[0]+200
            else:
                print('No peaks found in '+FileName)
                continue
            Peak0= Peak

        print(Peak)
        plt.figure(1)
        plt.plot(np.sum(-signalFiltered, axis= 1))
        plt.plot(Peak, np.sum(-signalFiltered, axis= 1)[Peak], 'x')

    try:
        AvgData= TotalData/shotN
        errorFn= lambda p1: np.ravel(shiftBG(*p1)- AvgData)
        p, success= optimize.leastsq(errorFn, np.array([20, 20]), epsfcn=2e-4)
        print(p)
        AvgData= shift(AvgData, (-p[0], -p[1]))
        signal= AvgData-AvgBG
        signalFiltered= gaussian_filter(signal, sigma=2)
        np.save(os.path.join(os.path.expanduser(ResultDir), '01_npy', DataNumber+'_NorSignal'), \
                (signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)))
        np.save(os.path.join(os.path.expanduser(ResultDir), '01_npy', DataNumber+'_Signal'), signal)
        np.save(os.path.join(os.path.expanduser(ResultDir), '01_npy', DataNumber+'_Data'), AvgData)
        np.save(os.path.join(os.path.expanduser(ResultDir), '01_npy', DataNumber+'_BG'), AvgBG)
        
        plt.figure()
        plt.pcolormesh((signalFiltered-np.amin(signalFiltered))/np.amax(signalFiltered-np.amin(signalFiltered)), \
                       cmap= 'bwr')
        plt.axis('scaled')
        plt.colorbar()
        plt.title(DataNumber+'_NorSignal')
        plt.savefig(os.path.expanduser(ResultDir)+'/02_Plots/'+DataNumber+'_NorSignal')
        plt.close()
        
        plt.figure()
        plt.pcolormesh(signal, cmap= 'PiYG')
        plt.axis('scaled')
        plt.colorbar()
        plt.title(DataNumber+'_Signal')
        plt.savefig(os.path.expanduser(ResultDir)+'/02_Plots/'+DataNumber+'_Signal')
        plt.close()
        
        plt.figure(1)
        plt.savefig(os.path.expanduser(ResultDir)+'/02_Plots/'+DataNumber+'_ManualCheck')
        plt.close()
        
        np.save()
    except:
        print('No Peaks found in Dataset number '+ DataNumber)
        