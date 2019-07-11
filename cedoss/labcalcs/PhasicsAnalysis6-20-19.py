#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 11:58:42 2019

Analysis of Phasics data from June 20, 2019.  This one is more sophisticated because
I take the background and average, then for each of the cases I average and
subtract the background.  The cases are then compared for qualitiative presenting

Also has some pretty cool analysis functions defined below.

Need to find some theory to compare results to, probably also quantify the
imaging system better.

@author: chris
"""

#from PIL import Image as im
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join

folder = '/home/chris/Desktop/LabData/Phasics/June202019Data/'

def getAverageSet(directory, case):
    setlist = [f for f in listdir(directory+case) if isfile(join(directory+case, f))]
    image = None
    for set_i in setlist:
        imagefile = directory + case + set_i
        if image == None:
            image = np.array(plt.imread(imagefile),dtype='double')
        else:
            image = np.array(plt.imread(imagefile),dtype='double') + image
    return image / len(setlist)


def averageSlices(image, centerLine, noSamples = 10, sampleFreq = 1):
    dataline = image[:,centerLine]
    samples = np.zeros(2*noSamples + 1)
    for i in range(len(dataline)):
        samples[0] = dataline[i]
        for jj in range(noSamples):
            j = jj + 1
            samples[2*j-1]=image[i,centerLine - sampleFreq*j]
            samples[2*j] = image[i,centerLine + sampleFreq*j]
        dataline[i] = np.average(samples)
    return dataline

background = getAverageSet(folder, '1mbar_background/')

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(background)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("Averaged Background Probe Intensity")
plt.show()

image_1mbar = getAverageSet(folder, '1mbar/')
image_2mbar = getAverageSet(folder, '2mbar/')
image_4mbar = getAverageSet(folder, '4mbar/')
image_8mbar = getAverageSet(folder, '8mbar/')
image_12mbar = getAverageSet(folder, '12mbar/')
image_16mbar = getAverageSet(folder, '16mbar/')
image_20mbar = getAverageSet(folder, '20mbar/')

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(image_20mbar-background)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("20mbar Ar Axicon Plasma - Averaged Probe Intensity Background Subtracted")
plt.show()

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(image_4mbar-background)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("4mbar Ar Axicon Plasma - Averaged Probe Intensity Background Subtracted")
plt.show()

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(image_1mbar-background)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("1mbar Ar Axicon Plasma - Averaged Probe Intensity Background Subtracted")
plt.show()

left20mbar = (image_20mbar-background)[:,30]

plt.figure(figsize=(20,10))
plt.plot(left20mbar)
plt.title("Single Slice of left side of 20mbar")
plt.show()

smoothed_left20mbar = averageSlices(image_20mbar - background, 40, 30, 1)
smoothed_right20mbar = averageSlices(image_20mbar - background, 600, 30, 1)
plt.figure(figsize=(20,10))
plt.plot(smoothed_left20mbar, label="Left of focus")
plt.plot(smoothed_right20mbar,label="Right of foucs")
plt.title("Smoothing achieved by averaging 61 slices")
plt.legend(); plt.show()

smoothed_left16mbar = averageSlices(image_16mbar - background, 40, 30, 1)
smoothed_right16mbar = averageSlices(image_16mbar - background, 600, 30, 1)
smoothed_left12mbar = averageSlices(image_12mbar - background, 40, 30, 1)
smoothed_right12mbar = averageSlices(image_12mbar - background, 600, 30, 1)
smoothed_left8mbar = averageSlices(image_8mbar - background, 40, 30, 1)
smoothed_right8mbar = averageSlices(image_8mbar - background, 600, 30, 1)
smoothed_left4mbar = averageSlices(image_4mbar - background, 40, 30, 1)
smoothed_right4mbar = averageSlices(image_4mbar - background, 600, 30, 1)
smoothed_left2mbar = averageSlices(image_2mbar - background, 40, 30, 1)
smoothed_right2mbar = averageSlices(image_2mbar - background, 600, 30, 1)
smoothed_left1mbar = averageSlices(image_1mbar - background, 40, 30, 1)
smoothed_right1mbar = averageSlices(image_1mbar - background, 600, 30, 1)

plt.figure(figsize=(20,10))
plt.plot(smoothed_left20mbar, label="20mbar")
plt.plot(smoothed_left16mbar, label="16mbar")
plt.plot(smoothed_left12mbar, label="12mbar")
plt.plot(smoothed_left8mbar, label="8mbar")
plt.plot(smoothed_left4mbar, label="4mbar")
plt.plot(smoothed_left2mbar, label="2mbar")
plt.plot(smoothed_left1mbar, label="1mbar")
plt.title("Diffraction vs Pressure - Left of Center")
plt.legend(); plt.show()

plt.figure(figsize=(20,10))
plt.plot(smoothed_right20mbar, label="20mbar")
plt.plot(smoothed_right16mbar, label="16mbar")
plt.plot(smoothed_right12mbar, label="12mbar")
plt.plot(smoothed_right8mbar, label="8mbar")
plt.plot(smoothed_right4mbar, label="4mbar")
plt.plot(smoothed_right2mbar, label="2mbar")
plt.plot(smoothed_right1mbar, label="1mbar")
plt.title("Diffraction vs Pressure - Right of Center")
plt.legend(); plt.show()

## Now we shift focus to the delay stage variation

midbackground = getAverageSet(folder, '20mbar_Middle_Background/')
midimage_20mbar = getAverageSet(folder, '20mbar_Middle/')
frontbackground = getAverageSet(folder, '20mbar_Farforward_Background/')
frontimage_20mbar = getAverageSet(folder, '20mbar_Farforward/')

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(midimage_20mbar - midbackground)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("20mbar Ar Axicon Plasma - Average Intensity Background Subtracted - Delay Stage at Middle")
plt.show()

plt.figure(figsize=(20,10))
plt.set_cmap('gray')
plt.imshow(frontimage_20mbar - frontbackground)
CB = plt.colorbar()
#plt.clim(0,1.)
plt.title("20mbar Ar Axicon Plasma - Average Intensity Background Subtracted - Delay Stage at Front")
plt.show()

smoothed_midleft20mbar = averageSlices(midimage_20mbar - midbackground, 40, 30, 1)
smoothed_midright20mbar = averageSlices(midimage_20mbar - midbackground, 600, 30, 1)
smoothed_frontleft20mbar = averageSlices(frontimage_20mbar - frontbackground, 40, 30, 1)
smoothed_frontright20mbar = averageSlices(frontimage_20mbar - frontbackground, 600, 30, 1)

xaxis = np.arange(0,len(smoothed_midleft20mbar),1)

plt.figure(figsize=(20,10))
plt.plot(xaxis+45, smoothed_left20mbar+35, label="Far Back")
plt.plot(xaxis, smoothed_midleft20mbar+10, label="Middle")
plt.plot(xaxis+18, smoothed_frontleft20mbar, label="Far Front")
plt.title("Diffraction vs Delay Stage - Left of Center")
plt.legend(); plt.show()

plt.figure(figsize=(20,10))
plt.plot(xaxis+46, smoothed_right20mbar+4, label="Far Back")
plt.plot(xaxis, smoothed_midright20mbar+15, label="Middle")
plt.plot(xaxis+18, smoothed_frontright20mbar, label="Far Front")
plt.title("Diffraction vs Delay Stage - Right of Center")
plt.legend(); plt.show()