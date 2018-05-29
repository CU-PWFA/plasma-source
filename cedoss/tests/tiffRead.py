#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:32:37 2018

Test loading .tiff images into python

@author: chris
"""

import matplotlib.pyplot as plt

folder = '/home/chris/Desktop/day_21'

imagefile = folder + '1805210001/Cam02_1805210001_0001.tiff'
image = plt.imread(imagefile)