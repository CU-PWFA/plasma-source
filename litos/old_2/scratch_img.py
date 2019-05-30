#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 15:08:56 2018

@author: mike
"""
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

im = Image.open('/home/mike/Desktop/Cam01_1804130002_0001.tiff')
#im.show()
imarray = np.array(im)
print(np.max(imarray))

img = plt.imread('/home/mike/Desktop/Cam01_1804130002_0001.tiff')