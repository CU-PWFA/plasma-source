# -*- coding: utf-8 -*-
"""
Created on Wed May 31 09:45:48 2017

@author: chris
"""
import numpy as np
I=np.zeros([3,3,3])
i=1
while i<4:
    j=1
    while j<4:
        k=1
        while k<4:
            I[i-1][j-1][k-1]=i*100+j*10+k
            k=k+1
        j=j+1
    i=i+1
    
h=I[:,0,:]
g=np.transpose(h)