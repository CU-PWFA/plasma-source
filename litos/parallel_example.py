#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 17:14:03 2017

@author: mike
"""

import numpy as np
from joblib import Parallel, delayed
import multiprocessing
     
# what are your inputs, and what operation do you want to
# perform on each input. For example...
inputs = range(10)

def processInput(derp,i):
    print(derp[i])
    return [i,i * i]

derp = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(processInput)(derp,i) for i in inputs)

print(num_cores)
print(results)
print(np.reshape(results,[10,2])[2,1])