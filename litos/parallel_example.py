#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 17:14:03 2017

@author: mike
"""

from joblib import Parallel, delayed
import multiprocessing
     
# what are your inputs, and what operation do you want to
# perform on each input. For example...
inputs = range(10)

def processInput(i):
    return i * i
 
num_cores = multiprocessing.cpu_count()
results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

print(num_cores)
print(results)