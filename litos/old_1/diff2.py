#!/usr/bin/python

#=======================
# calculate difference^2
#=======================

# common libraries
import sys
from numpy import *

def diff2(x,func):

    target = x[0]
    funcout = func(x[1:])
    diff = funcout-target
    return diff**2
