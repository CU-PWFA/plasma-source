#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:52:50 2017

@author: litos
"""

import numpy as np
from calc_rms import calc_rms


derp = [1,1,1,1,1,1,1]

deep = calc_rms(derp)

print(deep)