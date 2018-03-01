#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:26:53 2018

@author: chris
"""

import matplotlib.pyplot as plt
import numpy as np

kb_arr = [94.08, 210.37, 297.51]
sigm_arr= [1.48,1.72,1.71]
sigb_arr= [.0209,.0236,.0192]

za_arr = [-3.42, -3.77, -3.92]
zb_arr = [-3.30, -4.17, -4.67]
zc_arr = [0.0678, 0.0506, 0.0670]

#First is sigma vs beta
x_arr=np.linspace(kb_arr[0], kb_arr[-1], 200)

plt.scatter(kb_arr,sigm_arr,label="sig's m vs kb")
plt.grid(); plt.legend(); plt.show()

plt.scatter(kb_arr,sigb_arr, label="sig's b vs kb")
plt.grid(); plt.legend(); plt.show()

plt.scatter(kb_arr,za_arr, label="z's a vs kb")
plt.grid(); plt.legend(); plt.show()

plt.scatter(kb_arr,zb_arr, label="z's b vs kb")
plt.grid(); plt.legend(); plt.show()

plt.scatter(kb_arr,zc_arr, label="z's c vs kb")
plt.grid(); plt.legend(); plt.show()