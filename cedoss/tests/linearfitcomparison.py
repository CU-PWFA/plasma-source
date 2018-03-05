#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 10:55:48 2018

@author: chris
"""

import matplotlib.pyplot as plt
import numpy as np

m_arr = [-4.21e-1, -3.60e-1, -3.20e-1, -3.03e-1, -2.79e-1, -2.72e-1, -2.58e-1]
b_arr = [0.415e-2, 0.839e-2, 1.56e-2, 2.38e-2, 4.38e-2, 5.57e-2, 9.10e-2]
label_arr = [1,2.5,5,10,20,30,50]
x_arr=np.linspace(-1.0, 0, 300)

for i in range(len(m_arr)):
    plt.plot(x_arr, m_arr[i] * x_arr + b_arr[i], label = str(label_arr[i]))
plt.grid(); plt.legend(); plt.show()

plt.plot(x_arr, np.average(m_arr) * x_arr + np.average(b_arr), label = "average")
print("m_average:",np.average(m_arr))
print("b_average:",np.average(b_arr))
plt.grid(); plt.legend(); plt.show()