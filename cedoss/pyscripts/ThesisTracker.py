#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 13:45:50 2022

Plotting my bullshit thesis progress

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
import math


pagecount = np.array([3, 3, 12, 30, 35, 64, 69, 80, 95, 102, 113, 124, 134, 144, 178, 200])
#day count starts at day 1 at Dec 1 2022
daycount = np.array([1, 8, 15, 
                     31+18, 31+24, 
                     62+7, 62+14, 62+21,
                     90+21, 90+28,
                     121+4, 121+18, 121+25,
                     151+9, 151+16, 151+22])

#end of dec, end of jan ..., end of april
monthcount = np.array([0,31,62,90,121,151,182])
monthlabel = np.array(['Dec','Jan','Feb','Mar','Apr','May',''])

pagegoal = 150

xpad = 10
ypad = 10+55
extend = 30

dayrange = np.array([0,152])
dayaxis = np.arange(dayrange[0]-10,dayrange[-1]+10+extend)

#linfit = np.polyfit(daycount,pagecount,1)
linfit = np.polyfit(daycount[1:],pagecount[1:],1)

expected = (pagegoal - linfit[1])/linfit[0] - daycount[-1]

plt.plot([dayrange[0]-xpad,dayrange[-1]+xpad+extend],[pagegoal,pagegoal],ls = 'dashed',c='r',label='Goal = '+str(pagegoal))
plt.plot([dayrange[0]-xpad,dayrange[-1]+xpad+extend],[0,0],c='k')
for i in range(len(monthcount)):
    plt.plot([monthcount[i],monthcount[i]],[-5,pagegoal+ypad],ls = 'dotted',c='b')
    plt.text(monthcount[i]+10,pagegoal-20,monthlabel[i])
plt.plot([182+5,182+5],[-5,pagegoal+ypad],ls = 'dotted',c='black',label="June 5th Defense")
plt.plot([151+22,151+22],[-5,pagegoal+ypad],ls = 'dotted',c='orange',label="May 22nd Final Draft")
plt.scatter(daycount,pagecount,c='g',label='Current = '+str(pagecount[-1]))
plt.plot(dayaxis,dayaxis*linfit[0]+linfit[1],c='g')#,label='Exp. Days Remaining: '+str(math.ceil(expected)))
plt.ylim([-5,pagegoal+ypad])
plt.xlim([dayrange[0]-xpad,dayrange[-1]+xpad+extend])
plt.legend(loc=4);plt.show()























