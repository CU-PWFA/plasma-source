#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:39:43 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
im01 = Image.open('01/17583372_1911050001_0002.tiff')
im02 = Image.open('02/17583372_1911050002_0009.tiff')
im03 = Image.open('03/17583372_1911050003_0002.tiff')
im04 = Image.open('04/17583372_1911050004_0007.tiff')
im05 = Image.open('05/17583372_1911050005_0002.tiff')
im06 = Image.open('06/17583372_1911050006_0005.tiff')
im07 = Image.open('07/17583372_1911050007_0005.tiff')
im08 = Image.open('08/17583372_1911050008_0008.tiff')
im09 = Image.open('09/17583372_1911050009_0009.tiff')
im13 = Image.open('13/17583372_1911050013_0004.tiff')
im14 = Image.open('14/17583372_1911050014_0001.tiff')
im15 = Image.open('15/17583372_1911050015_0007.tiff')
im16 = Image.open('16/17583372_1911050016_0008.tiff')
im17 = Image.open('17/17583372_1911050017_0006.tiff')
im18 = Image.open('18/17583372_1911050018_0006.tiff')
im19 = Image.open('19/17583372_1911050019_0007.tiff')
im20 = Image.open('20/17583372_1911050020_0006.tiff')
im21 = Image.open('21/17583372_1911050021_0006.tiff')
im22 = Image.open('22/17583372_1911050022_0007.tiff')
im23 = Image.open('23/17583372_1911050023_0002.tiff')
im24 = Image.open('24/17583372_1911050024_0007.tiff')
im25 = Image.open('25/17583372_1911050025_0008.tiff')
im26 = Image.open('26/17583372_1911050026_0008.tiff')
im27 = Image.open('27/17583372_1911050027_0007.tiff')
im29 = Image.open('29/17583372_1911050029_0007.tiff')
im30 = Image.open('30/17583372_1911050030_0007.tiff')
im31 = Image.open('31/17583372_1911050031_0006.tiff')
im32 = Image.open('32/17583372_1911050032_0008.tiff')
im33 = Image.open('33/17583372_1911050033_0007.tiff')
im36 = Image.open('36/17583372_1911050036_0008.tiff')
im38 = Image.open('38/17583372_1911050038_0004.tiff')
im39 = Image.open('39/17583372_1911050039_0002.tiff')
im40 = Image.open('40/17583372_1911050040_0007.tiff')
im41 = Image.open('41/17583372_1911050041_0003.tiff')
#%%
PickOff={}
PickOff[1]= np.array(im01)
PickOff[2]= np.array(im02)
PickOff[3]= np.array(im03)
PickOff[4]= np.array(im04)
PickOff[5]= np.array(im05)
PickOff[6]= np.array(im06)
PickOff[7]= np.array(im07)
PickOff[8]= np.array(im08)
PickOff[9]= np.array(im09)
PickOff[10]= np.array(im09)#
PickOff[11]= np.array(im09)#
PickOff[12]= np.array(im09)#
PickOff[13]= np.array(im13)
PickOff[14]= np.array(im14)
PickOff[15]= np.array(im15)
PickOff[16]= np.array(im16)
PickOff[17]= np.array(im17)
PickOff[18]= np.array(im18)
PickOff[19]= np.array(im19)
PickOff[20]= np.array(im20)
PickOff[21]= np.array(im21)
PickOff[22]= np.array(im22)
PickOff[23]= np.array(im23)
PickOff[24]= np.array(im24)
PickOff[25]= np.array(im25)
PickOff[26]= np.array(im26)
PickOff[27]= np.array(im27)
PickOff[28]= np.array(im09)#
PickOff[29]= np.array(im29)
PickOff[30]= np.array(im30)
PickOff[31]= np.array(im31)
PickOff[32]= np.array(im32)
PickOff[33]= np.array(im33)
PickOff[34]= np.array(im09)#
PickOff[35]= np.array(im09)#
PickOff[36]= np.array(im36)
PickOff[37]= np.array(im09)#
PickOff[38]= np.array(im38)
PickOff[39]= np.array(im39)
PickOff[40]= np.array(im40)
PickOff[41]= np.array(im41)

#%%
NorInt={}
for n in range (1, 42):
    NorInt[n]=PickOff[n]/sum(sum(PickOff[n]))
#%%
pixel=3.45e-6
xaxis=np.linspace(-1023, 1024, 2048)*pixel*10**3
yaxis=np.linspace(-767, 768, 1536)*pixel*10**3
#%%
plt.figure(2)
plt.pcolormesh(xaxis, yaxis, NorInt[41])
plt.axis('scaled')
plt.title('Low Power at 67cm Down Stream')
plt.colorbar(label='Normalized Int')
plt.xlabel('mm')
plt.ylabel('mm')
plt.savefig('LowP67cm')

#%%

