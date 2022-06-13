#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:33:16 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
#P=4.9e17 and Off
im029 = Image.open('029/19423601_2001230029.tiff')
im030 = Image.open('030/19423601_2001230030.tiff')
im032 = Image.open('032/19423601_2001230032.tiff')
im033 = Image.open('033/19423601_2001230033.tiff')
im034 = Image.open('034/19423601_2001230034.tiff')
im035 = Image.open('035/19423601_2001230035.tiff')
im036 = Image.open('036/19423601_2001230036.tiff')
im037 = Image.open('037/19423601_2001230037.tiff')
im038 = Image.open('038/19423601_2001230038.tiff')
im039 = Image.open('039/19423601_2001230039.tiff')
im040 = Image.open('040/19423601_2001230040.tiff')
im041 = Image.open('041/19423601_2001230041.tiff')
im042 = Image.open('042/19423601_2001230042.tiff')
im043 = Image.open('043/19423601_2001230043.tiff')
im044 = Image.open('044/19423601_2001230044.tiff')
im045 = Image.open('045/19423601_2001230045.tiff')
im046 = Image.open('046/19423601_2001230046.tiff')
im047 = Image.open('047/19423601_2001230047.tiff')
#P=1e18
im049 = Image.open('049/19423601_2001230049.tiff')
im053 = Image.open('053/19423601_2001230053.tiff')
im054 = Image.open('054/19423601_2001230054.tiff')
im057 = Image.open('057/19423601_2001230057.tiff')
im058 = Image.open('058/19423601_2001230058.tiff')
im059 = Image.open('059/19423601_2001230059.tiff')
im060 = Image.open('060/19423601_2001230060.tiff')
im061 = Image.open('061/19423601_2001230061.tiff')
im062 = Image.open('062/19423601_2001230062.tiff')
im063 = Image.open('063/19423601_2001230063.tiff')
#P=1.9e17
im014 = Image.open('014/19423601_2001230014.tiff')
im017 = Image.open('017/19423601_2001230017.tiff')
im018 = Image.open('018/19423601_2001230018.tiff')
im021 = Image.open('021/19423601_2001230021.tiff')
im023 = Image.open('023/19423601_2001230023.tiff')
im026 = Image.open('026/19423601_2001230026.tiff')
im027 = Image.open('027/19423601_2001230027.tiff')
#P=8e16
im064 = Image.open('064/19423601_2001230064.tiff')
im067 = Image.open('067/19423601_2001230067.tiff')
im068 = Image.open('068/19423601_2001230068.tiff')
im069 = Image.open('069/19423601_2001230069.tiff')
im070 = Image.open('070/19423601_2001230070.tiff')
im071 = Image.open('071/19423601_2001230071.tiff')
im072 = Image.open('072/19423601_2001230072.tiff')
im073 = Image.open('073/19423601_2001230073.tiff')
#P=2.5e16
im074 = Image.open('074/19423601_2001230074.tiff')
im075 = Image.open('075/19423601_2001230075.tiff')
im076 = Image.open('076/19423601_2001230076.tiff')
im077 = Image.open('077/19423601_2001230077.tiff')
im078 = Image.open('078/19423601_2001230078.tiff')
im079 = Image.open('079/19423601_2001230079.tiff')
#P=5.1e16
im080 = Image.open('080/19423601_2001230080.tiff')
im081 = Image.open('081/19423601_2001230081.tiff')
im082 = Image.open('082/19423601_2001230082.tiff')
im083 = Image.open('083/19423601_2001230083.tiff')
im084 = Image.open('084/19423601_2001230084.tiff')
im085 = Image.open('085/19423601_2001230085.tiff')

#%%
im089 = Image.open('089/19423601_2001230089.tiff')
im090 = Image.open('090/19423601_2001230090.tiff')
im091 = Image.open('091/19423601_2001230091.tiff')
im092 = Image.open('092/19423601_2001230092.tiff')
im093 = Image.open('093/19423601_2001230093.tiff')
im094 = Image.open('094/19423601_2001230094.tiff')
im095 = Image.open('095/19423601_2001230095.tiff')
im096 = Image.open('096/19423601_2001230096.tiff')
im097 = Image.open('097/19423601_2001230097.tiff')
im098 = Image.open('098/19423601_2001230098.tiff')
im099 = Image.open('099/19423601_2001230099.tiff')
im100 = Image.open('100/19423601_2001230100.tiff')

#P=3.5e17
On_4e17_500_71= np.array(im089) #51
On_4e17_500_70= np.array(im090) #51
On_4e17_500_69= np.array(im091) #51
On_4e17_500_68= np.array(im092) #51
On_4e17_500_67= np.array(im093) #51
On_4e17_500_66= np.array(im094) #51
On_4e17_500_65= np.array(im095) #51
On_4e17_500_64= np.array(im096) #51
On_4e17_500_63= np.array(im097) #51
On_4e17_500_62= np.array(im100) #51

#%%
im015 = Image.open('015/19423601_2001230015.tiff')
test015= np.array(im015) #51

#%%
#Off
Off_500= np.array(im030) #51
Off_485= np.array(im076) #57
Off_470= np.array(im032) #58
Off_460= np.array(im079) #59
Off_455= np.array(im085) #60
Off_435= np.array(im036) #61
Off_400= np.array(im036) #62
Off_365= np.array(im039) #63
Off_340= np.array(im040) #64
Off_305= np.array(im043) #65
Off_260= np.array(im044) #66
Off_200= np.array(im047) #67
Off_180= np.array(im063) #68

#P=1e18
On_1e18_500= np.array(im049) #51
On_1e18_470= np.array(im053) #58
On_1e18_435= np.array(im054) #61
On_1e18_365= np.array(im057) #63
On_1e18_340= np.array(im058) #64
On_1e18_305= np.array(im059) #65
On_1e18_260= np.array(im060) #66
On_1e18_200= np.array(im061) #67
On_1e18_180= np.array(im062) #68

#P=4.9e17
On_5e17_500= np.array(im029) #51
On_5e17_470= np.array(im033) #58
On_5e17_435= np.array(im034) #61
On_5e17_400= np.array(im037) #62
On_5e17_365= np.array(im038) #63
On_5e17_340= np.array(im041) #64
On_5e17_305= np.array(im042) #65
On_5e17_260= np.array(im045) #66
On_5e17_200= np.array(im046) #67

#P=1.9e17
On_2e17_500= np.array(im014) #51
On_2e17_470= np.array(im017) #58
On_2e17_435= np.array(im018) #61
On_2e17_365= np.array(im021) #63
On_2e17_340= np.array(im023) #64
On_2e17_305= np.array(im026) #65
On_2e17_260= np.array(im027) #66

#P=8e16
On_8e16_500= np.array(im064) #51
On_8e16_470= np.array(im067) #58
On_8e16_435= np.array(im068) #61
On_8e16_400= np.array(im069) #62
On_8e16_365= np.array(im070) #63
On_8e16_340= np.array(im071) #64
On_8e16_305= np.array(im072) #65
On_8e16_260= np.array(im073) #66

#P=5.1e16
On_5e16_500= np.array(im080) #51
On_5e16_485= np.array(im081) #57
On_5e16_470= np.array(im082) #58
On_5e16_460= np.array(im083) #59
On_5e16_455= np.array(im084) #60

#P=2.5e16
On_3e16_500= np.array(im074) #51
On_3e16_485= np.array(im075) #57
On_3e16_470= np.array(im077) #58
On_3e16_460= np.array(im078) #59

#%%
#500->1
#485->0.97
#470->0.94
#460->0.92
#455->0.91
#435->0.87
#400->0.8
#365->0.73
#340->0.68
#305->0.61
#260->0.52
#200->0.4
#180->0.36
#%%
pixel=3.45e-6
Energy=1
NorE=On_4e17_500_62*Energy/sum(sum(On_4e17_500_62))*1e4
#NorE=test015*Energy/sum(sum(test015))*1e4

#NorInt=(NorE/30e-15)/((pixel)**2)

xaxis=np.linspace(-1023, 1024, 2048)*pixel*10**3/0.24
yaxis=np.linspace(-767, 768, 1536)*pixel*10**3/0.24

#NorInt=NorInt/10000/1e12
#%%
fig=plt.figure(1)
plt.pcolormesh(xaxis, yaxis, NorE)
plt.axis('scaled')
plt.colorbar(label='Normalized energy')
plt.xlabel('mm')
plt.ylabel('mm')
plt.title('$D=3.5e16cm^{-3}$ E=500mJ d_I=62cm')
#plt.title('PlasmaOff E=500mJ')
#plt.title('Signal-background_455')
#%%
fig.savefig('62cm.png')

#%%

Diff=abs(On_5e16_500-Off_500)

#%%
plt.figure(2)
plt.pcolormesh(On_5e16_500)
