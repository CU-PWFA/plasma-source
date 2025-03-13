#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 15:54:42 2022

@author: chris
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 13:13:13 2022

Looking at the 0 psi dataset, find ybeam vs z.  Also calculate the error in gamma

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
import sys
#from scipy import optimize
from scipy.optimize import curve_fit
import scipy.io as sio
import matplotlib.colors as colors

"""
def FitDataSomething(data, axis, function, guess = [0.,0.,0.]):
    errfunc = lambda p, x, y: function(p, x) - y
    p0 = guess
    p1, success = optimize.leastsq(errfunc,p0[:], args=(axis, data))
    return p1
"""

def Gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

#superpath = '/home/chris/Desktop/SLACData/'
superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'
dataset = 'E308_02493' #no gas
#dataset = 'E308_02512' #1 psi
#dataset = 'E308_02492' #6psi, lowres
#dataset = 'E308_02506' #9 bar over 0, 115.8 psi
#dataset = 'E308_02497' #6 psi, highres

dx = 30.3e-6*1e3
gasjetpos = 1993.3

path = superpath + day + dataset + '/' + dataset +'.h5'
f = h5.File(path,"r") 
data = f['data']
save_info = data['save_info']
params = data['params']
func_start = params['startVals'][0][0]
func_end = params['stopVals'][0][0]
camerapath = 'DTOTR2'
path = superpath + day + dataset + '/images/' + camerapath + '/'
image_list = os.listdir(path)

n_shot = int(params['n_shot'][0][0])
n_step = int(len(image_list)/n_shot)
f.close()

ybeam_arr = np.zeros((n_step,n_shot))
qbeam_arr = np.zeros((n_step,n_shot))
#step = 1
for s in range(n_step):
    step = s+1
    for n in range(n_shot):
        nimg = n
        
        path = superpath + day + dataset + '/' + dataset +'.h5'
        ###BACKGROUND FROM MAT FILE
        
        path = superpath + day + dataset + '/' + dataset +'.mat'
        mat = sio.loadmat(path)
        
        data = mat['data_struct']
        backs = data['backgrounds'][0][0]
        dtotr2_back = backs['DTOTR2'][0][0]
        dtotr2_back = np.transpose(dtotr2_back)
        #######################
        camerapath = 'DTOTR2'
        if step < 10:
            image_path = camerapath + "_data_step0" + str(step) + "_000" + str(nimg) + ".tif"
        else:
            image_path = camerapath + "_data_step" + str(step) + "_000" + str(nimg) + ".tif"
        
        path = superpath + day + dataset + '/images/' + camerapath + '/' + image_path
        
        image = plt.imread(path)
        
        im_arr = np.array(image)
        
        threshold = 20
        for i in range(im_arr.shape[0]):
            for j in range(im_arr.shape[1]):
                im_arr[i][j] = im_arr[i][j] - dtotr2_back[i][j]
                if im_arr[i][j] > 65000:
                    im_arr[i][j] = 0
                if im_arr[i][j] < threshold:
                    im_arr[i][j] = 0
        
        hori_range = np.array([0,im_arr.shape[1]])*dx
        en_arr = np.linspace(hori_range[0],hori_range[1],im_arr.shape[1])
        proj_arr = np.zeros(im_arr.shape[1])
        for i in range(im_arr.shape[1]):
            proj_arr[i] = np.sum(im_arr[:,i])
        int_arr = np.zeros(len(proj_arr))
        for i in range(len(proj_arr)):
            if i == 0:
                int_arr[0] = proj_arr[0]
            else:
                int_arr[i] = int_arr[i-1]+proj_arr[i]
        qbeam_arr[s,n]=int_arr[-1]
        int_arr = int_arr/int_arr[-1]
        cent_ind = (np.where(int_arr<0.5)[0])[-1]
        ybeam_arr[s,n]=en_arr[cent_ind]

ybeam_ave = np.zeros(n_step)
ybeam_std = np.zeros(n_step)
z_arr = np.linspace(func_start, func_end, n_step)-gasjetpos
for i in range(len(ybeam_ave)):
    ybeam_ave[i] = np.average(ybeam_arr[i,:])
    ybeam_std[i] = np.sqrt(np.var(ybeam_arr[i,:]))
    
qbeam_ave = np.zeros(n_step)
qbeam_std = np.zeros(n_step)
for i in range(len(qbeam_ave)):
    qbeam_ave[i] = np.average(qbeam_arr[i,:])
    qbeam_std[i] = np.sqrt(np.var(qbeam_arr[i,:]))

#lfit = np.polyfit(z_arr,ybeam_ave,1)

#plt.plot(z_arr,ybeam_ave)
plt.errorbar(z_arr,ybeam_ave,ybeam_std,label=dataset)
plt.xlabel("z (m)")
plt.ylabel("Energy Centroid on Screen (mm)")
plt.legend(); plt.grid(); plt.show()

plt.errorbar(z_arr,qbeam_ave,qbeam_std,label=dataset)
plt.xlabel("z (m)")
plt.ylabel("Total Signal on Screen (Charge?)")
plt.legend(); plt.grid(); plt.show()