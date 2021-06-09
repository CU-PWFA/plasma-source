#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  8 17:47:15 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import find_peaks
from sklearn.linear_model import LinearRegression

#%%
Den= [2, 5, 8, 10, 20, 40, 50, 80]
Pow= [5]
PWidth= [30, 35, 40, 45, 50]

#%%
PixelSize= 10

count= 0

XMatrix= []
BasePath_AG= '/media/valentinalee/TRANSCEND/AfterglowFinalPowerFix/ScanResults/AfterglowScan/'
BasePath_Sh= '/media/valentinalee/TRANSCEND/ShadowgraphyResults/FFResults2/'
#%%
for density in Den:
    for power in Pow:
        for pw in PWidth:
            try:
                
                PhotonMap= np.load(BasePath_AG+'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_map.npy')
                IntAG= PhotonMap[:, int(PhotonMap.shape[1]-1)]
                r= np.load(BasePath_AG+'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'_x.npy')
                xvar= r[int(len(r)/2):len(r)]
                yvar= IntAG[int(len(IntAG)/2):len(IntAG)]
                fitfn= lambda p, xvar: p[0]*np.exp(-(xvar**2/(2*p[1]**2))**p[2])
    
                errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar;
                p0= [1e11, 20, 2];
                p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))

                dataE= np.load(BasePath_Sh+'PW'+str(pw)+'_GP'+str(power)+'_n'+str(density)+'.npy')
                MidLineOutI= abs(dataE[:, int(dataE.shape[1]/2)])**2

                x= MidLineOutI
                peaks, _ = find_peaks(x)
                peak1idx, peak2idx= peaks[np.argsort(x[peaks])[-1]], peaks[np.argsort(x[peaks])[-2]]
                peak1, peak2= x[peak1idx], x[peak2idx]
                
                if (np.argsort(x[peaks])[-1])-(np.argsort(x[peaks])[-2])==2:
                   
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx, valley2idx= valleys[np.argsort(x[valleys])[0:2]]
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                else:
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx= valleys[np.argsort(x[valleys])[0]]
                    valley2idx= valley1idx
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                    
                if (np.mean(np.array([peak1, peak2])))-(np.mean(np.array([valley1, valley2])))>0.5:
                    
                    plt.plot(x)
                    plt.plot(np.array([peak1idx, peak2idx]), x[np.array([peak1idx, peak2idx])], 'x')
                    plt.show()     
                    plt.plot(np.array([valley1idx, valley2idx]), x[np.array([valley1idx, valley2idx])], 'x')
    
                    FirstPeakWidth= abs(peak1idx-peak2idx)*PixelSize
                    PTValue= 100*(np.mean(np.array([peak1, peak2]))-np.mean(np.array([valley1, valley2])))/ \
                                ((np.mean(np.array([peak1, peak2]))+np.mean(np.array([valley1, valley2])))/2)
                    XMatrix.append([density, pw, p1[1], p1[0], FirstPeakWidth, PTValue])
                    
                    count= count+1

            except:
                print(density, power, pw)
                print('Oh no!')
                pass
                
XMatrix= np.array(XMatrix)            
#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
NewXMatrix[:, 0]= (XMatrix[:, 2]-np.amin(XMatrix[:, 2]))/(np.amax(XMatrix[:, 2])-np.amin(XMatrix[:, 2]))
NewXMatrix[:, 1]= (XMatrix[:, 3]-np.amin(XMatrix[:, 3]))/(np.amax(XMatrix[:, 3])-np.amin(XMatrix[:, 3]))
NewXMatrix[:, 2]= (XMatrix[:, 4]-np.amin(XMatrix[:, 4]))/(np.amax(XMatrix[:, 4])-np.amin(XMatrix[:, 4]))
NewXMatrix[:, 3]= (XMatrix[:, 5]-np.amin(XMatrix[:, 5]))/(np.amax(XMatrix[:, 5])-np.amin(XMatrix[:, 5]))
NewXMatrix[:, 4]= (XMatrix[:, 2]*XMatrix[:, 3]-np.amin(XMatrix[:, 2]*XMatrix[:, 3]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 3])-np.amin(XMatrix[:, 2]*XMatrix[:, 3]))
NewXMatrix[:, 5]= (XMatrix[:, 2]*XMatrix[:, 4]-np.amin(XMatrix[:, 2]*XMatrix[:, 4]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 4])-np.amin(XMatrix[:, 2]*XMatrix[:, 4]))
NewXMatrix[:, 6]= (XMatrix[:, 2]*XMatrix[:, 5]-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 5])-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))
NewXMatrix[:, 7]= (XMatrix[:, 3]*XMatrix[:, 4]-np.amin(XMatrix[:, 3]*XMatrix[:, 4]))/ \
                    (np.amax(XMatrix[:, 3]*XMatrix[:, 4])-np.amin(XMatrix[:, 3]*XMatrix[:, 4]))
NewXMatrix[:, 8]= (XMatrix[:, 3]*XMatrix[:, 5]-np.amin(XMatrix[:, 3]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 3]*XMatrix[:, 5])-np.amin(XMatrix[:, 3]*XMatrix[:, 5]))
NewXMatrix[:, 9]= (XMatrix[:, 4]*XMatrix[:, 5]-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 4]*XMatrix[:, 5])-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))
NewXMatrix[:, 10]= (XMatrix[:, 2]**2-np.amin(XMatrix[:, 2]**2))/(np.amax(XMatrix[:, 2]**2)-np.amin(XMatrix[:, 2]**2))
NewXMatrix[:, 11]= (XMatrix[:, 3]**2-np.amin(XMatrix[:, 3]**2))/(np.amax(XMatrix[:, 3]**2)-np.amin(XMatrix[:, 3]**2))
NewXMatrix[:, 12]= (XMatrix[:, 4]**2-np.amin(XMatrix[:, 4]**2))/(np.amax(XMatrix[:, 4]**2)-np.amin(XMatrix[:, 4]**2))
NewXMatrix[:, 13]= (XMatrix[:, 5]**2-np.amin(XMatrix[:, 5]**2))/(np.amax(XMatrix[:, 5]**2)-np.amin(XMatrix[:, 5]**2))


#%%
#Do a multiple regression
n_reg= LinearRegression().fit(NewXMatrix, XMatrix[:, 0])
n_coef= n_reg.coef_
n_intercept= n_reg.intercept_
n_score= n_reg.score(NewXMatrix, XMatrix[:, 0])

w_reg= LinearRegression().fit(NewXMatrix, XMatrix[:, 1])
w_coef= w_reg.coef_
w_intercept= w_reg.intercept_
w_score= w_reg.score(NewXMatrix, XMatrix[:, 1])

#%%
for pop_idx in range(0, 14):
    NewXM_pop= np.delete(NewXMatrix, pop_idx, 1)
#    n_true_pop= np.delete(XMatrix[:, 0], pop_idx, 0)
#    w_true_pop= np.delete(XMatrix[:, 1], pop_idx, 0)

    n_reg_pop= LinearRegression().fit(NewXM_pop, XMatrix[:, 0])
    n_score_pop= n_reg_pop.score(NewXM_pop, XMatrix[:, 0])
    w_reg_pop= LinearRegression().fit(NewXM_pop, XMatrix[:, 1])
    w_score_pop= w_reg_pop.score(NewXM_pop, XMatrix[:, 1])
    
    print(n_score_pop, w_score_pop)
    if n_score_pop>= n_score:
        print(pop_idx, 'n')
    if w_score_pop>= w_score:
        print(pop_idx, 'w')
#%%
n_predict= n_coef[0]*NewXMatrix[:, 0]+ n_coef[1]*NewXMatrix[:, 1]+ n_coef[2]*NewXMatrix[:, 2]+ \
           n_coef[3]*NewXMatrix[:, 3]+ n_coef[4]*NewXMatrix[:, 4]+ n_coef[5]*NewXMatrix[:, 5]+ \
           n_coef[6]*NewXMatrix[:, 6]+ n_coef[7]*NewXMatrix[:, 7]+ n_coef[8]*NewXMatrix[:, 8]+ \
           n_coef[9]*NewXMatrix[:, 9]+ n_coef[10]*NewXMatrix[:, 10]+ n_coef[11]*NewXMatrix[:, 11]+\
           n_coef[12]*NewXMatrix[:, 12]+ n_coef[13]*NewXMatrix[:, 13]+ n_intercept

w_predict= w_coef[0]*NewXMatrix[:, 0]+ w_coef[1]*NewXMatrix[:, 1]+ w_coef[2]*NewXMatrix[:, 2]+ \
           w_coef[3]*NewXMatrix[:, 3]+ w_coef[4]*NewXMatrix[:, 4]+ w_coef[5]*NewXMatrix[:, 5]+ \
           w_coef[6]*NewXMatrix[:, 6]+ w_coef[7]*NewXMatrix[:, 7]+ w_coef[8]*NewXMatrix[:, 8]+ \
           w_coef[9]*NewXMatrix[:, 9]+ w_coef[10]*NewXMatrix[:, 10]+ w_coef[11]*NewXMatrix[:, 11]+ \
           w_coef[12]*NewXMatrix[:, 12]+ w_coef[13]*NewXMatrix[:, 13]+ w_intercept

#%%
errorM= np.zeros((n_predict.shape[0], 4))
errorM[:, 0]= XMatrix[:, 0]
errorM[:, 1]= XMatrix[:, 1]
errorM[:, 2]= abs(XMatrix[:, 0]-n_predict)/XMatrix[:, 0]
errorM[:, 3]= abs(XMatrix[:, 1]-w_predict)/XMatrix[:, 1]
#%%
errorM= np.round(errorM, 4)
#%%
plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 2], s= 500)
plt.colorbar()
plt.xscale('log')
plt.title('Density Error (14 terms)')
plt.xlabel('True Density (1e16 cm$^{-3}$)')
plt.ylabel('True Plasma Width ($\mu$m)')
#%%
plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 3], s= 500)
plt.colorbar()
plt.xscale('log')
plt.title('Width Error (14 terms)')
plt.xlabel('True Density (1e16 cm$^{-3}$)')
plt.ylabel('True Plasma Width ($\mu$m)')

#%%
#%%
n_s= np.zeros(14)
for k in range(14):
    n_reg= LinearRegression().fit(np.reshape(NewXMatrix[:, k], (len(NewXMatrix[:, k]), 1)), XMatrix[:, 0])
    n_score= n_reg.score(np.reshape(NewXMatrix[:, k], (len(NewXMatrix[:, k]), 1)), XMatrix[:, 0])
    n_s[k]= n_score
    print(n_score)
#%%Density R^2
n_importance_rank= np.argsort(abs(n_coef))[::-1]
n_importance_rank= np.argsort(abs(n_s))[::-1]

n_score= np.zeros(len(n_coef))
current_X_n= np.reshape(NewXMatrix[:, n_importance_rank[0]], (len(NewXMatrix[:, n_importance_rank[0]]), 1))
#%%
for idx in range (int(len(n_coef)-1)):
    
    n_reg= LinearRegression().fit(current_X_n, XMatrix[:, 0])
    n_score[idx]= n_reg.score(current_X_n, XMatrix[:, 0])

    current_X_n= np.append(current_X_n, np.reshape(NewXMatrix[:, n_importance_rank[int(idx+1)]], \
                         (len(NewXMatrix[:, n_importance_rank[int(idx+1)]]), 1)), axis= 1)
n_reg= LinearRegression().fit(current_X_n, XMatrix[:, 0])
n_score[int(idx+1)]= n_reg.score(current_X_n, XMatrix[:, 0])
n_intercept= n_reg.intercept_
#%%
plt.plot(n_score[1:], '.')
#%%
xvar= np.arange(len(n_score[1:]))+2
yvar= n_score[1:]
fitfn= lambda p, xvar: p[0]/(1+np.exp(-p[1]*(xvar-p[2])))+p[3]

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
p0= [0.3, 1, 4, -1]
p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
#%%
plot_x= np.linspace(np.amin(xvar), np.amax(xvar), 1000)
plt.plot(xvar, yvar, 'o')
plt.plot(plot_x, fitfn(p1, plot_x))
plt.xlabel('Terms')
plt.ylabel('R Square')
plt.title('Density Fitting')
#%%
Sat_value= p1[0]+p1[3]
Sat_term= np.log(p1[0]/(Sat_value-p1[3])-1)/-p1[1]+p1[2]

#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  7))

#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
#NewXMatrix[:, 0]= (XMatrix[:, 3]-np.amin(XMatrix[:, 3]))/(np.amax(XMatrix[:, 3])-np.amin(XMatrix[:, 3])) #7
#NewXMatrix[:, 1]= (XMatrix[:, 4]-np.amin(XMatrix[:, 4]))/(np.amax(XMatrix[:, 4])-np.amin(XMatrix[:, 4])) #3
#NewXMatrix[:, 2]= (XMatrix[:, 5]-np.amin(XMatrix[:, 5]))/(np.amax(XMatrix[:, 5])-np.amin(XMatrix[:, 5])) #4
#NewXMatrix[:, 3]= (XMatrix[:, 2]*XMatrix[:, 5]-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))/ \
#                    (np.amax(XMatrix[:, 2]*XMatrix[:, 5])-np.amin(XMatrix[:, 2]*XMatrix[:, 5])) #6
#NewXMatrix[:, 4]= (XMatrix[:, 3]*XMatrix[:, 4]-np.amin(XMatrix[:, 3]*XMatrix[:, 4]))/ \
 #                   (np.amax(XMatrix[:, 3]*XMatrix[:, 4])-np.amin(XMatrix[:, 3]*XMatrix[:, 4])) #8
#NewXMatrix[:, 7]= (XMatrix[:, 4]*XMatrix[:, 5]-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))/ \
#                    (np.amax(XMatrix[:, 4]*XMatrix[:, 5])-np.amin(XMatrix[:, 4]*XMatrix[:, 5])) #2
NewXMatrix[:, 8]= (XMatrix[:, 4]**2-np.amin(XMatrix[:, 4]**2))/(np.amax(XMatrix[:, 4]**2)-np.amin(XMatrix[:, 4]**2)) #1
#NewXMatrix[:, 3]= (XMatrix[:, 5]**2-np.amin(XMatrix[:, 5]**2))/(np.amax(XMatrix[:, 5]**2)-np.amin(XMatrix[:, 5]**2)) #5

#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
NewXMatrix[:, 0]= (XMatrix[:, 3]-np.amin(XMatrix[:, 3]))/(np.amax(XMatrix[:, 3])-np.amin(XMatrix[:, 3])) #7
NewXMatrix[:, 1]= (XMatrix[:, 4]-np.amin(XMatrix[:, 4]))/(np.amax(XMatrix[:, 4])-np.amin(XMatrix[:, 4])) #3n
NewXMatrix[:, 2]= (XMatrix[:, 5]-np.amin(XMatrix[:, 5]))/(np.amax(XMatrix[:, 5])-np.amin(XMatrix[:, 5])) #4n
NewXMatrix[:, 3]= (XMatrix[:, 2]*XMatrix[:, 5]-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 5])-np.amin(XMatrix[:, 2]*XMatrix[:, 5])) #6
NewXMatrix[:, 4]= (XMatrix[:, 4]*XMatrix[:, 5]-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 4]*XMatrix[:, 5])-np.amin(XMatrix[:, 4]*XMatrix[:, 5])) #2n
NewXMatrix[:, 5]= (XMatrix[:, 4]**2-np.amin(XMatrix[:, 4]**2))/(np.amax(XMatrix[:, 4]**2)-np.amin(XMatrix[:, 4]**2)) #1n
NewXMatrix[:, 6]= (XMatrix[:, 5]**2-np.amin(XMatrix[:, 5]**2))/(np.amax(XMatrix[:, 5]**2)-np.amin(XMatrix[:, 5]**2)) #5n

#%%
#Do a multiple regression
n_reg= LinearRegression().fit(NewXMatrix, XMatrix[:, 0])
n_coef= n_reg.coef_
n_intercept= n_reg.intercept_
n_score= n_reg.score(NewXMatrix, XMatrix[:, 0])

#%%
w_s= np.zeros(14)
for k in range(14):
    w_reg= LinearRegression().fit(np.reshape(NewXMatrix[:, k], (len(NewXMatrix[:, k]), 1)), XMatrix[:, 1])
    w_coef= w_reg.coef_
    w_intercept= w_reg.intercept_
    w_score= w_reg.score(np.reshape(NewXMatrix[:, k], (len(NewXMatrix[:, k]), 1)), XMatrix[:, 1])
    w_s[k]= w_score
    print(w_score)
#%%
#%%Width R^2
w_importance_rank= np.argsort(abs(w_coef))[::-1]
w_importance_rank= np.argsort(abs(w_s))[::-1]

w_score= np.zeros(len(w_coef))
current_X_w= np.reshape(NewXMatrix[:, w_importance_rank[0]], (len(NewXMatrix[:, w_importance_rank[0]]), 1))
#%%
for idx in range (int(len(w_coef)-1)):
    print(idx)
    w_reg= LinearRegression().fit(current_X_w, XMatrix[:, 1])
    w_score[idx]= w_reg.score(current_X_w, XMatrix[:, 1])

    current_X_w= np.append(current_X_w, np.reshape(NewXMatrix[:, w_importance_rank[int(idx+1)]], \
                         (len(NewXMatrix[:, w_importance_rank[int(idx+1)]]), 1)), axis= 1)
w_reg= LinearRegression().fit(current_X_w, XMatrix[:, 1])
w_score[int(idx+1)]= w_reg.score(current_X_w, XMatrix[:, 1])
#%%
plt.plot(w_score[1:], '.')
#%%
xvar= np.arange(len(w_score[1:]))+2
yvar= w_score[1:]
fitfn= lambda p, xvar: p[0]/(1+np.exp(-p[1]*(xvar-p[2])))+p[3]

errfunc = lambda p, xvar, yvar: fitfn(p, xvar) - yvar
p0= [1, 1, 4, 0]
p1, success = optimize.leastsq(errfunc, p0[:], args=(xvar, yvar))
#%%
plot_x= np.linspace(np.amin(xvar), np.amax(xvar), 1000)
plt.plot(xvar, yvar, 'o')
plt.plot(plot_x, fitfn(p1, plot_x))
plt.xlabel('Terms')
plt.ylabel('R Square')
plt.title('Width Fitting')

#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
NewXMatrix[:, 0]= (XMatrix[:, 2]-np.amin(XMatrix[:, 2]))/(np.amax(XMatrix[:, 2])-np.amin(XMatrix[:, 2])) #1
NewXMatrix[:, 1]= (XMatrix[:, 3]-np.amin(XMatrix[:, 3]))/(np.amax(XMatrix[:, 3])-np.amin(XMatrix[:, 3]))
NewXMatrix[:, 2]= (XMatrix[:, 4]-np.amin(XMatrix[:, 4]))/(np.amax(XMatrix[:, 4])-np.amin(XMatrix[:, 4]))
NewXMatrix[:, 3]= (XMatrix[:, 5]-np.amin(XMatrix[:, 5]))/(np.amax(XMatrix[:, 5])-np.amin(XMatrix[:, 5]))
NewXMatrix[:, 4]= (XMatrix[:, 2]*XMatrix[:, 3]-np.amin(XMatrix[:, 2]*XMatrix[:, 3]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 3])-np.amin(XMatrix[:, 2]*XMatrix[:, 3])) #6
NewXMatrix[:, 5]= (XMatrix[:, 2]*XMatrix[:, 4]-np.amin(XMatrix[:, 2]*XMatrix[:, 4]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 4])-np.amin(XMatrix[:, 2]*XMatrix[:, 4])) #3
NewXMatrix[:, 6]= (XMatrix[:, 2]*XMatrix[:, 5]-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 5])-np.amin(XMatrix[:, 2]*XMatrix[:, 5])) #4
NewXMatrix[:, 7]= (XMatrix[:, 3]*XMatrix[:, 4]-np.amin(XMatrix[:, 3]*XMatrix[:, 4]))/ \
                    (np.amax(XMatrix[:, 3]*XMatrix[:, 4])-np.amin(XMatrix[:, 3]*XMatrix[:, 4]))
NewXMatrix[:, 8]= (XMatrix[:, 3]*XMatrix[:, 5]-np.amin(XMatrix[:, 3]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 3]*XMatrix[:, 5])-np.amin(XMatrix[:, 3]*XMatrix[:, 5])) #7
NewXMatrix[:, 9]= (XMatrix[:, 4]*XMatrix[:, 5]-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 4]*XMatrix[:, 5])-np.amin(XMatrix[:, 4]*XMatrix[:, 5]))
NewXMatrix[:, 10]= (XMatrix[:, 2]**2-np.amin(XMatrix[:, 2]**2))/(np.amax(XMatrix[:, 2]**2)-np.amin(XMatrix[:, 2]**2)) #2
NewXMatrix[:, 11]= (XMatrix[:, 3]**2-np.amin(XMatrix[:, 3]**2))/(np.amax(XMatrix[:, 3]**2)-np.amin(XMatrix[:, 3]**2)) #5
NewXMatrix[:, 12]= (XMatrix[:, 4]**2-np.amin(XMatrix[:, 4]**2))/(np.amax(XMatrix[:, 4]**2)-np.amin(XMatrix[:, 4]**2))
NewXMatrix[:, 13]= (XMatrix[:, 5]**2-np.amin(XMatrix[:, 5]**2))/(np.amax(XMatrix[:, 5]**2)-np.amin(XMatrix[:, 5]**2))

#%%
NewXMatrix= np.zeros((XMatrix.shape[0],  14))
#NewXMatrix[:, 0]= (XMatrix[:, 2]-np.amin(XMatrix[:, 2]))/(np.amax(XMatrix[:, 2])-np.amin(XMatrix[:, 2])) #1
NewXMatrix[:, 0]= (XMatrix[:, 2]*XMatrix[:, 4]-np.amin(XMatrix[:, 2]*XMatrix[:, 4]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 4])-np.amin(XMatrix[:, 2]*XMatrix[:, 4])) #3
NewXMatrix[:, 3]= (XMatrix[:, 2]*XMatrix[:, 5]-np.amin(XMatrix[:, 2]*XMatrix[:, 5]))/ \
                    (np.amax(XMatrix[:, 2]*XMatrix[:, 5])-np.amin(XMatrix[:, 2]*XMatrix[:, 5])) #4
NewXMatrix[:, 1]= (XMatrix[:, 2]**2-np.amin(XMatrix[:, 2]**2))/(np.amax(XMatrix[:, 2]**2)-np.amin(XMatrix[:, 2]**2)) #2
NewXMatrix[:, 2]= (XMatrix[:, 3]**2-np.amin(XMatrix[:, 3]**2))/(np.amax(XMatrix[:, 3]**2)-np.amin(XMatrix[:, 3]**2)) #5

#%%
w_reg= LinearRegression().fit(NewXMatrix, XMatrix[:, 1])
w_coef= w_reg.coef_
w_intercept= w_reg.intercept_
w_score= w_reg.score(NewXMatrix, XMatrix[:, 1])

#%%
#%%
n_predict= n_coef[8]*NewXMatrix[:, 8]+ n_intercept
#%%
n_predict= n_coef[0]*NewXMatrix[:, 0]+ n_coef[1]*NewXMatrix[:, 1]+ n_coef[2]*NewXMatrix[:, 2]+ \
           n_coef[3]*NewXMatrix[:, 3]+ n_coef[4]*NewXMatrix[:, 4]+ n_coef[5]*NewXMatrix[:, 5]+ \
           n_coef[6]*NewXMatrix[:, 6]+  n_intercept
#%%
w_predict= w_coef[0]*NewXMatrix[:, 0]+ w_coef[1]*NewXMatrix[:, 1]+ w_coef[2]*NewXMatrix[:, 2]+ \
           w_coef[3]*NewXMatrix[:, 3]+ w_intercept
#%%
n_predict= n_coef[1]*NewXMatrix[:, 1]+ n_coef[3]*NewXMatrix[:, 3]+ n_coef[4]*NewXMatrix[:, 4]+  \
           n_coef[6]*NewXMatrix[:, 6]+ n_coef[8]*NewXMatrix[:, 8]+ n_coef[10]*NewXMatrix[:, 10]+ \
           n_coef[11]*NewXMatrix[:, 11]+ n_intercept

w_predict= w_coef[0]*NewXMatrix[:, 0]+ w_coef[2]*NewXMatrix[:, 2]+ w_coef[4]*NewXMatrix[:, 4]+ \
           w_coef[6]*NewXMatrix[:, 6]+ w_coef[7]*NewXMatrix[:, 7]+ w_coef[8]*NewXMatrix[:, 8]+ \
           w_coef[9]*NewXMatrix[:, 9]+ w_coef[12]*NewXMatrix[:, 12]+  w_intercept
#%%
errorM= np.zeros((n_predict.shape[0], 4))
errorM[:, 0]= XMatrix[:, 0]
errorM[:, 1]= XMatrix[:, 1]
errorM[:, 2]= abs(XMatrix[:, 0]-n_predict)/XMatrix[:, 0]
errorM[:, 3]= abs(XMatrix[:, 1]-w_predict)/XMatrix[:, 1]
#%%
errorM= np.round(errorM, 4)
#%%
plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 2], s= 500, cmap= 'nipy_spectral')
plt.colorbar()
plt.xscale('log')
plt.title('Density Error (9 terms)')
plt.xlabel('True Density (1e16 cm$^{-3}$)')
plt.ylabel('True Plasma Width ($\mu$m)')
#%%
plt.scatter(errorM[:, 0], errorM[:, 1], c=errorM[:, 3], s= 500, cmap= 'nipy_spectral')
plt.colorbar()
plt.xscale('log')
plt.title('Width Error (1245)')
plt.xlabel('True Density (1e16 cm$^{-3}$)')
plt.ylabel('True Plasma Width ($\mu$m)')

