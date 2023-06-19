#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 15:00:29 2022

Magnification calculations from the imaging spectrometer

@author: chris
"""

import numpy as np
import scipy.io as sio

def Rmat(l,k=None):
    ## GENERAL TRANSFER MATRIX (DRIFT AND QUADRUPOLE)
#     [matrix] = Rmat(l,k)
    # drift matrix (if k==0) or quad matrix

    if k is None:
        matrix = np.array([[1,l,0,0],[0,1,0,0],[0,0,1,l],[0,0,0,1]])
    else:
        matrix = np.real(np.array([[np.cos(np.emath.sqrt(k)*l),np.sin(np.emath.sqrt(k)*l)/np.emath.sqrt(k),0,0],
                                    [-np.sin(np.emath.sqrt(k)*l)*np.emath.sqrt(k),np.cos(np.emath.sqrt(k)*l),0,0],
                                    [0,0,np.cosh(np.emath.sqrt(k)*l),np.sinh(np.emath.sqrt(k)*l)/np.emath.sqrt(k)],
                                    [0,0,np.sinh(np.emath.sqrt(k)*l)*np.emath.sqrt(k),np.cosh(np.emath.sqrt(k)*l)]]))
    return matrix
"""
def Rmat(l,k=None):
    if k is None:
        matrix = np.array([[1,l],[0,1]])
    elif k > 0: #Q1
        matrix = np.array([[np.cos(np.sqrt(k)*l),1/np.sqrt(k)*np.sin(np.sqrt(k)*l)],
            [-np.sqrt(k)*np.sin(np.sqrt(k)*l),np.cos(np.sqrt(k)*l)]])
    elif k < 0: #Q0 and Q2
        matrix = np.array([[np.cosh(np.sqrt(-k)*l),1/np.sqrt(-k)*np.sinh(np.sqrt(-k)*l)],
            [np.sqrt(-k)*np.sinh(np.sqrt(-k)*l),np.cosh(np.sqrt(-k)*l)]])
    return matrix
"""
def FACET_Espec_Matrix(B_Vector_IN,Gamma_Beam,Location):

    #Quad lengths
    LQ0=1;
    LQ1=1;
    LQ2=1;

    c0=299792458;
    qE=1.6022e-19;
    mE=9.1094e-31;

    kGauss2Telsa=10**-1;

    B_Vector=np.array([B_Vector_IN[0]/LQ0,B_Vector_IN[1]/LQ1,B_Vector_IN[2]/LQ2])*kGauss2Telsa; # T/m
    K_Vector=B_Vector*qE/(Gamma_Beam*mE*c0);

    #print(" K Vector:");print(K_Vector)

    PosFilament=1993.2737; #position of gas jet.
    PENT=PosFilament+0.5; # 0.5 m DS of fils or so .Just a guess...
    PosQ0D=1996.98249;  #  position of Q0 (center)
    PosQ1D=1999.206615; # position of Q1 (center)
    PosQ2D=2001.431049; # position of Q2 (center)
    PosB5D36=2005.940051; # position of 
    PosDTOTR=2015.2599;  # position of DT OTR screen
    PosLFOV=2015.63; # Gamma 1 is at the same location
    PosPRDump=2017.53; # location of PR dump

    ObjectPlaneLocation = Location
    #ObjectPlaneLocation=PosFilament;
    ImagePlaneLocation=PosDTOTR;

    #Drift lengths
    LD1=(PosQ0D-LQ0/2)-ObjectPlaneLocation; # filament IP to Q0
    LD2=(PosQ1D-LQ1/2)-(PosQ0D+LQ0/2); # Q0 to Q1
    LD3=(PosQ2D-LQ2/2)-(PosQ1D+LQ1/2); # Q2 to Q1
    LD4=ImagePlaneLocation-(PosQ2D+LQ2/2); # Q2 to LFOV screen

    M1=Rmat(LD1); # filament IP to Q0
    #print(" M1: ");print(M1)
    M2=Rmat(LQ0,K_Vector[0]); # Q0
    #print(" M2: ");print(M2)
    M3=Rmat(LD2); # Q0 to Q1
    #print(" M3: ");print(M3)
    M4=Rmat(LQ1,K_Vector[1]); # Q1
    #print(" M4: ");print(M4)
    M5=Rmat(LD3); # Q1 to Q2
    #print(" M5: ");print(M5)
    M6=Rmat(LQ2,K_Vector[2]); # Q2
    #print(" M6: ");print(M6)
    M7=Rmat(LD4); # Q2 to screen
    #print(" M7: ");print(M7)
    MatrixOut = np.matmul(M7,np.matmul(M6,np.matmul(M5,np.matmul(M4,np.matmul(M3,np.matmul(M2,M1))))))
    
    #print(" Result:")
    return MatrixOut

def getSpectrometerArrays(mat_path):
    #Load the mat file
    mat = sio.loadmat(mat_path)
    data = mat['data_struct']
    
    #get the step array
    scalars = data['scalars'][0][0]
    steps = scalars['steps'][0][0]
    steps = np.array(steps[:,0])
    
    #get q01,q02, and q03 arrays
    nonbsa_li20mag = scalars['nonBSA_List_S20Magnets'][0][0]
    q0_vals = nonbsa_li20mag['LI20_LGPS_3141_BACT'][0][0]
    q0_vals = np.array(q0_vals[:,0])
    q1_vals = nonbsa_li20mag['LI20_LGPS_3261_BACT'][0][0]
    q1_vals = np.array(q1_vals[:,0])
    q2_vals = nonbsa_li20mag['LI20_LGPS_3091_BACT'][0][0]
    q2_vals = np.array(q2_vals[:,0])
    
    #convert each q array to an array vs steps
    step_total = max(steps)
    q0_arr = np.zeros(step_total)
    q1_arr = np.zeros(step_total)
    q2_arr = np.zeros(step_total)
    for i in range(len(steps)):
        q0_arr[steps[i]-1] = q0_vals[i]
        q1_arr[steps[i]-1] = q1_vals[i]
        q2_arr[steps[i]-1] = q2_vals[i]
        
    return q0_arr, q1_arr, q2_arr

def getObjectPlaneLocArr(mat_path):
    #Load the mat file
    mat = sio.loadmat(mat_path)
    data = mat['data_struct']
    params = data['params'][0][0]
    startVals = params['startVals'][0][0][0][0]
    stopVals = params['stopVals'][0][0][0][0]
    totalSteps = params['totalSteps'][0][0][0][0]
    
    zobj_arr = np.linspace(startVals,stopVals,totalSteps)
    return zobj_arr

def getM11ForAllSteps(superpath,day,dataset,gamma):
    q0, q1, q2 = getSpectrometerArrays(superpath+day+dataset+'/'+dataset+'.mat')
    zobj_arr = getObjectPlaneLocArr(superpath+day+dataset+'/'+dataset+'.mat')
    
    m_arr = np.zeros(len(q0))
    for i in range(len(m_arr)):
        m_full = FACET_Espec_Matrix([q0[i],q1[i],q2[i]],gamma,zobj_arr[i])
        m_arr[i] = m_full[0,0]
        
    return m_arr

if __name__ == '__main__':
    """
    superpath = '/media/chris/New Volume/SLACData/'
    day = '20220812/'
    dataset = 'E308_02493'
    gamma = 19569.5
    m_arr = getM11ForAllSteps(superpath,day,dataset,gamma)
    print(m_arr)
    """
    
    Q01 = -1.020636901855469e+02
    Q02 = 1.635670928955078e+02
    Q03 = -1.020625762939453e+02
    
    gamma = 19569.5*9.8/10
    loc = 1992
    
    test = FACET_Espec_Matrix([Q01,Q02,Q03],gamma,loc)
    print(test)