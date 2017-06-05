# -*- coding: utf-8 -*-
"""
Contains functions for building a Gaussian Wave, obtaining an analyzing its
parameters (such as Rayleigh length, waist size, etc).  Also contians functions
for propagating a gaussian wave through optical systems using ABCD matrices and
the complex q parameter

@author: chris
"""

import numpy as np

#Calculates a  gaussian beam in terms of intensity  (most useful)
#  I_0 - Intensity at beam waist
#  k - wavenumber
#  w0 - waist size of beam
#  z0 - location of beam waist
#  z - Location along propagation axis for calculation
#  r - location transverse to beam for calculation
def GaussianBeamIntensity(I0,k,w0,z0,z,r):
    zR=Gauss_zR(w0,2*np.pi/k)
    w=Gauss_w(z,z0,zR,w0)
    return I0*np.power(w0/w,2)*np.exp(-2*np.power(r/w,2))

#Calculates a  gaussian beam in terms of intensity
# this one assumes you are given an array of spot sizes, rather than
# calculating one assuming free space propagation (ie: w obtained from
# using the ABCD propagation functions with the complex q parameter)
#  I_0 - Intensity at the initial beam waist
#  w - array of spot sizes
#  w0 - waist size of beam initially (corresponding to I0)
#  r - location transverse to beam for calculation
def GaussianBeamIntensity_SpotArray(I0,w,w0,r):
    return I0*np.power(w0/w,2)*np.exp(-2*np.power(r/w,2))

#Calculates a  gaussian beam in terms of intensity
# this one assumes you are given an array of spot sizes, rather than
# calculating one assuming free space propagation (ie: w obtained from
# using the ABCD propagation functions with the complex q parameter)
#  I_0 - Intensity at the initial beam waist
#  w - array of spot sizes
#  w0 - waist size of beam initially (corresponding to I0)
#  r - location transverse to beam for calculation
def GaussianBeamIntensity_SpotArray_2D(I0,w_x,w_y,w0,x,y):
    return I0*np.power(w0,2)/(w_x*w_y)*np.exp(-2*np.power(x/w_x,2)-2*np.power(y/w_y,2))

#Calculates the spot size (Transerve spread of beam)
#  z - distance along propagation axis
#  z0 - location of beam waist
#  zR - Rayleigh length
#  w0 - spot size of beam waist
def Gauss_w(z,z0,zR,w0):
    return w0*np.sqrt(1+np.power((z-z0)/zR,2))

#Calculates the radius of curvature of the beam front
#  z - distance along propagation axis
#  z0 - location of beam waist
#  zR - Rayleigh length
def Gauss_R(z,z0,zR):
    return z-z0+np.power(zR,2)/(z-z0)

#Calculates the Rayleigh length (Longitudinal width of beam waist)
#  w0 - spot size of beam waist
#  wavelength - wavelength of beam
def Gauss_zR(w0,wavelength):
    return np.pi*np.power(w0,2)/wavelength

#Calculates the Guoy phase shift
#  z - distance along propagation axis
#  z0 - location of beam waist
#  zR - Rayleigh length
def Gauss_alpha(z,z0,zR):
    return np.arctan((z-z0)/zR)

#Initializes a list of q parameters to be used for propagation
# NOTE: this list is actually a list of three elements:
# complex q parameter, z-coordinate, index of refraction
#  wavelength - wavelength of the beam
#  w0 - beam waist size
#  z0 - starting point relative to location of beam waist
#  n - starting index of refraction
def Prop_Init_q(wavelength,w0,z0,n):
    zR = Gauss_zR(w0, wavelength)
    offset = z0
    #Because R is inverse of z0, I just make sure it is not zero in order
    # to avoid dividing by zero
    if offset == 0:
        offset = 1e-40
    q_init = [complex(offset,zR),z0,n]
    q=[q_init]
    return q

#Propagates the q parameter using an ABCD matrix and returns a new q parameter
#  A,B,C,D - the four elements on 2x2 matrix
#  q_start - the initial q parameter
def Prop_Evolve(A,B,C,D,q_start):
    q_in=q_start[0]
    inv_q = (C+D/q_in)/(A+B/q_in)
    return 1/inv_q

#Propagates q through free space, returning the updated list of q parameters
#  q - list of q paramters
#  d - total propagation distance
#  step - incrememts of distance to caluclate through (example:  total distance
#          is 1 m, but our step can be 1 cm for a total of 100 propagations)
def Prop_FreeSpace(q,d,step):
    x = step
    while (x <= d):
        q_new=[Prop_Evolve(1,step,0,1,q[-1]),q[-1][1]+step,q[-1][2]]
        q.append(q_new)
        x=x+step
    return q

def Prop_Cylindrical_FreeSpace(q_P,q_T,d,step):
    Prop_FreeSpace(q_P,d,step)
    Prop_FreeSpace(q_T,d,step)

#Propagates q through a thin lens, returning list of q paramters with topmost
# q updates due to the lens
#  q - list of q parameters
#  f - focal length of thin lens
def Prop_ThinLens(q,f):
    q[-1][0]=Prop_Evolve(1,0,-1/f,1,q[-1])
    return q

#Propagates q through a cylindrcal thin lens, returning list of q paramters 
# with topmost q updates due to the lens
#  q_P - list of q parameters Parallel to cyclindrical axis
#  q_T - list of q parameters Transverse to cylindrical axis
#  R - radius of curvature in the radial direction
def Prop_CylindricalLens(q_P,q_T,R):
    q_T[-1][0]=Prop_Evolve(1,0,-2/R,1,q_T[-1])

#Propagates q through a flat interface (normal to axis of propagation)
#  q - list of q parameters
#  n2 - new index of refraction
def Prop_FlatInterface(q,n2):
    n1=q[-1][2]
    q[-1][0]=Prop_Evolve(1,0,0,n1/n2,q[-1])
    q[-1][2]=n2
    return q

#Propagates q thorugh a spherically curved interface
#  q - list of q parameters
#  n2 - new index of refraction
#  R - radius of curvature of the interface
def Prop_CurvedInterface(q,n2,R):
    n1=q[-1][2]
    q[-1][0]=Prop_Evolve(1,0,(n1-n2)/(R*n2),n1/n2,q[-1])
    q[-1][2]=n2
    return q

#Propagates q thorugh a cylindrically curved interface
#  q_P - list of q parameters Parallel to cyclindrical axis
#  q_T - list of q parameters Transverse to cylindrical axis
#  n2 - new index of refraction
#  R - radius of curvature of the interface
def Prop_CylindricalInterface(q_P,q_T,n2,R):
    Prop_FlatInterface(q_P,n2)
    Prop_CurvedInterface(q_T,n2,R)

#Propagates q thorugh two curved interfaces and free space within
# (for a convex lens, R1 must be positive and R2 must be negative)
#  q - list of q parameters
#  n2 - index of refraction within lens
#  d - thickness of lens
#  R1 - radius of curvature of the first interface
#  R2 - radius of curvature of the second interface
#  step - step size of propagation within the lens
def Prop_ThickLens(q,n2,d,R1,R2,step):
    n1=q[-1][2]
    q=Prop_CurvedInterface(q,n2,R1)
    q=Prop_FreeSpace(q,d,step)
    q=Prop_CurvedInterface(q,n1,R2)
    return q

#Propagates the q parameter using a custom ABCD matrix and returns the updated 
# list of q paramters
#  q - list of q parameters
#  A,B,C,D - the four elements on 2x2 matrix
#  offset - length of optical element
#  n_final - index of refraction at the end of optical element
def Prop_Custom(q,A,B,C,D,offset,n_final):
    q_final=Prop_Evolve(A,B,C,D,q[-1])
    if n_final == -1:
        n_final = q[-1][2]
    
    if offset == 0:
        q[-1][0] = q_final
        q[-1][2] = n_final
    else:
        q_new=[q_final,q[-1][1]+offset,n_final]
        q.append(q_new)
    return q

#Returns the list for only the complex q parameters
#  q - list of q parameters
def Prop_GetParameter(q):
    return [param[0] for param in q]

#Returns the list for only the z-coordinates
#  q - list of q parameters
def Prop_GetRange(q):
    return [pos[1] for pos in q]

#Returns the list for only the index of refraction
#  q - list of q parameters
def Prop_GetIndex(q):
    return [index[2] for index in q]

#Returns the list for radius of curvature
#  q - list of q parameters
def Prop_RadiusList(q):
    q_param=Prop_GetParameter(q)
    R=list(q_param)
    for i in range(len(q_param)):
        x=1/q_param[i]
        R[i]=1/(x.real)
    return R

#Returns the list for spot size of beam
#  q - list of q parameters
def Prop_SpotList(q,lmba):
    q_param=Prop_GetParameter(q)
    n=Prop_GetIndex(q)
    W=list(q_param)
    for i in range(len(q_param)):
        x=1/q_param[i]
        x=-1/(x.imag)
        W[i]=np.sqrt(lmba*x/np.pi/n[i])
    return W

#Unused thus far

#Calculates a static Gausssian beam in terms of electric field strength
# (unused)
#  E_0 - Electric Field Strength at waist
#  k - wavenumber
#  w0 - waist size of beam
#  z0 - location of beam waist
#  z - Location along propagation axis for calculation
#  r - location transverse to beam for calculation
def GaussianBeamStatic(E0,k,w0,z0,z,r):
    zR=Gauss_zR(w0,2*np.pi/k)
    w=Gauss_w(z,z0,zR,w0)
    alpha=Gauss_alpha(z,z0,zR)
    R=Gauss_R(z,z0,zR)
    return 1j*E0/zR*np.exp(k*z*1j)*w0/w* \
        np.exp(-np.power(r/w,2))* \
        np.exp(1j*k*np.power(r,2)/2/R)* \
        np.exp(-1j*alpha)

#Caluclates a gaussian beam with time dependency in terms of electric field strength
# (unused)
#  E_0 - Electric Field Strength at waist
#  k - wavenumber
#  w0 - waist size of beam
#  z0 - location of beam waist
#  w - frequency of wave
#  z - Location along propagation axis for calculation
#  r - location transverse to beam for calculation
def GaussianBeamTime(E0,k,w0,z0,w,z,r,t):
    static = GaussianBeamStatic(E0,k,w0,z0,z,r)
    E = static*np.exp(-1j*w*t)
    return E