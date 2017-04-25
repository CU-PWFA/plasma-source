# -*- coding: utf-8 -*-
"""
@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt
#import ADK_Combined as adk

def GaussianBeamStatic(E0,k,w0,z0,z,r):
    zR=Gauss_zR(w0,2*np.pi/k)
    w=Gauss_w(z,z0,zR,w0)
    alpha=Gauss_alpha(z,z0,zR)
    R=Gauss_R(z,z0,zR)
    return 1j*E0/zR*np.exp(k*z*1j)*w0/w* \
        np.exp(-np.power(r/w,2))* \
        np.exp(1j*k*np.power(r,2)/2/R)* \
        np.exp(-1j*alpha)

def GaussianBeamTime(E0,k,w0,z0,w,z,r,t):
    static = GaussianBeamStatic(E0,k,w0,z0,z,r)
    E = static*np.exp(-1j*w*t)
    return E

def GaussianBeamIntensity(I0,k,w0,z0,z,r):
    zR=Gauss_zR(w0,2*np.pi/k)
    w=Gauss_w(z,z0,zR,w0)
    return I0*np.power(w0/w,2)*np.exp(-2*np.power(r/w,2))

def Gauss_w(z,z0,zR,w0):
    return w0*np.sqrt(1+np.power((z-z0)/zR,2))

def Gauss_R(z,z0,zR):
    return z-z0+np.power(zR,2)/(z-z0)

def Gauss_zR(w0,lmba):
    return np.pi*np.power(w0,2)/lmba

def Gauss_alpha(z,z0,zR):
    return np.arctan((z-z0)/zR)



def Prop_Evolve(A,B,C,D,q_in):
    inv_q = (C+D/q_in)/(A+B/q_in)
    return 1/inv_q

def Prop_FreeSpace(q,d,step):
    x = step
    while (x <= d):
        q.append(Prop_Evolve(1,step,0,1,q[-1]))
        x=x+step
    return q

def Prop_ThinLens(q,f):
    q[-1]=Prop_Evolve(1,0,-1/f,1,q[-1])
    return q

def Prop_FlatInterface(q,n1,n2):
    q[-1]=Prop_Evolve(1,0,0,n1/n2,q[-1])
    return q
        
def Prop_RadiusList(q):
    p=list(q)
    for i in range(len(q)):
        x=1/q[i]
        p[i]=1/(x.real)
    return p

def Prop_SpotList(q,lmba,n):
    m=list(q)
    for i in range(len(q)):
        x=1/q[i]
        x=-1/(x.imag)
        m[i]=np.sqrt(lmba*x/np.pi/n)
    return m

"""
size = 500
domain = np.arange(-size,size)
counter = domain+size
r=domain/size*1e-4
z=(domain/size)*2e-1
z0 = 0
I0 = 1

wavelength = 100e-9
k = np.pi*2/(wavelength)
w0 = np.sqrt(0.159*wavelength/2)*0.4
zR = np.power(w0,2)*np.pi/wavelength

print("Waist:",w0)
print("Rayleigh:",Gauss_zR(w0,wavelength))
print("Wavelength:",wavelength)

I=np.empty([size*2,size*2])
for i in counter:
    for j in counter:
        I[i][j]=GaussianBeamIntensity(I0,k,w0,z0,z[i],r[j])
"""
"""
H=np.empty([size,size])
for i in domain:
    H[i] = adk.f(adk.Elefield(I[i]*1.4e14),13.6,1)
"""
"""
plt.imshow(I, interpolation="none", origin="lower",
           extent=[r[0],r[2*size-1],z[0],z[2*size-1]], aspect='.001')
CB=plt.colorbar()
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()

plt.figure()
levels = np.arange(0.05,0.95,0.1)
CS=plt.contour(r*10000,z,I,levels)
CB=plt.colorbar(CS,extend='both')
CB.set_label('I/I0')
plt.ylabel('z (m)')
plt.xlabel('r (10^-4 m)')
plt.title('Relative Intensity of a Gaussian Beam')
plt.show()
"""
"""
plt.imshow(H, interpolation="none", origin="lower",
           extent=[r[0],r[size-1],z[0],z[size-1]], aspect='.001')
plt.colorbar()
plt.show()
"""