"""
Contains functions for analyzing the ionization levels of
a uniform gas with an applied electric field.  Currently 
only supports sqaure wave pulses.  Also contains functions
for converting to/from electric fields and intensities

@author: chris
"""

import numpy as np
from scipy.special import gamma as Gamma

#Defining constants and conversion factors
c = 2.99e8
eps_0 = 8.854e-12
GV_V = 1e9
m_cm = 1e-2

#Ionization energies [eV] for species
def Get_chi_H():
    return 13.6
def Get_chi_He():
    return 24.6
def Get_chi_Ar1():
    return 15.8
def Get_chi_Ar2():
    return 27.6

#Effective principal quantum number
#  chi - ionization energy
#  z   - ionization level
def nn(chi,z):
    return 3.68859*z/np.sqrt(chi)

#Fraction of atoms ionized per second, ADK
#  t   - Electric field strength
#  chi - ionization energy
#  z   - ionization level
def DiffProb(t,chi,z):
    return 1.51927e15*np.power(4,nn(chi,z))*chi/nn(chi,z)/Gamma(2*nn(chi,z))* \
                    np.power(20.4927/t*np.power(chi,1.5),2*nn(chi,z)-1)* \
                    np.exp(-6.83089/t*np.power(chi,1.5))

#Total fraction of atoms ionized during full pulse
#  e_fields - array of electric field strengths
#  delt_t - time interval of square pulse
#  t   - Electric field strength
#  chi - ionization energy
#  z   - ionization level
def TotalProb(e_fields,delt_t,chi,z):
    return 1.0-np.exp(-DiffProb(e_fields,chi,z)*delt_t)

#Convertion from/to an electric field magnitude to its equivalent 
# intensity.  Intensities returned are in W/cm^2
#  E_field - Electric field strength (V/m)
#  Intensity - Laser Intensity in W/cm^2
def Intensity(E_field):
    return np.power(m_cm,2)*.5*c*eps_0*np.power(E_field*GV_V,2)
def Elefield(Intensity):
    return 1e9 * 27.4492*np.sqrt(Intensity/1e14)
    #return 1/GV_V/m_cm*np.sqrt(Intensity*2/c/eps_0)

#Prints out information on intensities needed for various values
# of ionization fractions, and information about how ionization
# fractions change with a 1% change in Intensity
#  e_fields - array of electric field strengths
#  delt_t - time interval of square pulse
#  chi - ionization energy of species
#  z   - levle of ionization
def get_Table(e_fields,delt_t,chi,z):
    print("==============")
    print("Values for",chi,"eV")
    dist = TotalProb(e_fields,delt_t,chi,z)
    diff = (TotalProb(e_fields+e_fields*.00495,delt_t,chi,z)- \
            TotalProb(e_fields,delt_t,chi,z))/TotalProb(e_fields,delt_t,chi,z)
    vals = [.01,.05,.1,.15,.25,.5,.75,.85,.9,.95,.99,.999]
    for x in vals:
        find_where(e_fields,delt_t,dist,diff,chi,x,z)
    print("==============")
    return 

#Given an array of ionization fractions, finds where that array exceeds
# a given value and prints out the intensity and electric field of
# the corresponding ionization fraction.  Also gives information about
# how ionization fractions change with a 1% change in Intensity
#  e_fields - array of electric field strengths
#  delt_t - time interval of square pulse
#  dist - array of ionization fractions
#  diff - array of %change in io. per %change in I
#  chi  - ionization energy of species
#  x    - cutoff ionization fraction
#  z    - level of ionization
def find_where(e_fields,delt_t,dist,diff,chi,x,z):
    top = np.where(dist > x)
    print(x*100,"% Ionized:")
    ele = e_fields[top][0]
    fra = diff[top][0]
    print("E --> I: ",ele,"  -->  ",Intensity(ele))
    print("% io. per % I: ",fra)
    print("io. range: ",TotalProb(ele,delt_t,chi,z)*(1-fra),"<--->", \
          TotalProb(ele,delt_t,chi,z)*(1+fra))
    print()
    return