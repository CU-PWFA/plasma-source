#Uses the ADK Model of Ionization for H, He, Ar, Ar+2 to find
# ionization levels as a function of laser intensity, approximated
# as a square wave in this case.  Plots ionization curves on a
# single graph and outputs table data for each species

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as Gamma

#Defining constants and conversion factors
c = 2.99e8
eps_0 = 8.854e-12
GV_V = 1e9
m_cm = 1e-2

#Ionization energies [eV] for species
chi_H = 13.6
chi_He = 24.6
chi_Ar = 15.8
chi_Ar2 = 27.6

#Effective principal quantum number
#  chi - ionization energy
#  z   - ionization level
def nn(chi,z):
    return 3.68859*z/np.sqrt(chi)

#We are using a square pulse of width 100fs
delt_t = 100e-15

#Fraction of atoms ionized per second, ADK
#  t   - Electric field strength
#  chi - ionization energy
#  z   - ionization level
def W(t,chi,z):
    return 1.51927e15*np.power(4,nn(chi,z))*chi/nn(chi,z)/Gamma(2*nn(chi,z))* \
                    np.power(20.4927/t*np.power(chi,1.5),2*nn(chi,z)-1)* \
                    np.exp(-6.83089/t*np.power(chi,1.5))

#Total fraction of atoms ionized during full pulse
#  t   - Electric field strength
#  chi - ionization energy
#  z   - ionization level
def f(t,chi,z):
    return 1.0-np.exp(-W(t,chi,z)*delt_t)

#Convertion from an electric field magnitude to its equivalent 
# intensity.  Intensities returned are in W/cm^2
#  E_field - Electric field strength (V/m)
def Intensity(E_field):
    return np.power(m_cm,2)*.5*c*eps_0*np.power(E_field*GV_V,2)

#Setup the arrays.  t1 is for cosmetics and t2 is more precise
start = 17
end = 117
step = (end-start)/40

t1 = np.arange(start, end, step)
t2 = np.arange(start, end, step/100)

#Given an array of ionization fractions, finds where that array exceeds
# a given value and prints out the intensity and electric field of
# the corresponding ionization fraction.  Also gives information about
# how ionization fractions change with a 1% change in Intensity
#  dist - array of ionization fractions
#  diff - array of %change in io. per %change in I
#  chi  - ionization energy of species
#  x    - cutoff ionization fraction
#  z    - level of ionization
def find_where(dist,diff,chi,x,z):
    top = np.where(dist > x)
    print(x*100,"% Ionized:")
    ele = t2[top][0]
    fra = diff[top][0]
    print("E --> I: ",ele,"  -->  ",Intensity(ele))
    print("% io. per % I: ",fra)
    print("io. range: ",f(ele,chi,z)*(1-fra),"<--->",f(ele,chi,z)*(1+fra))
    print()
    return

#Prints out information on intensities needed for various values
# of ionization fractions, and information about how ionization
# fractions change with a 1% change in Intensity
#  chi - ionization energy of species
#  z   - levle of ionization
def get_Table(chi,z):
    print("==============")
    print("Values for",chi,"eV")
    dist = f(t2,chi,z)
    diff = (f(t2+t2*.00495,chi,z)-f(t2,chi,z))/f(t2,chi,z)
    vals = [.01,.05,.1,.15,.25,.5,.75,.85,.9,.95,.99,.999]
    for x in vals:
        find_where(dist,diff,chi,x,z)
    print("==============")
    return 

#Get all the information on all the species
get_Table(chi_H,1)
get_Table(chi_He,1)
get_Table(chi_Ar,1)
get_Table(chi_Ar2,2)

#Plots the ionization curves for all species
plt.figure(1)
plt.plot(Intensity(t1), f(t1,chi_H,1), 'bo', Intensity(t2), f(t2,chi_H,1), 'k',
         Intensity(t1), f(t1,chi_He,1), 'ro', Intensity(t2), f(t2,chi_He,1), 'k',
         Intensity(t1), f(t1,chi_Ar,1), 'go', Intensity(t2), f(t2,chi_Ar,1), 'k',
         Intensity(t1), f(t1,chi_Ar2,2), 'yo', Intensity(t2), f(t2,chi_Ar2,2), 'k')
plt.title("ADK Ionization of H(Blue), He(Red), Ar(Green 1st, Yellow 2nd)")
plt.ylabel("Ionization Fraction")
plt.xscale('log')
plt.xlabel("Laser Intensity (W/cm^2)")
plt.show()