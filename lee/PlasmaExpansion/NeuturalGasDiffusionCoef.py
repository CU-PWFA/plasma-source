# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
x= np.linspace(-1000e-6, 1000e-6, 1000)
y= np.linspace(-1000e-6, 1000e-6, 1000)
X, Y= np.meshgrid(x, y)
r=np.sqrt(X**2+Y**2)

n0= 1e23 #m^-3
n= n0-n0*np.exp(-r**2/(200e-6)**2)
#%%
d= 0.49e-10 #m
sigma= d*1e10 #Angstrom
omega= 1
P= 0.0041 #atm = 4.2mbar= 1e17cm^-3
M= 4 #g/mol
T= 293
A= 1.859e-3 #weird unit see notebook
D= A*T**(3/2)*np.sqrt(2/M)/(P*sigma**2*omega) #cm^2/s
#%%
kb=1.38064852e-23  # J/K 
m= 6.646e-27
l=1/(np.sqrt(2)*np.pi*d**2*n0)
v=np.sqrt(3*kb*T/m)
D1= 1/3*(l*v)
#%%
D2= 3/(8*n0*d**2)*(np.sqrt(kb*T/(2*np.pi)*(2/m)))
#%%
plt.pcolormesh(D)
plt.axis("scaled")
plt.colorbar()
#%%Collision frequency
cross= 5e-15 #(cm^2)
#kb=1.6e-12  # erg/eV
kb=1.38064852e-23  # J/K 
n0= 1e23 #m^-3
EV = 11604.5
Ti = 0.025*EV
Te = 2*EV
#me= 9.11e-28
#mi=4*1.67e-24
me= 9.11e-31
mi=4*1.67e-27
e= 1.6e-19
epsilon=8.85e-12
f_He_e= n0*cross*np.sqrt(kb*Te/me)
f_He_i= n0*cross*np.sqrt(kb*Ti/mi)

eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me


print(f_He_e)
print(f_He_i)
print(f_e_i)
#%%Collision frequency OURS
cross= 5e-19 #(m^2)
kb=1.38064852e-23  # J/K 
n0= 1e23 #m^-3
EV = 11604.5
Ti = 0.025*EV
Te = 2*EV
me= 9.11e-31
mi=4*1.67e-27
e= 1.6e-19
epsilon=8.85e-12
f_He_e= n0*cross*np.sqrt(kb*Te/me)
f_He_i= n0*cross*np.sqrt(kb*Ti/mi)
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me
f_e_e= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(me)*(kb*Te)**1.5))*10
f_i_i= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(mi)*(kb*Ti)**1.5))*10

print("{:.2e}".format(f_He_e))
print("{:.2e}".format(f_He_i))
print("{:.2e}".format(f_e_i))
print("{:.2e}".format(f_e_e))
print("{:.2e}".format(f_i_i))

#%%Collision frequency HOFI
cross= 5e-19 #(m^2)
kb=1.38064852e-23  # J/K 
n0= 1e24 #m^-3
EV = 11604.5
Ti = 0.025*EV
Te = 13.7*EV
me= 9.11e-31
mi=1*1.67e-27
e= 1.6e-19
epsilon=8.85e-12
f_H_e= n0*cross*np.sqrt(kb*Te/me)
f_H_i= n0*cross*np.sqrt(kb*Ti/mi)
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me
f_e_e= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(me)*(kb*Te)**1.5))*10
f_i_i= (1/(3*np.sqrt(np.pi)))*n0*((e**2/(4*np.pi*epsilon))**2)*(4*np.pi/(np.sqrt(mi)*(kb*Ti)**1.5))*10

print("{:.2e}".format(f_H_e))
print("{:.2e}".format(f_H_i))
print("{:.2e}".format(f_e_i))
print("{:.2e}".format(f_e_e))
print("{:.2e}".format(f_i_i))



#%%Gauss unit
me= 9.11e-28
mi= 4*1.67e-24
Te= 2
Ti= 0.025
#%%
#%%SI unit
me= 9.11e-31
mi= 4*1.67e-27
EV= 11604.5
Te= 2*EV
Ti= 0.025*EV
#%%Collision frequency
cross= 5e-15 #(cm^2)
kb=1.6e-12  # erg/eV
n0= 1e17 #cm^-3
Ti = 0.025
Te = 2
me= 9.11e-28
mi=4*1.67e-24
f_He_e= n0*cross*np.sqrt(kb*Te/me)
f_He_i= n0*cross*np.sqrt(kb*Ti/mi)
print("{:.2e}".format(f_He_e))
print("{:.2e}".format(f_He_i))
#%%
#%%SI unit
kb=1.38064852e-23  # J/K 
me= 9.11e-31 #kg
mi= 4*1.67e-27 #kg
EV= 11604.5
Te= 2*EV
Ti= 0.025*EV
e= 1.6e-19
epsilon=8.85e-12
eta= np.pi*e**2*np.sqrt(me)/((4*np.pi*epsilon)**2*(kb*Te)**(3/2))*10
f_e_i= n0*e**2*eta/me
print("{:.2e}".format(f_e_i))
