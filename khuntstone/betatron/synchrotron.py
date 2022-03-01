# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, m_e, mu_0


from scipy.special import kv

import thomas as th
inv_c = 1 / c
def analytic_synch(w, gamma, p, angle):
    """
    Function to compute the analytic solution for bending magnet radiation.
    
    Parameters:
    -----------
    w     : array_like
            Frequency array rad/s
    gamma : float
            Electron lorentz factor
    p     : float
            Bending radius, m
    angle : float
            x-z angle of observation
    """
    
    t1 = (mu_0 * e**2 * c * w**2 / (12 * np.pi**3)) * (p / c)**2
    t2 = ((1 / gamma**2) + angle**2)**2
    xi = w * p / (3 * c * gamma**3 * (1 + angle**2 * gamma**2)**(1.5))
    t3 = kv(2/3, xi)**2 + ((angle**2 * kv(1/3, xi)**2) / ((1 / gamma**2) \
                           + angle**2))
    return t1 * t2 * t3

# Initial parameters from paper
Nw = int(1e3)
w0    = 2.5e15
w     = np.linspace(w0, w0*1e7, Nw)

gamma = 1000
gamma2 = gamma*gamma
inv_gamma2 = 1 / gamma2

angle    = 0

beta  = np.sqrt(1 - 1/gamma2)
B     = m_e * w0 / e
p     = gamma * m_e * c / (e * B)
# Time step
dtau  = np.pi*1e-4 / w0
dt   = gamma * dtau

# Time for 1 rotation
t_rot = (2 * np.pi * p) / (0.999 * c)
# Useful variables
inv_pi = 1 / np.pi
inv_p = 1/p

t    = np.arange(0, 0.1*t_rot, dt)
Nt   = 10000
t    = np.linspace(-0.5*Nt*dt, 0.5*Nt*dt, Nt)
tau  = t / gamma
# Particle position
x    = p * np.sin(beta * c * t * inv_p)
y    = p * (1 - np.cos(beta * c * t * inv_p))
z    = np.zeros(len(x))
x4   = np.array([c*t, x, y, z], dtype = "longdouble")

# Particle velocity
v00 = th.get_dtau(tau, x)[0]
v01 = th.get_dtau(tau, y)[0]
v02 = th.get_dtau(tau, z)[0]
v0  = inv_c*np.array([v00, v01, v02], dtype = "longdouble")

# Interpolation of position 
# ct, have to do manually
ct1   = np.zeros(len(tau)) + gamma * c
ct2   = np.zeros(len(tau))
# x
x1, x2 = th.get_dtau(tau, x4[1])
x2      = 0.5 * x2
# y 
y1, y2 = th.get_dtau(tau, x4[2])
y2      = 0.5 * y2
# z
z1, z2 = th.get_dtau(tau, x4[3])
z2      = 0.5 * z2


x41   = np.array([ct1, x1, y1, z1], dtype = "longdouble")
x42   = np.array([ct2, x2, y2, z2], dtype = "longdouble")

# Interpolation of velocity 
# v_ct
vct1 = np.zeros(len(tau))
# vx
vx1  = th.get_dtau(tau, v0[0])[0]
# vy
vy1  = th.get_dtau(tau, v0[1])[0]
# vz
vz1  = th.get_dtau(tau, v0[2])[0]
v1   = inv_c * np.array([vx1, vy1, vz1], dtype = "longdouble")

# Computation of on-axis radiation
# Set observation vector and angle
s_hat = np.array([1, 0, 0], dtype = "longdouble") # on-axis
theta = 90 * np.pi / 180

# Parameter for switching between Frensel and Taylor
small     = 1e-3
# Compute radiation
d2I_comp  = th.get_d2I(w, s_hat, gamma, x4, x41, x42, v0, v1, dtau,\
                       theta, small)
    
# Plot and compare
w_ana    = np.linspace(w0, w0*1e7, int(1e6))
d2I_ana  = analytic_synch(w_ana, gamma, p, angle)
fig = plt.figure(figsize = (8, 4))
ax1 = fig.add_subplot(121)
ax1.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$ $\left[J \cdot rad^{-1}' \
               + r'$\cdot s \cdot sr^{-1}\right]$')
ax1.set_xlabel(r'$\omega/\omega_0$')
ax1.set_title('Simulated')
ax1.semilogx(w/w0, np.real(d2I_comp), '-')
ax2 = fig.add_subplot(122)
ax2.semilogx(w_ana/w0, np.real(d2I_ana))
ax2.set_xlabel(r'$\omega/\omega_0$')
ax2.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$ $\left[J \cdot rad^{-1} '\
               + r'$\cdot s \cdot sr^{-1}\right]$')
ax2.set_title("Analytic")
ax1.set_xlim([1e0, 1e7])
ax2.set_xlim([1e0, 1e7])
plt.show()
fig2 = plt.figure(figsize = (8, 4), dpi=200)
ax3  = fig2.gca()
ax3.semilogx(w/w0, np.real(d2I_comp)/max(np.real(d2I_comp)), \
                           label = "Computed")
ax3.semilogx(w_ana/w0, np.real(d2I_ana)/max(np.real(d2I_ana)), '--', \
                               label = "Analytic")
ax3.legend()
ax3.set_xlabel(r'$\omega /\omega_0$')
ax3.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$ $\left[J \cdot rad^{-1} '\
               +'$\cdot s \cdot sr^{-1}\right]$')
ax3.set_title("Normalized Comparison")
plt.show()

