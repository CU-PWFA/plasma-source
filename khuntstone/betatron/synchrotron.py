# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, e, m_e, mu_0
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.special import kv

import thomas as th


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

# Functions for interpolating 
def get_ders(tau, var):
    s = ius(tau, var)
    xp = s.derivative(n = 1)
    xpp = s.derivative(n = 2)
    return xp(tau), xpp(tau)

def quad_interp(var, var1n, var2n, tau, tau_int):
    var_int = np.zeros(len(tau_int))
    for i in range(len(tau_int)):
        ind = np.argmin(abs(tau_int[i] - tau))
        taun = tau[ind]
        var_int[i] = var[ind]\
                     + var1n[ind] * (tau_int[i] - taun) \
                     + var2n[ind] * (tau_int[i] - taun)**2
    return var_int

def lin_interp(var, var1n, tau, tau_int):
    var_int = np.zeros(len(tau_int))
    for i in range(len(tau_int)):
        ind  = np.argmin(abs(tau_int[i] - tau))
        taun = tau[ind]
        var_int[i] = var[ind] + var1n[ind] * (tau_int[i] - taun)
    return var_int

# Thomas Algorithm implemented below
    
# Initial parameters from paper
w0    = 2.5e15
w     = np.linspace(1e0, 1e7, 10000) * w0
B     = m_e * w0 / e
dtau  = np.pi*1e-4 / w0
gamma = 1000
beta  = np.sqrt(1 - 1/gamma**2)
p     = gamma * m_e * c / (e * B)
angle    = 0
beta = np.sqrt(1 - (1/gamma**2))
dt   = gamma * dtau
N    = 1000
t    = np.linspace(-N * dt, N * dt, N*2)
tau  = t / gamma

# Particle position
x    = p * np.sin(beta * c * t / p)
y    = p * (1 - np.cos(beta * c * t / p))
z    = np.zeros(len(x))
x4   = (c*t, x, y, z)

# Particle velocity
v40 = np.zeros(len(tau)) + c * gamma
v41 = get_ders(tau, x)[0]
v42 = get_ders(tau, y)[0]
v43 = get_ders(tau, z)[0]
v4  = (v40, v41, v42, v43)

# Interpolation of position 
tau_int = np.linspace(tau[0], tau[-1], 10000)
# ct, have to do manually
ct1n   = np.zeros(len(tau)) + gamma * c
ct2n   = np.zeros(len(tau))
ct_int = quad_interp(c*t, ct1n, ct2n, tau, tau_int)
# x
x1n, x2n = get_ders(tau, x4[1])
x2n      = 0.5 * x2n
x_int    = quad_interp(x, x1n, x2n, tau, tau_int)
# y 
y1n, y2n = get_ders(tau, x4[2])
y2n      = 0.5 * y2n
y_int    = quad_interp(y, y1n, y2n, tau, tau_int)
# z
z1n, z2n = get_ders(tau, x4[3])
z2n      = 0.5 * z2n
z_int    = quad_interp(z, z1n, z2n, tau, tau_int)

x41n   = (ct1n, x1n, y1n, z1n)
x42n   = (ct2n, x2n, y2n, z2n)
x4_int = (ct_int, x_int, y_int, z_int)

# Interpolation of velocity 
# v_ct
vct1n   = np.zeros(len(tau))
vct_int = lin_interp(v40, vct1n, tau, tau_int)
# vx
vx1n = get_ders(tau, v4[1])[0]
vx_int = lin_interp(v41, vx1n, tau, tau_int)
# vy
vy1n = get_ders(tau, v4[2])[0]
vy_int = lin_interp(v42, vy1n, tau, tau_int)
# vz
vz1n = get_ders(tau, v4[3])[0]
vz_int = lin_interp(v43, vz1n, tau, tau_int)

v4_int = (vct_int, vx_int, vy_int, vz_int)
v1n   = (vx1n, vy1n, vz1n)

# Computation of on-axis radiation
# Set observation vector and angle
s_hat = [1, 0, 0] # on-axis
theta = 90 * np.pi / 180
# Compute radiation (computationally and analytically)
d2I_comp = th.get_d2I(w, s_hat, x4, x41n, x42n, v4, v1n, dtau, theta)
d2I_ana  = analytic_synch(w, gamma, p, angle)
# Plot & Compare
fig1 = plt.figure()
ax1  = fig1.add_subplot(121)
ax2  = fig1.add_subplot(122)
ax1.set_xlabel(r'$\omega / \omega_0$')
ax1.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax1.set_title("Thomas Algorithm")
ax2.set_xlabel(r'$\omega / \omega_0$')
ax2.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax2.set_title("Analytical Solution")
ax1.semilogx(w / w0, d2I_comp)
ax2.semilogx(w / w0, d2I_ana)
plt.show()

fig2 = plt.figure()
ax3  = fig2.gca()
ax3.set_xlabel(r'$\omega / \omega_0$')
ax3.set_ylabel(r'$\frac{d^2I}{d\omega d\Omega}$')
ax3.semilogx(w / w0, d2I_comp / max(d2I_comp), label = "Computational")
ax3.semilogx(w / w0, d2I_ana / max(d2I_ana), label = "Analytic")
ax3.legend()
plt.show()