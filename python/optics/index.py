# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:55:43 2022

@author: Robert
"""

import numpy as np
from scipy.constants import physical_constants
c = physical_constants['speed of light in vacuum'][0]
eps_0 = physical_constants['vacuum electric permittivity'][0]

# Coefficients for the Sellmeier equation, at 20deg C
# Retrieved from the CVI materials catalog and refractiveindex.info
disp_coef = {
    'fused_silica' : {
        'B1' : 0.6961663,
        'B2' : 0.4079426,
        'B3' : 0.8974794,
        'C1' : 0.0046791,
        'C2' : 0.0135121,
        'C3' : 97.9340025,
    },
    'CaF2' : {
        'B1' : 0.5675888,
        'B2' : 0.4710914,
        'B3' : 3.8484723,
        'C1' : 0.00252643,
        'C2' : 0.01007833,
        'C3' : 1200.5560,
    },
    'MgF2' : {
        'B1' : 0.48755108,
        'B2' : 0.39875031,
        'B3' : 2.3120353,
        'C1' : 0.001882178,
        'C2' : 0.008951888,
        'C3' : 566.13559,
    },
    'N-BK7' : {
        'B1' : 1.03961212,
        'B2' : 0.231792344,
        'B3' : 1.01046945,
        'C1' : 0.00600069867,
        'C2' : 0.0200179144,
        'C3' : 103.560653,
    },
}

# Nonlinear index of refraction, in m^2/W
# Convert from esu n_2[m^2/W]=40*pi*n_2[esu]/(c*n)
# An analytic formula exists to approximate n_2 in crystals, Boling 1987
n2 = {
    'fused_silica' : 2.19e-20,
    'CaF2' : 1.71e-20,
    'MgF2' : 0.76e-20,
    'N-BK7' : 3.18e-20,
}

def index_coef(lam, B1, B2, B3, C1, C2, C3):
    """ Calculate the index of refraction of using the Sellmeier equation.

    Parameters
    ----------
    lam : float
        The wavelength of the light in vacuum.
    B1 : float
        Sellmeier coefficient B1.
    B2 : float
        Sellmeier coefficient B2.
    B3 : float
        Sellmeier coefficient B3.
    C1 : float
        Sellmeier coefficient C1.
    C2 : float
        Sellmeier coefficient C2.
    C3 : float
        Sellmeier coefficient C3.

    Returns
    -------
    n : float
        Index of refraction of the material.

    """
    lam_um = lam*1e6
    n2 = 1+B1*lam_um**2/(lam_um**2-C1)+B2*lam_um**2/(lam_um**2-C2)+B3*lam_um**2/(lam_um**2-C3)
    return np.sqrt(n2)

def index_coef_f(f, B1, B2, B3, C1, C2, C3):
    """ Calculate the index of refraction using the Sellmeier equation.

    Parameters
    ----------
    f : float
        The frequency of the light in vacuum.
    B1 : float
        Sellmeier coefficient B1.
    B2 : float
        Sellmeier coefficient B2.
    B3 : float
        Sellmeier coefficient B3.
    C1 : float
        Sellmeier coefficient C1.
    C2 : float
        Sellmeier coefficient C2.
    C3 : float
        Sellmeier coefficient C3.

    Returns
    -------
    n : float
        Index of refraction of the material.

    """
    f = f*1e-6/c
    n2 = np.array(1+B1/(1-C1*f**2)+B2/(1-C2*f**2)+B3/(1-C3*f**2), dtype='complex128')
    return np.sqrt(n2)

def index_material(lam, mat):
    """ Calculate the index of refraction of using the Sellmeier equation.

    Parameters
    ----------
    lam : float
        The wavelength of the light in vacuum.
    mat : string
        Name of the material in the disp_coef dictionary defined in this module.

    Returns
    -------
    n : float
        Index of refraction of the material.

    """
    C = disp_coef[mat]
    return index_coef(lam, C['B1'], C['B2'], C['B3'], C['C1'], C['C2'], C['C3'])

def index_material_f(f, mat):
    """ Calculate the index of refraction of using the Sellmeier equation.

    Parameters
    ----------
    f : float
        The frequency of the light in vacuum.
    mat : string
        Name of the material in the disp_coef dictionary defined in this module.

    Returns
    -------
    n : float
        Index of refraction of the material.

    """
    C = disp_coef[mat]
    return index_coef_f(f, C['B1'], C['B2'], C['B3'], C['C1'], C['C2'], C['C3'])

def group_vg_f(f, B1, B2, B3, C1, C2, C3):
    """ Calculate the group velocity using the Sellmeier equation.

    Parameters
    ----------
    f : float
        The frequency of the light in vacuum.
    B1 : float
        Sellmeier coefficient B1.
    B2 : float
        Sellmeier coefficient B2.
    B3 : float
        Sellmeier coefficient B3.
    C1 : float
        Sellmeier coefficient C1.
    C2 : float
        Sellmeier coefficient C2.
    C3 : float
        Sellmeier coefficient C3.

    Returns
    -------
    vg : float
        Group velocity in the material.
    """
    n = index_coef_f(f, B1, B2, B3, C1, C2, C3)
    # Convert Ci coefficients to m
    C1 *= 1e-12
    C2 *= 1e-12
    C3 *= 1e-12
    vg = c/n - f**2/(c*n**3)*(B1*C1/(1-C1*f**2/c**2)**2+B2*C2/(1-C2*f**2/c**2)**2+B3*C3/(1-C3*f**2/c**2)**2)
    return vg

def vg_material_f(f, mat):
    """ Calculate the group velocity using the Sellmeier equation.

    Parameters
    ----------
    f : float
        The frequency of the light in vacuum.
    mat : string
        Name of the material in the disp_coef dictionary defined in this module.

    Returns
    -------
    vg : float
        Group velocity in the material.

    """
    C = disp_coef[mat]
    return group_vg_f(f, C['B1'], C['B2'], C['B3'], C['C1'], C['C2'], C['C3'])
    