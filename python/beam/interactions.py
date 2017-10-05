#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 14:56:11 2017

@author: robert
"""

import numpy as np
import beam.beams as beams
import beam.elements as elements
import beam.calc.plasma as pcalc


def pulse_uniform_gas(pulse, plasma):
    """ Propagates a pulse through a gas, ionizing and refracting as it goes.
    
    Parameters
    ----------
    pulse : Pulse class
        The laser pulse to propagate through the plasma.
    plasma : Plasma class
        The gas to propagate the laser pulse through, must be uniform.
    """
    
    pulse.e = pcalc.plasma_refraction(pulse.e, pulse.x, pulse.y,
                      plasma.z, pulse.t, pulse.lam, plasma.n0, pulse.z[-1],
                      pulse.fft, pulse.ifft, pulse.save_field, 
                      plasma.save_density, plasma.atom)
    