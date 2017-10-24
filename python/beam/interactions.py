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


def pulse_plasma(pulse, plasma):
    """ Propagates a pulse through a gas, ionizing and refracting as it goes.
    
    Parameters
    ----------
    pulse : Pulse class
        The laser pulse to propagate through the plasma.
    plasma : Plasma class
        The gas to propagate the laser pulse through.
    """
    pulse.e = pcalc.plasma_refraction(pulse.e, pulse.x, pulse.y,
                      plasma.z, pulse.t, pulse.lam, plasma.n0, pulse.z[-1],
                      pulse.fft, pulse.ifft, pulse.save_field, 
                      plasma.save_plasma_density, plasma.atom, 
                      plasma.load_num_den, plasma.load_plasma_den)


def beam_phase(beam, phase):
    """ Applies a phase mask to a optical beam, either a pulse or laser.
    
    Parameters
    ----------
    beam : Pulse or Laser class
        The optical beam the apply the phase mask to.
    phase : Phase class
        The phase mask to apply to the beam.
    """
    beam.set_field(beam.e * np.exp(1j*phase.phi))
    