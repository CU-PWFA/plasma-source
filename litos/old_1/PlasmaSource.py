#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:58:21 2017

@author: mike
"""

class PlasmaSource(object):
    """Entire PWFA Plasma Source"""
    def __init__(self):
        self.main = PlasmaMain()
        self.lens = PlasmaLens()
        
class PlasmaMain(object):
    """Main PWFA Plasma Source"""
    def __init__(self):
        self.up_ramp = UpRamp()
        self.dn_ramp = DnRamp()
        self.bulk    = Bulk()

class UpRamp(object):
    """Up-Ramp Part of Main Plasma Source"""
    def __init__(self):
        self.shape = 0
        self.npl0  = 0
        self.hwhm  = 0

class DnRamp(object):
    """Down-Ramp Part of Main Plasma Source"""
    def __init__(self):
        self.shape = 0
        self.npl0  = 0
        self.hwhm  = 0

class Bulk(object):
    """Bulk Flat-Top Part of Main Plasma Source"""
    def __init__(self):
        self.npl0  = 0
        self.Lp    = 0

class PlasmaLens(object):
    """Thin Plasma Lens Part of PWFA Plasma Source"""
    def __init__(self):
        self.npl0 = 0
        self.Lp   = 0