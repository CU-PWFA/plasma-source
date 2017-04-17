#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 09:25:20 2017

@author: litos
"""

import numpy as np
from calc_Bmag import calc_Bmag

def calc_M(Tb,Tl):
    """ Calculates mismatch parameter M.
    
        Bmag measures matching of a beam to a lattice.
    
        Parameters
        ----------
        Tb
            Twiss parameter array for beam:
                Tb[0] = beta (m)
                Tb[1] = alpha
                Tb[2] = gamma (1/m)
        Tl
            Twiss parameters for lattice:
                Tl[0] = beta (m)
                Tl[1] = alpha
                Tl[2] = gamma (1/m)

        Returns
        -------
        Bmag
            matching parameter >= 1
                =1 -> matched beam
                >1 -> mismatched beam
    """
    Bmag = calc_Bmag(Tb,Tl)
    
    M = Bmag + np.sqrt(Bmag**2 - 1)

    return M
