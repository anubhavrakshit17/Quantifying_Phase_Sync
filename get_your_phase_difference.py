# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 19:10:42 2022

@author: Anubhav
"""

import numpy as np
import scipy.signal as sig

def hilbert_example_phase(y1,y2):
    sig1_hill=sig.hilbert(y1)
    sig2_hill=sig.hilbert(y2)
    pdt=(np.inner(sig1_hill,np.conj(sig2_hill))/(np.sqrt(np.inner(sig1_hill,
               np.conj(sig1_hill))*np.inner(sig2_hill,np.conj(sig2_hill)))))
    phase = np.angle(pdt)

    return phase
def hilbert_phase(y1,y2):
    sig1_hill=sig.hilbert(y1)
    sig2_hill=sig.hilbert(y2)
    phase = np.angle(sig2_hill-sig1_hill)
    
    return phase

    