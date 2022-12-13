# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 23:08:31 2022

@author: Anubhav
"""
from scipy import signal
import numpy as np
# filter_aplha = np.float32(signal.firwin(400, [0.0004, 0.0010], 
#                                   pass_zero=False))
# filter_beta = np.float32(signal.firwin(400, [0.0010, 0.001], 
#                                   pass_zero=False))

# ii_filt_ipsi = np.float32(signal.convolve(ii_raw_ipsi, filter, 
#                                        mode='same'))


def filtering(y1,y2,sampling_freq=10000):
    
    band_input = int(input("Which Band would you like to filter ?\nFor delta (1 – 4 Hz) Press 1\nFor theta (4 – 10 Hz) Press 2\nFor alpha (10 – 13 Hz) Press 3\nFor beta (13 – 30 Hz) Press 4\nFor gamma (30 – 100 Hz) Press 5\n"))
    if band_input ==1:
        frequency_band =  np.array([1,4])
    elif band_input ==2:
        frequency_band =  np.array([4,10])
        # For our case we used 4,10 as theta
    elif band_input ==3:
        frequency_band =  np.array([10,13])
    elif band_input ==4:
        frequency_band =  np.array([13,30])
    elif band_input ==5:
        frequency_band =  np.array([30,100])
    else:
        print("Invalid Input")
        
    filter = np.float32(signal.firwin(400, [frequency_band[0]/sampling_freq,frequency_band[1]/sampling_freq], 
                                      pass_zero=False))
    
    # First we convolve then we convert to float32 because 
    #float64 hangs the prog.
    
    filt_ii_ipsi = np.float32(signal.convolve(y1, filter, 
                                           mode='same'))
    filt_ii_contra = np.float32(signal.convolve(y2, 
                                             filter, mode='same'))
    return filt_ii_ipsi,filt_ii_contra