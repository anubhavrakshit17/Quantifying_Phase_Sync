# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 19:49:21 2022

@author: Anubhav
"""

import numpy as np
import matplotlib.pyplot as plt
import os

F = 5.e2          # No. of cycles per second, F = 500 Hz
T = 5.e-3         # Time period, T = 2 ms
Fs = 50.e3        # No. of samples per second, Fs = 50 kHz
Ts = 1./Fs        # Sampling interval, Ts = 20 us
N = int(T/Ts)     # No. of samples for 2 ms, N = 100
t1 = np.pi/1
t2 = np.pi/1.5
t3 = np.pi/2
t4 = np.pi/3
t5 = np.pi/4
t = np.linspace(0, T, N)
y1 = np.sin(2*np.pi*F*t)
y2 = np.sin(2*np.pi*F*t+t2)
y3 = np.sin(2*np.pi*F*t+t3)
y4 = np.sin(2*np.pi*F*t+t4)
y5 = np.sin(2*np.pi*F*t+t5)
plt.plot(t, y1, label='Original Signal',linewidth= 2)
plt.plot(t, y2, label='Shifted by +3pi/2')
plt.plot(t, y3,label='Shifted by +pi/2')
plt.plot(t, y4,label='Shifted by +pi/3')
plt.plot(t, y5,label='Shifted by +pi/4')
plt.xlabel('Time (s)',fontsize=14)
plt.ylabel('Voltage (V)',fontsize=14)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title("Phase Difference",fontsize=18)
plt.show()
os.chdir(r'E:\Student\EAfiles\AR')
import phase_difference as pdiff
phase_diff1 = pdiff.hilbert_phase(y1,y2)
phase_diff2 = pdiff.hilbert_phase(y1,y3)
phase_diff3 = pdiff.hilbert_phase(y1,y4)
phase_diff4 = pdiff.hilbert_phase(y1,y5)
# plv1 = pdiff.phase_locking_value(y1,y2)
# plv2 = pdiff.phase_locking_value(y1,y3)
# plv3 = pdiff.phase_locking_value(y1,y4)
# plv4 = pdiff.phase_locking_value(y1,y5)
# print(plv1)
# print(plv2)
# print(plv3)
# print(plv4)


