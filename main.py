# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 15:20:28 2022

@author: Anubhav
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 20:00:27 2022

@author: Anubhav
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy
import sys
import os
import neo
import numpy as np
import array
from scipy.signal import firwin, lfilter, filtfilt, hilbert
from scipy.io import loadmat
from pylab import *
from scipy import signal
from colorama import Fore
EAdetection_path = 'E:\Student\EAfiles\EAdetection'
os.chdir(EAdetection_path)
from core import ea_analysis as eana
from core import ea_management as eam

#%% Loading the raw smr files using neo module
filename= "E:\Student\EEG\PK85-86-87_38d000.smr"

# create a reader
reader = neo.io.Spike2IO(filename)
# read the block
data = reader.read(lazy=False)[0]
seg = reader.read_segment(time_slice=None)
raw_all_channels= np.squeeze(np.asarray(seg.analogsignals))
 
#%% 
#The raw_all_channels contains all the channels from various mouse, 
# with each mouse having two channels each (ipsi and contra)
# so if there are 8 channels, then there are 4 mouse
# Or there are 3 mouse and one trigger input 

while True:
    try:
        mouse_select= int(input(Fore.RED+"Select which Mouse you want to select: "))
    except ValueError:
        print(Fore.RED+"Sorry, I didn't understand that.")
        continue
    if mouse_select not in range(1,(int(shape(raw_all_channels)[1]/2))+1):
        # shape(raw_all_channels) = (36000578,8)
        # shape(raw_all_channels)[1] = 8
        # shape(raw_all_channels)[1]/2 = 4.0
        # int(shape(raw_all_channels)[1]/2 = 4
        # range(1,4) = 1,2,3
        # range(1,5) = 1,2,3,4
        # range(1,4+1)= 1,2,3,4
        # range(1,(int(shape(raw_all_channels)[1]/2))+1) = 1,2,3,4
        
        
        print("Not an appropriate choice.")
    else:
        #age was successfully parsed!
        #we're ready to exit the loop.
            break
if mouse_select ==1:
    raw = raw_all_channels[:,0:2]
    
elif mouse_select == 2:
    raw = raw_all_channels[:,2:4]
    
elif mouse_select == 3:
    raw= raw_all_channels[:,4:6]

elif mouse_select == 4:
    raw = raw_all_channels[:,6:8]

#%% Make appropriate timevector

sampling_freq= 10000
t = np.float32(np.arange(0.0, ((len(raw))/sampling_freq), 1/sampling_freq))
t=np.around(t, decimals=4, out=None)
#%% Run Only when required
#Raw contains two channels, 1st is from ipsilateral electrode, 2nd is from contralateral electrode
raw_ipsi = raw[:,0]
raw_contra = raw[:,1]
#%% Plot the figure to see if channels from smr files are correctly imported

# How to subplot 
# fig, axs = plt.subplots(number of rows, no. of columns)
# fig.suptitle('Vertically stacked subplots')
# axs[0].plot(x, y)
# axs[1].plot(x, -y)

fig,axs = plt.subplots(2,1,facecolor='w',
                       figsize=(16,7),sharey=True)
# sharey creates common y axis markings for each channel 
# for comparative visualization of two channels properly
fig.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
axs[0].plot(t,raw_ipsi,linewidth=0.5,color='black')
axs[1].plot(t,raw_contra,linewidth=0.5,color='black')

#%% Getting the timepoints from yml and h5 code



ymlfile = 'E:\Student\EAfiles\YML\yml_setup\PK86_38d1_HCi1_2__runparams.yml'

# Load recording object
#ymlfile = ymlfiles[3]
aRec = eam.Rec(init_ymlpath=ymlfile)
#%%
# Compute States:
SA = eana.StateAnalysis()
SA.run(aRec)

# now you can plot also states!:
aRec.plot(['raw','bursts','states']) #CHECK, are you happy with the states?
#%%

ii_phases = [S for S in aRec.states if S.state in ('IIP', 'void')] #ii_phases=interictal phases
starts=[S.start for S in ii_phases]
stops=[S.stop for S in ii_phases]

start= np.asarray(np.where(t == starts[1])).flatten()
stop= np.asarray(np.where(t == stops[1])).flatten()
del starts, stops
ss= np.append(start,stop) # Put start stop to one single array
#isolate interictal time-series

ii_raw_ipsi=raw_ipsi[ss[0]:ss[1]]
ii_raw_contra = raw_contra[ss[0]:ss[1]]
iit= t[ss[0]:ss[1]]

#%% Time domain - Raw interictal

f,ax = plt.subplots(2,1,facecolor='w',figsize=(16,3),sharex= True) # sharex zooms together both subplots 
f.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
ax[0].plot(iit,ii_raw_ipsi,linewidth=0.5,color='black')
ax[0].set_ylim(-5,5)
ax[1].plot(iit,ii_raw_contra,linewidth=0.5,color='black')
ax[1].set_ylim(-5,5)
#%% frequency domain - Raw interictal

x= ii_raw_ipsi
x = x - x.mean()          # Subtract the mean from the data.
N = len(x)                 # Define the number of points in the trial.
xh = hanning(N) * x       # Multiply data by Hanning window,
xf = rfft(xh)[:-1]        # ... compute Fourier transform,

dt = t[1] - t[0]           # Define the sampling interval.
T = len(iit) / sampling_freq # Define the duration of the trial.
df = 1 / T;               # Determine frequency resolution.
fNQ = 1 / dt / 2;         # Determine Nyquist frequency.
Sxx = 2 * dt ** 2 / T * (xf * xf.conj());  #... and compute spectrum.

faxis = arange(0, fNQ, df) # Construct frequency axis
faxis= np.delete(faxis, 0) #dropping the first element, because of length mismatch

plt.figure()
semilogx(faxis, 10 * log10(Sxx.real))   # Plot decibels vs frequency,
xlim([df, 100])           # ... in limited frequency range,
xlabel('Frequency [Hz]')  # ... with axes labeled.
ylabel('Power [dB]')
title('Single-trial spectrum')
grid('on')

#%% For MATLAB Conversion [Optional]
# raw_ipsi_m = transpose(np.reshape(raw_ipsi, (len(raw),1)))
# raw_contra_m = transpose(np.reshape(raw_contra, (len(raw),1)))
# ii_raw_ipsi_m = transpose(np.reshape(ii_raw_ipsi, (len(ii_raw_ipsi),1)))
# ii_raw_contra_m = transpose(np.reshape(ii_raw_contra, (len(ii_raw_contra),1)))
# os.chdir(r'E:\Student\EAfiles\ar_matlab_data')

# matlab_save_select= int(input(Fore.RED+"Select 1 for raw and 2 for interictal export : "))
# if matlab_save_select ==1:
#     scipy.io.savemat('raw_ipsi_contra.mat', {'raw_ipsi': raw_ipsi,'raw_contra':raw_contra})
# elif matlab_save_select == 2:
#     scipy.io.savemat('ii_raw_ipsi_contra.mat', {'ii_raw_ipsi': ii_raw_ipsi,'ii_raw_contra':ii_raw_contra})
#%%

del raw_all_channels

#%% Filtering (Not used in python but in MATLAB, subject of investigation)
# Next step is Filtering.  Filtering is convolution of a window and the wave. 
# Visualize the window of the filter

filter = np.float32(signal.firwin(400, [0.0004, 0.0010], 
                                  pass_zero=False))

# First we convolve then we convert to float32 because 
#float64 hangs the prog.

ii_filt_ipsi = np.float32(signal.convolve(ii_raw_ipsi, filter, 
                                       mode='same'))
ii_filt_contra = np.float32(signal.convolve(ii_raw_contra, 
                                         filter, mode='same'))
#%% Time domain - Filtered interictal
plt.figure()
plt.plot(iit,ii_filt_ipsi,linewidth=0.5)

#%% Frequency Domain - Filtered interictal

fil_x= ii_filt_ipsi
fil_x = fil_x - fil_x.mean()          # Subtract the mean from the data.
fil_N = len(fil_x)                 # Define the number of points in the trial.
fil_xh = hanning(fil_N) * fil_x       # Multiply data by Hanning window,
fil_xf = rfft(fil_xh)[:-1]        # ... compute Fourier transform,

fil_dt = t[1] - t[0]           # Define the sampling interval.
fil_T = len(iit) / sampling_freq # Define the duration of the trial.
fil_df = 1 / fil_T;               # Determine frequency resolution.
fil_fNQ = 1 / fil_dt / 2;         # Determine Nyquist frequency.
fil_Sxx = 2 * fil_dt ** 2 / fil_T * (fil_xf * fil_xf.conj());  #... and compute spectrum.

fil_faxis = arange(0, fil_fNQ, fil_df) # Construct frequency axis
fil_faxis= np.delete(fil_faxis, 0) #dropping the first element, because of length mismatch
# Visualizing in frequency domain

plt.figure()
semilogx(fil_faxis, 10 * log10(fil_Sxx.real))   # Plot decibels vs frequency,
xlim([fil_df, 100])           # ... in limited frequency range,
xlabel('Frequency [Hz]')  # ... with axes labeled.
ylabel('Power [dB]')
title('Single-trial spectrum')
grid('on')

#%% Comparing both filtered and raw ipsilateral segments
ff,axf = plt.subplots(2,1)
axf[0].semilogx(faxis, 10 * log10(Sxx.real))
axf[1].semilogx(fil_faxis, 10 * log10(fil_Sxx.real))

#%% Hilbert Transform to extract phase information
analytic_ipsi = hilbert(ii_filt_ipsi)
analytic_contra = hilbert(ii_filt_contra)
#%% Extracting phase angles
instantaneous_phase_ipsi = np.angle(analytic_ipsi)
instantaneous_phase_contra = np.angle(analytic_contra)

#%%
plt.plot(iit, instantaneous_phase_ipsi)
plt.plot(iit, instantaneous_phase_contra)
plt.xlim(180,190)
#%% Phase difference
# Phase diff = angle (hilbert_ipsi - hilbert_contra)
# NOT abs(angle(hilbert_ipsi) - angle(hilbert_contra))
instantaneous_phase_diff_ipsi_contra = np.angle(analytic_contra-analytic_ipsi)
plt.plot(iit,instantaneous_phase_diff_ipsi_contra)
plt.xlim(180,190)

#%%
