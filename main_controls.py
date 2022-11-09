# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 19:59:03 2022

@author: Anubhav
"""
import os
import numpy as np
import array
from scipy.signal import firwin, lfilter, filtfilt, hilbert
from scipy.io import loadmat
from pylab import *
from scipy import signal
from colorama import Fore
main_path = 'E:\\Student\\' # CHANGE REQUIRED FOR NEW SETUP IN main_path
EAfiles_path = main_path+'\EAfiles\\' 
EAdetection_path = EAfiles_path+'\\EAdetection'
yml_path = EAfiles_path+'\YML\\yml_setup\\'
EEG_path = main_path+'\EEG\\'
os.chdir(EAdetection_path)
from core import ea_analysis as eana
from core import ea_management as eam
from core import helpers as hp

##%% Loading the raw smr files using neo module
smrfile = input ("Enter the SMR file: \n")
#Paste values from the excel file
filename= (EEG_path+smrfile)
#Paste values from the excel file 
yml = input("Enter the YML file : \n")
##%%
# # create a reader
# reader = neo.io.Spike2IO(filename) # Used for Epileptic PK85-86-87_38d000.smr
# #reader = neo.io.Spike2IO(filename,try_signal_grouping=False)
# # read the block
# data = reader.read(lazy=False)[0]
# seg = reader.read_segment(time_slice=None)
# raw_all_channels= np.squeeze(np.asarray(seg.analogsignals))
# print(shape(raw_all_channels))
#%%
si= hp.extract_smrViaNeo(filename)
#%%
# Change the channels
raw_ipsi = (si['HCi1_3']['trace'])
raw_contra = (si['HCc_3']['trace'])
# plt.plot(raw_ipsi)
#%%
sampling_freq = 10000
t = np.float32(np.arange(0.0, ((len(raw_ipsi))/sampling_freq), 1/sampling_freq))
t=np.around(t, decimals=4, out=None)
# plt.plot(t,raw_ipsi,linewidth=0.5)
# plt.grid('on')
#%%
# Define t index after seeing the signal portion having less spikes
# Manually see the control plot and put the data.
# start_input = int(input("Start Index from the graph: \n"))
# stop_input = int(input("Stop Index from the graph: \n"))
# start = start_input*sampling_freq
# stop =  stop_input *sampling_freq
start = 1336*sampling_freq
stop =  1950 *sampling_freq
#%% Confirmation
ii_ipsi = raw_ipsi[start:stop]
iit = t[start:stop]
plt.plot(t,raw_ipsi,linewidth=0.5)
plt.plot(iit,ii_ipsi,linewidth=0.5)
plt.grid('on')
#%%
# plt.figure()
ii_contra = raw_contra[start:stop]
iit = t[start:stop]
# plt.plot(t,raw_contra,linewidth=0.5)
# plt.grid('on')
# plt.plot(iit,ii_contra,linewidth=0.5)
# plt.grid('on')

#%%
filter = np.float32(signal.firwin(400, [0.0004, 0.0010], 
                                  pass_zero=False))
# filter = np.float32(signal.firwin(400, [0.0030, 0.010], 
#                                   pass_zero=False))

# First we convolve then we convert to float32 because 
#float64 hangs the prog.

filt_ii_ipsi = np.float32(signal.convolve(ii_ipsi, filter, 
                                       mode='same'))
filt_ii_contra = np.float32(signal.convolve(ii_contra, 
                                         filter, mode='same'))
#%%
analytic_ipsi = hilbert(filt_ii_ipsi)
analytic_contra = hilbert(filt_ii_contra)
instantaneous_phase_ipsi = np.angle(analytic_ipsi)
instantaneous_phase_contra = np.angle(analytic_contra)
#%%
os.chdir(r'E:\Student\EAfiles\AR')
import phase_difference as pdiff
phase_diff = pdiff.hilbert_phase(filt_ii_ipsi,filt_ii_contra)
#plt.plot(iit_m1,phase_diff)
plv = pdiff.phase_locking_value(filt_ii_ipsi,filt_ii_contra)
print(plv)
