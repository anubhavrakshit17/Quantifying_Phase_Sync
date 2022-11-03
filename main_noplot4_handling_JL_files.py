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
raw_ipsi = (si['HCi1']['trace'])
raw_contra = (si['HCc']['trace'])
plt.plot(raw_ipsi)
#%%
ymlfile = (yml_path+yml+'__runparams.yml')
aRec = eam.Rec(init_ymlpath=ymlfile)
# Compute States:
SA = eana.StateAnalysis()
SA.run(aRec)
ii_phases = [S for S in aRec.states if S.state in ('IIP', 'void')] #ii_phases=interictal phases
starts=[S.start for S in ii_phases]
stops=[S.stop for S in ii_phases]
aRec.plot(['raw','bursts','states']) #CHECK, are you happy with the states?

#%%import pandas as pd

#create DataFrame
import pandas as pd
df = pd.DataFrame({'starts': starts,
                   'stops': stops})
a = df.to_numpy()

#%%
sampling_freq = 10000
t = np.float32(np.arange(0.0, ((len(raw_ipsi))/sampling_freq), 1/sampling_freq))
t=np.around(t, decimals=4, out=None)

#%%
# t_idx =[]
# for i in range(1,len(a)):
#    # start = np.asarray(np.where(t==a[i,0])).flatten()[0]
#     start = np.asarray(np.where(t==a[i,0])).flatten()   
#     t_idx = np.float32(np.append(t_idx,start)) 
#     stop = np.asarray(np.where(t==a[i,1])).flatten()
#     # stop =  np.asarray(np.where(t==a[i,1])).flatten()[0]
#     t_idx = np.int32(np.append(t_idx,stop))
    #%%
t_idx =[]
for i in range(1,len(a)):
    start = a[i,0]*sampling_freq
    t_idx = np.int32(np.append(t_idx,start))
    stop = a[i,1]* sampling_freq
    t_idx = np.int32(np.append(t_idx,stop))
#%%
#%%
#isolate interictal time-series
ii_ipsi = np.float32([])
for i in range(0,len(starts),2):
    ii_ipsi = np.append(ii_ipsi,raw_ipsi[t_idx[i]:t_idx[i+1]])
ii_contra = np.float32([])
for i in range(0,len(starts),2):
    ii_contra = np.append(ii_contra,raw_contra[t_idx[i]:t_idx[i+1]])
        
#%%
iit = np.float32(np.arange(0.0, ((len(ii_ipsi))/sampling_freq), 1/sampling_freq))
iit = np.around(iit, decimals=4, out=None)

#%%
filter = np.float32(signal.firwin(400, [0.0004, 0.0010], 
                                  pass_zero=False))

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

