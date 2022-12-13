# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 23:56:26 2022

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
group = int(input('Press 1 for analyzing saline group and 2 for epileptic group\n'))
smrfile = input ("Enter the SMR file: \n")
#Paste values from the excel file
filename= (EEG_path+smrfile)
#Paste values from the excel file 
yml = input("Enter the YML file : \n")
yml_channel = int(yml[-1:])
print(yml_channel)
if yml_channel == 1:
    ch_ipsi = 'HCi1'
    ch_contra = 'HCc'
elif yml_channel == 2:
    ch_ipsi = 'HCi1_2'
    ch_contra = 'HCc_2'
elif yml_channel == 2:
    ch_ipsi = 'HCi1_2'
    ch_contra = 'HCc_2'
elif yml_channel == 3:
    ch_ipsi = 'HCi1_3'
    ch_contra = 'HCc_3'
elif yml_channel == 4:
    ch_ipsi = 'HCi1_4'
    ch_contra = 'HCc_4'

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
raw_ipsi = (si[ch_ipsi]['trace'])
raw_contra = (si[ch_contra]['trace'])

sampling_freq = 10000
t = np.float32(np.arange(0.0, ((len(raw_ipsi))/sampling_freq), 1/sampling_freq))
t=np.around(t, decimals=4, out=None)
#%%
if group ==1:
    start = 1336*sampling_freq
    stop =  1950 *sampling_freq
    #%% Confirmation
    ii_ipsi = raw_ipsi[start:stop]
    iit = t[start:stop]
    ii_contra = raw_contra[start:stop]
    iit = t[start:stop]
    plt.plot(t,raw_ipsi,linewidth=0.5)
    plt.plot(iit,ii_ipsi,linewidth=0.5)
    plt.grid('on')

else :
    ymlfile = (yml_path+yml+'__runparams.yml')
    aRec = eam.Rec(init_ymlpath=ymlfile)
    # Compute States:
    SA = eana.StateAnalysis()
    SA.run(aRec)
    ii_phases = [S for S in aRec.states if S.state in ('IIP', 'depr')] #ii_phases=interictal phases
    starts=[S.start for S in ii_phases]
    stops=[S.stop for S in ii_phases]
    aRec.plot(['raw','bursts','states']) #CHECK, are you happy with the states?
    
    #create DataFrame
    import pandas as pd
    a_df = pd.DataFrame({'starts': starts,
                       'stops': stops})
    a = a_df.to_numpy()
    if len(a)==0:
        ii_phases = [S for S in aRec.states if S.state in ('IIP', 'depr','void')] #ii_phases=interictal phases
        starts=[S.start for S in ii_phases]
        stops=[S.stop for S in ii_phases]
        a_df = pd.DataFrame({'starts': starts,
                           'stops': stops})
        a = a_df.to_numpy()
    
    
    #%% Newer approach
    
    ii_ipsi = np.float32([])
    ii_contra = np.float32([])
    
    for i in range(0,len(a)):
            ii_ipsi = np.append(ii_ipsi,raw_ipsi[int(a[i,0]*sampling_freq):int(a[i,1]*sampling_freq)])
            ii_contra = np.append(ii_contra,raw_contra[int(a[i,0]*sampling_freq):int(a[i,1]*sampling_freq)])
    #%% Getting interictal timevector (iit)
    iit = np.float32(np.arange(0.0, ((len(ii_ipsi))/sampling_freq), 1/sampling_freq))
    iit = np.around(iit, decimals=4, out=None)
iit_len= len(iit)/sampling_freq
#%% Filtering 
os.chdir(r'E:\Student\EAfiles\AR')
import band_pass_filter as fil
filtering_sig = fil.filtering(ii_ipsi,ii_contra)
2#%% Filtered signal to analytic signal

filt_ii_ipsi = filtering_sig[0][:,]
filt_ii_contra = filtering_sig[1][:,]
analytic_ipsi = hilbert(filt_ii_ipsi)
analytic_contra = hilbert(filt_ii_contra)

#%% Getting phase angles

instantaneous_phase_ipsi = np.angle(analytic_ipsi)
instantaneous_phase_contra = np.angle(analytic_contra)

#%% Obtaining the PLV 

import phase_difference as pdiff
phase_diff = pdiff.hilbert_phase(filt_ii_ipsi,filt_ii_contra)
#plt.plot(iit_m1,phase_diff)
plv = pdiff.phase_locking_value(filt_ii_ipsi,filt_ii_contra)
print(plv)

#%% Plots (Takes more computational time, run when required)

# fig,axs = plt.subplots(2,1,facecolor='w',
#                        figsize=(16,7),sharey=True)
# # sharey creates common y axis markings for each channel 
# # for comparative visualization of two channels properly
# fig.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
# axs[0].plot(t,raw_ipsi,linewidth=0.5,color='black')
# axs[1].plot(t,raw_contra,linewidth=0.5,color='black')
# #%%
# f,ax = plt.subplots(2,1,facecolor='w',figsize=(16,3),sharex= True) # sharex zooms together both subplots 
# f.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
# ax[0].plot(iit,ii_ipsi,linewidth=0.5,color='black')
# ax[0].set_ylim(-5,5)
# ax[1].plot(iit,ii_contra,linewidth=0.5,color='black')
# ax[1].set_ylim(-5,5)
# #%%
# fig,axs = plt.subplots(2,1,facecolor='w',
#                        figsize=(16,7),sharey=True)
# # sharey creates common y axis markings for each channel 
# # for comparative visualization of two channels properly
# fig.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
# axs[0].plot(t,raw_ipsi,linewidth=0.5,color='black')
# axs[1].plot(t,raw_contra,linewidth=0.5,color='black')

# f,ax = plt.subplots(2,1,facecolor='w',figsize=(16,3),sharex= True) # sharex zooms together both subplots 
# f.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
# ax[0].plot(iit,ii_ipsi,linewidth=0.5,color='black')
# ax[0].set_ylim(-5,5)
# ax[1].plot(iit,ii_contra,linewidth=0.5,color='black')
# ax[1].set_ylim(-5,5)
# #%%
# # plt.figure(figsize=(16,7))
# # plt.plot(iit,filt_ii_ipsi)
# # plt.plot(iit,ii_ipsi,alpha=0.4)
# # plt.xlim(180,190)
# # plt.ylim(-2,2)

# f,ax = plt.subplots(2,1,facecolor='w',figsize=(16,3),sharex= True) # sharex zooms together both subplots 
# f.subplots_adjust(left=0.04,right=0.99,top=0.98,bottom=0.2)
# ax[0].plot(iit,ii_ipsi,linewidth=0.5,color='black',alpha =0.3)
# ax[0].plot(iit,filt_ii_ipsi,linewidth=1,color='magenta')
# ax[0].set_ylim(-2,2)
# ax[1].plot(iit,ii_contra,linewidth=0.5,color='black',alpha =0.3)
# ax[1].plot(iit,filt_ii_contra,linewidth=1,color='blue')
# ax[1].set_ylim(-2,2)
