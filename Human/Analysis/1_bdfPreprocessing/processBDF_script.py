#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 11:13:02 2023

@author: samhauser
Converting AS jupyter notebook to regular script.
"""

import os
import sys
import warnings

import mne
from anlffr.helper import biosemi2mne as bs
from matplotlib import pyplot as plt
import numpy as np
from bdf_preproc import poolBDF

from anlffr.preproc import find_blinks
from mne import compute_proj_epochs

import scipy.io

subj = 'S356';
cond = 'MANH';
fmod = 223;
EFR = 1; #change to 1 if you want to look at EFRs

isAndrew = 0;
local = 0;

if isAndrew:
    if local:
        #Local Storage
        measure_dir = '/mnt/ECF4E22CF4E1F92A/Research/Data_Backup/Pitch_Study/F30_Full_Data/ACC/SNAPLab/';

    else:
        #Ext Drive
        measure_dir = '/media/sivaprakasaman/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/EFR_RAM/Human/';
else:
    measure_dir = 'D:/THESIS/Pitch_Diagnostics_Data/EFR_RAM/Human/'
    # print('please define your own data_dir and out_loc');
    sys.path.append('C:/Users/saman/Desktop/Code/mne-python/')
    sys.path.append('C:/Users/saman/Desktop/Code/ANLffr/') 

## Start

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['figure.dpi']  = 120


pwd = os.getcwd();

trialID = '*RAM*' + str(fmod) +'*';

data_dir = measure_dir+cond+'/'+subj;
out_loc = data_dir + '/Preprocessed';

os.chdir(data_dir+'/Raw');

filtband = [60,4000];
fs_new = 2*filtband[1];

## Reference to Earlobes
refchans = ['EXG1','EXG2'];
raw, eves, files = poolBDF(trialID, refchans, filtband, fs_new = fs_new)
fname = subj+'_RAM_' + str(fmod) +'_EFR_preProcessed.mat';

## Reference to tiptrodes
refchans2 = ['EXG3','EXG4'];
raw2, eves2, files2 = poolBDF(trialID, refchans2, filtband, fs_new = fs_new)
fname2 = subj+'_RAM_' + str(fmod) +'_EFR_preProcessed_tiptrodes.mat';

os.chdir(pwd);

## 
#raw.plot(duration = 25.0, scalings = dict(eeg=80e-6), n_channels = 37);

bad_chans = []; #A29 bad for S366
raw.drop_channels(bad_chans);
raw.info

if not EFR:
    blinks = find_blinks(raw)
    epochs_blinks = mne.Epochs(raw, blinks, event_id=998, baseline=(-0.25, 0.25),
                               reject=dict(eeg=500e-6), tmin=-0.25, tmax=0.25)

    blink_proj = compute_proj_epochs(epochs_blinks, n_eeg=1)
    raw.add_proj(blink_proj)

tbounds = [-0.2,1.2];
bsline = (-.2,0);

pos = 1;
neg = 2;

to_proj = bool(~EFR);

epochs_p = mne.Epochs(raw, eves, event_id=pos, baseline=bsline, proj=to_proj,
                    tmin=tbounds[0], tmax=tbounds[1], reject=dict(eeg=200e-5), verbose = 'ERROR')

epochs_n = mne.Epochs(raw, eves, event_id=neg, baseline=bsline, proj=to_proj,
                    tmin=tbounds[0], tmax=tbounds[1], reject=dict(eeg=200e-5), verbose = 'ERROR');

chan_names = epochs_p.ch_names;
efr_chans = ['A32'];
cort_chans= efr_chans;

if EFR:
#     chans2pick = mne.pick_channels_regexp(chan_names,'A.');
    chans2pick = mne.pick_channels(chan_names,efr_chans);
else:
    chans2pick = mne.pick_channels(chan_names,cort_chans);

all_epochs_mean_cap = epochs_p.get_data();
all_epochs_mean_cap = all_epochs_mean_cap[:,chans2pick,:];
all_epochs_mean_p = np.mean(all_epochs_mean_cap,1);

cap_p = np.mean(all_epochs_mean_p,0);
cap_p_std_err = np.std(all_epochs_mean_p,0)/np.sqrt(np.size(all_epochs_mean_p,0));

all_epochs_mean_cap = epochs_n.get_data();
all_epochs_mean_cap = all_epochs_mean_cap[:,chans2pick,:];
all_epochs_mean_n = np.mean(all_epochs_mean_cap,1);

cap_n = np.mean(all_epochs_mean_n,0);
cap_n_std_err = np.std(all_epochs_mean_n,0)/np.sqrt(np.size(all_epochs_mean_n,0));

chan_id = np.array(chans2pick);
selected_chans = str([chan_names[index] for index in chan_id]);
tmin = tbounds[0];
t_vect = np.arange(0,np.size(cap_p,0))/fs_new;
t_vect = t_vect + tmin;

plt.rcParams['figure.figsize'] = [12, 8]
#plot params
xlims = tbounds;
# xlims = [0.1,0.8];
ylims = [-2e-6,1e-6];
# ylims = [-9e-7,9e-7];


buff = 0;
# plt.figure();
# plt.subplot(2,1,1);

# plt.plot(t_vect, cap_p-buff,linewidth=1,color = 'purple', label = "Mean of EEG chans "+ selected_chans)
# plt.fill_between(t_vect,cap_p+cap_p_std_err-buff,cap_p-cap_p_std_err-buff, color = 'purple', alpha=0.2)

# plt.xlim(xlims)
# plt.ylim(ylims)
# plt.legend(loc=3)
# plt.ylabel('Amplitude (V)')
# plt.title("+ polarity")

# plt.subplot(2,1,2);
# plt.plot(t_vect, cap_n-buff,linewidth=1,color = 'blue', label = "Mean of EEG chans "+ selected_chans)
# plt.fill_between(t_vect,cap_n+cap_n_std_err-buff,cap_n-cap_n_std_err-buff, color = 'blue', alpha=0.2)

# plt.xlim(xlims)
# plt.ylim(ylims)
# plt.legend(loc=3)
# plt.ylabel('Amplitude (V)')
# plt.title("- polarity")
# plt.rcParams.update({'font.size': 20})
# plt.show()

erp_p = epochs_p.average();
erp_n = epochs_n.average();

comb = mne.combine_evoked([erp_p,erp_n],[.5,.5]);
comb.plot(picks = chans2pick);

comb_export = comb.get_data(picks = chans2pick);

# topo_fig = plt.figure();
# topo_fig.clear();
# plt.rcParams.update({'figure.figsize': (4,4)})
# plt.rcParams.update({'lines.linewidth': 1})
# topo_fig = plt.figure(dpi = 300)
# ax = plt.gca();
# # topo_fig = erp_up.plot_topo(ylim = dict(eeg=[-4,4]),legend=False, axes = ax,title = 'ACCs', color = 'purple');
# topo_fig = comb.plot_topo(ylim = dict(eeg=[-2,2.]),legend=False, axes = ax, title = 'ACCs', color = 'black');
# # topo_fig.show()
# plt.show()



os.chdir(out_loc);
scipy.io.savemat(fname, {'time':t_vect,'fs':fs_new,'filt_band':filtband,
                         'mean_p_cap':cap_p, 'std_p_cap':cap_p_std_err,
                         'mean_n_cap':cap_n, 'std_n_cap':cap_n_std_err,
                         'cap_chan_ids':selected_chans,
                         'all_epochs_pos':all_epochs_mean_p,
                         'all_epochs_neg':all_epochs_mean_n,
                        'combined_mean':comb_export});
os.chdir(pwd)

## Again for ref2

if not EFR:
    blinks = find_blinks(raw2)
    epochs_blinks = mne.Epochs(raw2, blinks, event_id=998, baseline=(-0.25, 0.25),
                               reject=dict(eeg=500e-6), tmin=-0.25, tmax=0.25)

    blink_proj = compute_proj_epochs(epochs_blinks, n_eeg=1)
    raw.add_proj(blink_proj)


epochs_p = mne.Epochs(raw2, eves2, event_id=pos, baseline=bsline, proj=to_proj,
                    tmin=tbounds[0], tmax=tbounds[1], reject=dict(eeg=200e-5), verbose = 'ERROR')

epochs_n = mne.Epochs(raw2, eves2, event_id=neg, baseline=bsline, proj=to_proj,
                    tmin=tbounds[0], tmax=tbounds[1], reject=dict(eeg=200e-5), verbose = 'ERROR');


all_epochs_mean_cap = epochs_p.get_data();
all_epochs_mean_cap = all_epochs_mean_cap[:,chans2pick,:];
all_epochs_mean_p = np.mean(all_epochs_mean_cap,1);

cap_p = np.mean(all_epochs_mean_p,0);
cap_p_std_err = np.std(all_epochs_mean_p,0)/np.sqrt(np.size(all_epochs_mean_p,0));

all_epochs_mean_cap = epochs_n.get_data();
all_epochs_mean_cap = all_epochs_mean_cap[:,chans2pick,:];
all_epochs_mean_n = np.mean(all_epochs_mean_cap,1);

cap_n = np.mean(all_epochs_mean_n,0);
cap_n_std_err = np.std(all_epochs_mean_n,0)/np.sqrt(np.size(all_epochs_mean_n,0));


plt.rcParams['figure.figsize'] = [12, 8]
#plot params
xlims = tbounds;
# xlims = [0.1,0.8];
ylims = [-2e-6,1e-6];
# ylims = [-9e-7,9e-7];


buff = 0;
# plt.figure();
# plt.subplot(2,1,1);

# plt.plot(t_vect, cap_p-buff,linewidth=1,color = 'purple', label = "Mean of EEG chans "+ selected_chans)
# plt.fill_between(t_vect,cap_p+cap_p_std_err-buff,cap_p-cap_p_std_err-buff, color = 'purple', alpha=0.2)

# plt.xlim(xlims)
# plt.ylim(ylims)
# plt.legend(loc=3)
# plt.ylabel('Amplitude (V)')
# plt.title("+ polarity")

# plt.subplot(2,1,2);
# plt.plot(t_vect, cap_n-buff,linewidth=1,color = 'blue', label = "Mean of EEG chans "+ selected_chans)
# plt.fill_between(t_vect,cap_n+cap_n_std_err-buff,cap_n-cap_n_std_err-buff, color = 'blue', alpha=0.2)

# plt.xlim(xlims)
# plt.ylim(ylims)
# plt.legend(loc=3)
# plt.ylabel('Amplitude (V)')
# plt.title("- polarity")
# plt.rcParams.update({'font.size': 20})
# plt.show()

erp_p = epochs_p.average();
erp_n = epochs_n.average();

comb = mne.combine_evoked([erp_p,erp_n],[.5,.5]);
comb.plot(picks = chans2pick);

comb_export = comb.get_data(picks = chans2pick);

os.chdir(out_loc);
scipy.io.savemat(fname2, {'time':t_vect,'fs':fs_new,'filt_band':filtband,
                         'mean_p_cap':cap_p, 'std_p_cap':cap_p_std_err,
                         'mean_n_cap':cap_n, 'std_n_cap':cap_n_std_err,
                         'cap_chan_ids':selected_chans,
                         'all_epochs_pos':all_epochs_mean_p,
                         'all_epochs_neg':all_epochs_mean_n,
                        'combined_mean':comb_export});
