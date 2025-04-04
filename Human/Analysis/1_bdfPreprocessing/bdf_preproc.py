"""
Created on Sat Sep 17 2022

@author: Andrew Sivaprakasam
EEG pre-processing from BDF to pool data/events across files 
based on a trial identifier

(based on Rav's EEGpp.py)
"""

import glob
from anlffr.helper import biosemi2mne as bs
from mne import concatenate_raws
import numpy as np

def poolBDF(trialID, refchans, filtband, fs_new = 8e3, exclude=[], add_subderms = False):
    files = glob.glob(trialID+'.bdf');
    files.sort();
    files.insert(0,files.pop(len(files)-1))
    print(files)
    raw = [];
    eves = [];
    
    for f in files:
        raw_temp, eves_temp = bs.importbdf(bdfname = f, refchans = refchans, exclude = exclude)
        #raw_temp, eves_temp = bs.importbdf(f)
        
        #Sometimes used simultaneously with chinchilla EEG Cap
        if add_subderms:
            raw_temp.set_channel_types({'EXG3':'eeg'});
            raw_temp.set_channel_types({'EXG4':'eeg'});
            raw_temp.set_channel_types({'EXG5':'eeg'});
            
        raw_temp.filter(l_freq = filtband[0], h_freq = filtband[1], picks = 'all')
        if fs_new:
            print('Resampling to ' + str(fs_new) + 'Hz and updating event indices')
            eves_temp[:,0] = np.round(eves_temp[:,0]/raw_temp.info['sfreq']*fs_new).astype('int')
            raw_temp.resample(fs_new)
            
        raw.append(raw_temp)
        eves.append(eves_temp)
        
        
    EEG_full, eves_full = concatenate_raws(raw,events_list=eves)
    
    return EEG_full, eves_full, files 

