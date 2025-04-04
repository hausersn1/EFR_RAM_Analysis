%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

uname = 'sivaprakasaman';

condition = 'YNH';
subj = 'S363';

prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/EFR_RAM/Human/'];
suffix = [condition,'/',subj,'/Preprocessed'];
datapath = [prefix,suffix];


processSubject;