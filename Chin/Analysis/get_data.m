function [t, T_env, PKS, PLV_env, f, LOCS] = get_data(prefix, subj, condition)

suffix = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
    filesep, 'Chin', filesep, condition, filesep, subj, filesep, 'Processed'];
datapath = [prefix,suffix];

% Import Data
cwd = pwd;
cd(datapath)
datafile = dir(fullfile(cd,[ subj, '_RAM_EFR_223_' condition, '*.mat']));
if length(datafile) < 1
    fprintf(sprintf('No file for %s %s...Filling with NaNs!\n', subj, condition));
    t=NaN(1,5601); 
    T_env=NaN(1,5601); 
    f=NaN(1,2801); 
    LOCS=NaN(1,16); 
    PKS=NaN(1,16); 
    PLV_env=NaN(1,2801);
elseif size(datafile,1) > 1
    checkDIR =uigetfile('.mat');
    load(checkDIR);
    file = checkDIR;
    load(file);
else
    load(datafile(1).name);
    file = datafile(1).name;
    load(file);
end



cd(cwd);