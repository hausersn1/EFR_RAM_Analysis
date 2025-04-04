%Author: Andrew Sivaprakasam
%Updated: July, 2023
%Purpose: Script to import/plot/apply additional processing to RAM_EFR
%files (chin version)

fmod = 223;
file = 'p*RAM_223*.mat';
harmonics = 16;

fs = 8e3; %fs to resample to
t_win = [.2,.9]; %signal window, ignoring onset/offset effects
filts = [60,4000];
frames = round(t_win*fs);

%% Import data
if exist(datapath, 'dir')
    cd(datapath)
else
    fprintf('No such datapath - Subject: %s, Group: %s\n', subj, condition);
    return
end

datafile = dir(fullfile(cd,'*/p*.mat'));
if length(datafile) < 1
    fprintf('No file...Quitting!\n');
    return
end

cd(cwd);
numOfFiles = size(datafile,1);

for fileNum = 1:numOfFiles

    % set rng seed so that it's the same realization every time:
    load([cwd,filesep,'s.mat'])
    rng(s)

    cd(datafile(fileNum).folder)
    load(datafile(fileNum).name)
    cd(cwd)

    %% Date and time
    expdate = data.General.date;
    exptime = data.General.time;
    exptime(3:3:6) = '-';

    %% Data analysis & plotting:
    fs_orig = data.Stimuli.RPsamprate_Hz;

    all_dat = cell2mat(data.AD_Data.AD_All_V{1,1}');
    all_dat = all_dat';

    [b,a] = butter(4,filts./(fs_orig/2));
    all_dat = filtfilt(b,a,all_dat);

    all_dat = resample(all_dat,fs,round(fs_orig));
    all_dat = all_dat(frames(1):frames(2),:);

    pos = all_dat(:,1:2:end)*1e6/data.AD_Data.Gain; %+ polarity
    neg = all_dat(:,2:2:end)*1e6/data.AD_Data.Gain; %- polarity

    %% Get PLV spectra/Time domain waveform:

    %params for random sampling with replacement
    subset = 80;
    k_iters = 30;

    %only output things we want to look at
    [f, P1_env_raw, ~, PLV_env, ~, ~, T_env] = helper.getSpectAverage(pos,neg, fs, subset, k_iters);
    [floorx, floory] = helper.getNoiseFloor(pos,neg, fs, k_iters);
    t = (1:length(T_env))/fs;

    P1_env = P1_env_raw - floory';

    %% Get Peaks

    [PKS_PLV,LOCS_PLV] = helper.getPeaks(f,PLV_env,fmod,harmonics);
    [PKS_DFT,LOCS_DFT] = helper.getPeaks(f,P1_env,fmod,harmonics);

    % %% Plot: PLV
    % blck = [0.25, 0.25, 0.25];
    % rd = [0.8500, 0.3250, 0.0980, 0.5];
    % figure;
    %
    % %Spectral Domain
    % hold on;
    % title([subj,' | RAM - 25% Duty Cycle | F_{mod}',num2str(fmod),'|',condition],'FontSize',14);
    % plot(f,PLV_env,'Color',blck,'linewidth',1.5);
    % plot(LOCS,PKS,'*','Color',rd,'MarkerSize',10,'LineWidth',2);
    %
    % hold off;
    % ylim([0,1])
    % ylabel('PLV','FontWeight','bold')
    % xlabel('Frequency(Hz)','FontWeight','bold')
    %
    % %Time Domain
    % xstart = .6;
    % xend = .9;
    % ystart = 0.6;
    % yend = .9;
    %
    % axes('Position',[xstart ystart xend-xstart yend-ystart])
    % box on
    % hold on
    % plot(t, T_env,'Color',blck, 'LineWidth',2);
    % xlim([0.3,.4]);
    % ylim([-1,1]);
    % yticks([-1,0,1])
    % xlabel('Time(s)','FontWeight','bold');
    % ylabel('Amplitude \muV','FontWeight','bold')
    % hold off
    %
    % set(gcf,'Position',[50 100 560 420])

    %% Plot:
    if yes_plot
        blck = [0.25, 0.25, 0.25];
        rd = [0.8500, 0.3250, 0.0980, 0.5];
        figure;

        %Spectral Domain
        hold on;
        title([subj,' | RAM - 25% Duty Cycle | F_{mod}',num2str(fmod),'|',condition],'FontSize',14);
        plot(f,P1_env,'Color',blck,'linewidth',1.5);
        plot(LOCS_DFT,PKS_DFT,'*','Color',rd,'MarkerSize',10,'LineWidth',2);

        hold off;
        %ylim([0,10])
        ylabel('EFR Magnitude (\muV)','FontWeight','bold')
        xlabel('Frequency(Hz)','FontWeight','bold')

        %Time Domain
        xstart = .6;
        xend = .9;
        ystart = 0.6;
        yend = .9;

        axes('Position',[xstart ystart xend-xstart yend-ystart])
        box on
        hold on
        plot(t, T_env,'Color',blck, 'LineWidth',2);
        xlim([0.3,.4]);
        ylim([-1,1]);
        yticks([-1,0,1])
        xlabel('Time(s)','FontWeight','bold');
        ylabel('Amplitude \muV','FontWeight','bold')
        hold off

        set(gcf,'Position',[50 100 560 420])
    end
    %% Calculations based on the resulting data
    EFR_pk2nf = PKS_DFT;
    [EFR_pk2bl, ~] = helper.getPeaks(f,P1_env_raw,fmod,harmonics);

    clear data;
    %% Set a results structure
    resp.all_epochs.neg = neg;
    resp.all_epochs.pos = pos;
    info.subjgroup = condition;
    resp.numOfHarms = harmonics;
    resp.fmod = fmod;
    resp.filt_band = filts;
    resp.fs = fs;
    info.subj = subj;
    resp.time = t;
    result.t = t;
    result.T_env = T_env;
    result.f = f;
    result.PLV_env = PLV_env;
    result.PLV_PKS = PKS_PLV;
    result.PLV_LOCS = LOCS_PLV;
    result.DFT_PKS = PKS_DFT;
    result.DFT_LOCS = LOCS_DFT;
    result.DFT_raw = P1_env_raw;
    result.DFT_nf = floory';
    result.DFT_corrected = P1_env;
    result.sumDFTpk2nf = sum(EFR_pk2nf(1:4));

    data.info = info;
    data.resp = resp;
    data.result = result;

    %% Export:
    if export
        cd(datapath);
        cd ..
        cd('Processed')
        fname = [subj,'_RAM_',num2str(fmod), '_',condition, '_', expdate, '_', exptime];
        if yes_plot
            print(gcf,[fname,'_figure'],'-dpng','-r300');
        end
        save(fname,'data', 't','T_env','f','PLV_env','PKS_PLV','LOCS_PLV', 'P1_env', 'PKS_DFT', 'LOCS_DFT')
        cd(cwd);
        fprintf('Saved Data for subject %s\n', subj)
    end
end

cd(cwd);
