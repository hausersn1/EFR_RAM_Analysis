%Author: Andrew Sivaprakasam, Samantha Hauser
%Updated: January 2025
%Purpose: Script to import/plot/apply additional processing to RAM_EFR
%files
%Helpful Info: Be sure to define datapath so Import data section works.
%see my example setup_AS file.

fmod = 223;
harmonics = 8;

t_win = [.2,.9]; %signal window, ignoring onset/offset effects
filts = [60,4000];

%% Import data
if exist(datapath, 'dir')
    cd(datapath)
else
    fprintf('No such datapath - Subject: %s, Group: %s\n', subj, condition);
    return
end

datafile = dir(fullfile(cd,['*',num2str(fmod),'*preProcessed*.mat']));
if length(datafile) < 1
    fprintf('No file...Quitting!\n');
    return
end

numOfFiles = size(datafile,1); % also processes *_tiprodes.mat too

for fileNum = 1:numOfFiles

    % set rng seed so that it's the same realization every time:
    load([cwd,filesep,'s.mat'])
    rng(s)

    % change to the datapath and load the file (must be preprocessed, ie not BDF)
    cd(datapath)
    load(datafile(fileNum).name);
    file = datafile(fileNum).name;
    cd(cwd);

    %% Data analysis & plotting:
    fs = double(fs);
    frames_shift = round((t_win-min(time))*fs); %shift past baseline period in MNE

    pos = all_epochs_pos(:,frames_shift(1):frames_shift(2))'*1e6; %+ polarity
    neg = all_epochs_neg(:,frames_shift(1):frames_shift(2))'*1e6; %- polarity

    %% Get PLV spectra/Time domain waveform:

    %params for random sampling with replacement
    subset = 100;
    k_iters = 30;

    %only output things we want to look at
    [f, P1_env_raw, ~, PLV_env, ~, ~, T_env] = helper.getSpectAverage(pos,neg, fs, subset, k_iters);
    [floorx, floory] = helper.getNoiseFloor(pos,neg, fs, k_iters);
    t = (1:length(T_env))/fs;

    P1_env = P1_env_raw - floory';

    %% Get Peaks
    [PKS_PLV,LOCS_PLV] = helper.getPeaks(f,PLV_env,fmod,harmonics);
    [PKS_DFT,LOCS_DFT] = helper.getPeaks(f,P1_env,fmod,harmonics);

    %% Plot: PLV
    if yes_plot
        blck = [0.25, 0.25, 0.25];
        rd = [0.8500, 0.3250, 0.0980, 0.5];
        figure;

        %Spectral Domain
        hold on;
        title([subj,' | RAM - 25% Duty Cycle | F_{mod}',num2str(fmod),'|',condition],'FontSize',14);
        plot(f,PLV_env,'Color',blck,'linewidth',1.5);
        plot(LOCS_PLV,PKS_PLV,'*','Color',rd,'MarkerSize',10,'LineWidth',2);

        hold off;
        ylim([0,1])
        ylabel('PLV','FontWeight','bold')
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

        % Plot:
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
    resp.all_epochs.neg = all_epochs_neg;
    resp.all_epochs.pos = all_epochs_pos;
    resp.cap_chan_ids = cap_chan_ids;
    resp.combined_mean = combined_mean;
    info.subjgroup = condition;
    resp.numOfHarms = harmonics;
    resp.fmod = fmod;
    resp.filt_band = filt_band;
    resp.fs = fs;
    info.subj = subj;
    resp.time = time;
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
    result.sumDFTpk2nf_4 = sum(EFR_pk2nf(1:4));

    data.info = info;
    data.resp = resp;
    data.result = result;
    %% Export:
    if export
        cd(datapath);
        cd ..
        cd('Processed')
        if contains(file, 'tiptrodes')
            tiptrode = '_tiptrodes';
        else
            tiptrode = "";
        end
        fname = strcat(subj,'_RAM_',num2str(fmod), '_',condition, tiptrode);
        if yes_plot
            print(gcf,strcat(fname,'_figure'),'-dpng','-r300');
        end
        save(fname,'data', 't','T_env','f','PLV_env','PKS_PLV','LOCS_PLV', 'P1_env', 'PKS_DFT', 'LOCS_DFT')
        fprintf('Saved Data for subject %s\n', subj)
    end
    %% Add the data to the CSV

    if contains(file, 'tiptrodes')
        ref_chans = 'canals';
    else
        ref_chans = 'lobes';
    end

    cellData = {subj, condition, ref_chans, sum(data.result.DFT_PKS(1,1:4)), sum(data.result.DFT_PKS(1,:))...
        PKS_DFT(1), PKS_DFT(2), PKS_DFT(3), PKS_DFT(4), PKS_DFT(5), PKS_DFT(6), PKS_DFT(7), PKS_DFT(8), ...
        PKS_PLV(1), PKS_PLV(2), PKS_PLV(3), PKS_PLV(4), PKS_PLV(5), PKS_PLV(6), PKS_PLV(7), PKS_PLV(8) ...
        };
    M = [M; cellData];
    [C, ci] = unique(M);
    M = M(ci, :);

    cd(cwd);

end