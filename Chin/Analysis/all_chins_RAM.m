
% Get all humans
cd(prefix)
all_chins = dir(sprintf('Baseline/Q*/Processed/Q*_RAM_223_*_*.mat'));
numFiles = numel(all_chins); 

% initialize data matricies
chin_ID = cell(numel(all_chins),1); 
exposure = cell(numel(all_chins),1);

PLV_baseline = nan(numFiles,2801);
f_baseline= nan(numFiles,2801);
DFT_baseline = nan(numFiles,2801);
LOCS_baseline = nan(numFiles,16);
PLV_pks_baseline = nan(numFiles,16);
DFT_pks_baseline = nan(numFiles,16);
t_baseline = nan(numFiles,5601);
t_env_baseline = nan(numFiles,5601);

PLV_exp = nan(numFiles,2801);
f_exp= nan(numFiles,2801);
DFT_exp = nan(numFiles,2801);
LOCS_exp = nan(numFiles,16);
PLV_pks_exp = nan(numFiles,16);
DFT_pks_exp = nan(numFiles,16);
t_exp = nan(numFiles,5601);
t_env_exp = nan(numFiles,5601);


% Load data
for z = 1:numFiles
    subj = all_chins(z).name(1:4);
    cd(prefix)
    clear cond_cell
    chins_cond = dir(sprintf('*/%s/Processed/*-*', subj));
    
    if isempty(chins_cond)
        continue
    end
    
    for folder_i = 1:size(chins_cond,1)
        cond_cell(folder_i,1) = extractBetween(chins_cond(folder_i).folder,sprintf('%sChin%s',filesep,filesep),sprintf('%s%s%s',filesep, subj,filesep));
    end
    
    [conds,~,ind] = unique(cond_cell);
    cd(cwd);

    for c = 1:size(conds,1)
        condition = conds{c};
        suffix = [condition, filesep,subj,filesep,'Processed'];
        datapath = [prefix,suffix];
        cd(datapath)
        datafiles = dir(fullfile(cd,[subj, '_RAM_223_*_*.mat']));

       % Load the file
        if size(datafiles,1) > 1
            file = uigetfile('*_RAM_*.mat');
            load(file)
        elseif size(datafiles,1) == 1
            load(datafiles.name)
        else
            disp("File not found")
            break
        end

        if contains(condition,'Baseline')
            chin_ID{z,1} = subj;
            PLV_baseline(z,:) = data.result.PLV_env';
            f_baseline(z,:) = data.result.f;
            DFT_baseline(z,:) = data.result.DFT_corrected'; 
            LOCS_baseline(z,:) = data.result.PLV_LOCS; 
            PLV_pks_baseline(z,:) = data.result.PLV_PKS; 
            DFT_pks_baseline(z,:) = data.result.DFT_PKS; 
            t_baseline(z,:) = data.result.t; 
            t_env_baseline(z,:) = data.result.T_env'; 
        else
            exposure{z,1} = condition;
            PLV_exp(z,:) = data.result.PLV_env';
            f_exp(z,:) = data.result.f;
            DFT_exp(z,:) = data.result.DFT_corrected';
            LOCS_exp(z,:) = data.result.PLV_LOCS;
            PLV_pks_exp(z,:) = data.result.PLV_PKS;
            DFT_pks_exp(z,:) = data.result.DFT_PKS;
            t_exp(z,:) = data.result.t;
            t_env_exp(z,:) = data.result.T_env';
        end
    end
end
%% Plot

blck = [0.25, 0.25, 0.25];
rd = [194 106 119]./255; %TTS
blu = [148 203 236]./255; %CA
yel = [220 205 125]./255; %PTS
gre = [093 168 153]./255; %GE

i_rd = [194 106 119 ]./255; %TTS
i_blu = [148 203 236 ]./255; %CA
i_yel = [220 205 125 ]./255; %PTS
i_gre = [93 168 153]./255; %GE

i_cols = [ i_rd; i_blu; i_yel; i_gre];
cols = [blck; rd; blu; yel; gre];
groups = {'TTS', 'CA', 'PTS', 'GE'};
subp = [1 3 2 4]';

figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 12, 12])
hold on; 
for j = 1:numel(all_chins)
    if ~isempty(exposure{j,1})
        grp = strcmp(groups, extractBefore(exposure{j,1}, '_'));
        if sum(grp) > 0
            subplot(2,2,grp * subp)
            hold on;

            % Plot DFT
            plot(f_baseline(j,:)./1e3, DFT_baseline(j,:), '-','Color',blck(1,1:3),'linewidth',1.5)
            alpha(.3)
            plot(LOCS_baseline(j,1:8)./1e3, DFT_pks_baseline(j,1:8), '*','Color',blck(1,1:3),'linewidth',1.5, 'HandleVisibility','off')
            alpha(.3)
            % Plot DFT
            plot(f_exp(j,:)./1e3, DFT_exp(j,:), '-','Color',grp * i_cols,'linewidth',1.5)
            alpha(.3)
            plot(LOCS_exp(j,1:8)./1e3, DFT_pks_exp(j,1:8), '*','Color',grp * i_cols,'linewidth',1.5, 'HandleVisibility','off')
            alpha(.3)
        
        
        end
    end
end
legend(conds)

% TTS
subplot(2,2,1)
maxy = max(DFT_baseline,[],"all"); 
ylim([-0.05, maxy+.05])
xlim([.15,2])
xticks([.25:0.25:2])
yticks([-.05:0.05:.5])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Synaptopathy', 'FontSize', 18, 'Color', rd);
set(gca, 'XScale', 'linear', 'FontSize', 14)

% Carbo
subplot(2,2,3)
ylim([-0.05, maxy+.05])
xlim([.15,2])
xticks([.25:0.25:2])
yticks([-.05:0.05:.5])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('IHC Dysfunction', 'FontSize', 18, 'Color', blu);
set(gca, 'XScale', 'linear', 'FontSize', 14)

% PTS
subplot(2,2,2)
ylim([-0.05, maxy+.05])
xlim([.15,2])
xticks([.25:0.25:2])
yticks([-.05:0.05:.5])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Complex SNHL', 'FontSize', 18, 'Color', yel);
set(gca, 'XScale', 'linear', 'FontSize', 14)

% GE
subplot(2,2,4)
ylim([-0.05, maxy+.05])
xlim([.15,2])
xticks([.25:0.25:2])
yticks([-.05:0.05:.5])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('OHC Dysfunction', 'FontSize', 18, 'Color', gre);
set(gca, 'XScale', 'linear', 'FontSize', 14)




%% Boxplot of Amplitudes

xaxis = [0:1.5:5]
figure
hold on; 
for j = 1:numel(all_chins)
    if ~isempty(exposure{j,1})
        grp = strcmp(groups, extractBefore(exposure{j,1}, '_'));
        plot([grp*xaxis', grp*xaxis'+.5], [sum(DFT_pks_baseline(j,1:4),2), sum(DFT_pks_exp(j,1:4),2)], 'o-', 'Color', grp * i_cols, 'linew',4)
    end
end

xlim([-.5,5.5])
xticks([.25, 1.75, 3.25, 4.75])
xticklabels({'Syn','IHC','Complex', 'OHC'})
ylabel('Threshold (dB FPL)')
set(gca, 'FontSize', 20)
title('EFR Amp', 'FontSize', 24)
grid on;

%% Boxplot of PLV Ratio

xaxis = [0:1.5:5]

PLVratio_baseline = sum(PLV_pks_baseline(:, 3:8),2) ./ sum(PLV_pks_baseline(:, 1:2),2); 
PLVratio_exp = sum(PLV_pks_exp(:,3:8),2) ./ sum(PLV_pks_exp(:, 1:2),2); 
diff = PLVratio_exp - PLVratio_baseline; 
figure
hold on; 
for j = 1:numel(all_chins)
    if ~isempty(exposure{j,1})
        grp = strcmp(groups, extractBefore(exposure{j,1}, '_'));
        plot([grp*xaxis', grp*xaxis'+.5], [PLVratio_baseline(j,1), PLVratio_exp(j,1)], 'o-', 'Color', grp * i_cols, 'linew',4)
    end
end

figure
hold on; 
for j = 1:numel(all_chins)
    if ~isempty(exposure{j,1})
        grp = strcmp(groups, extractBefore(exposure{j,1}, '_'));
        plot([grp*xaxis'], [diff(j,1)], 'o-', 'Color', grp * i_cols, 'linew',4)
    end
end

figure
hold on; 
for j = 1:numel(all_chins)
    if ~isempty(exposure{j,1})
        grp = strcmp(groups, extractBefore(exposure{j,1}, '_'));
        plot(1:16, [PLV_pks_exp(j,:)], '-', 'Color', grp * i_cols, 'linew',2)
    end
end

xlim([-.5,5.5])
xticks([.25, 1.75, 3.25, 4.75])
xticklabels({'Syn','IHC','Complex', 'OHC'})
ylabel('Threshold (dB FPL)')
set(gca, 'FontSize', 20)
title('EFR Amp', 'FontSize', 24)
grid on;


%%
% Set up for boxplot
gr_code = []; 
for k = 1:length(exposure)
    if isempty(exposure{k,1})
        continue
    elseif contains(exposure{k,1}, 'TTS')
        gr_code(1,k)=1; 
    elseif contains(exposure{k,1}, 'CA')
        gr_code(1,k)=2; 
    elseif contains(exposure{k,1}, 'PTS')
        gr_code(1,k)=3; 
    elseif contains(exposure{k,1}, 'GE')
        gr_code(1,k)=4; 
    end
end
gr_code = [zeros(size(gr_code)), gr_code]; 
figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 6, 4])
boxplot([sum(DFT_pks_baseline(:,1:4),2)', sum(DFT_pks_exp(:,1:4),2)'], [gr_code], 'PlotStyle', 'Traditional', 'ColorGroup', gr_code, 'Colors',[blck; rd;blu;yel; gre])
hold on; 
scatter(ones(size(chin_ID)).*(1+(rand(size(chin_ID))-0.5)/10),sum(DFT_pks_baseline(:,1:4),2), 'filled', 'MarkerEdgeColor',blck, 'MarkerFaceColor',blck)
gr_code = gr_code(size(gr_code, 2)/2+1:end); 
scatter(2.*ones(size(gr_code(gr_code==1))).*(1+(rand(size(gr_code(gr_code==1)))-0.5)/10),sum(DFT_pks_exp(gr_code==1,1:4),2), 'filled', 'MarkerEdgeColor',rd, 'MarkerFaceColor',rd)
scatter(3.*ones(size(gr_code(gr_code==2))).*(1+(rand(size(gr_code(gr_code==2)))-0.5)/10),sum(DFT_pks_exp(gr_code==2,1:4),2), 'filled', 'MarkerEdgeColor',blu, 'MarkerFaceColor',blu)
scatter(4.*ones(size(gr_code(gr_code==3))).*(1+(rand(size(gr_code(gr_code==3)))-0.5)/10),sum(DFT_pks_exp(gr_code==3,1:4),2),'filled', 'MarkerEdgeColor',yel, 'MarkerFaceColor',yel)
scatter(5.*ones(size(gr_code(gr_code==4))).*(1+(rand(size(gr_code(gr_code==4)))-0.5)/10),sum(DFT_pks_exp(gr_code==4,1:4),2),'filled', 'MarkerEdgeColor',gre, 'MarkerFaceColor',gre)
set(gca, 'XTickLabel', {'Baseline', 'SYN', 'IHC', 'COMPLEX', 'OHC'})
xlabel('Group')
ylabel('\SigmaEFR_{1:4} Amplitude (\muV)')
ylim([0,1.25])
title('RAM EFR', 'FontWeight', 'bold', 'FontSize', 16)
set(gca, "FontSize", 14,'FontWeight','bold')


%%
% Set up for boxplot Difference
gr_code = []; 
for k = 1:length(exposure)
    if isempty(exposure{k,1})
        gr_code(1,k)=nan; 
    elseif contains(exposure{k,1}, 'TTS')
        gr_code(1,k)=0; 
    elseif contains(exposure{k,1}, 'CA')
        gr_code(1,k)=1; 
    elseif contains(exposure{k,1}, 'PTS')
        gr_code(1,k)=2; 
    elseif contains(exposure{k,1}, 'GE')
        gr_code(1,k)=3; 
    end
end

figure;

hold on;
plot([0,5], [0, 0], '--', 'Color', [.8, .8, .8], 'LineWidth',1)
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 4, 3])
boxplot(sum(DFT_pks_exp(:,1:4),2)'-sum(DFT_pks_baseline(:,1:4),2)', gr_code, 'PlotStyle', 'Traditional', 'ColorGroup', gr_code, 'Colors',[rd;blu;yel; gre])

hold on; 


for k = 1:length(exposure)
    if isempty(exposure{k,1})
        gr_code(1,k)=nan; 
    elseif contains(exposure{k,1}, 'TTS')
        gr_code(1,k)=0; 
    elseif contains(exposure{k,1}, 'CA')
        gr_code(1,k)=1; 
        %   if strcmp(chin_ID{k},'Q441')
        %     gr_code(1,k) = 5; 
        % end
    elseif contains(exposure{k,1}, 'PTS')
        gr_code(1,k)=2; 
        % if strcmp(chin_ID{k},'Q422') | strcmp(chin_ID{k},'Q443') | strcmp(chin_ID{k},'Q428')
        %     gr_code(1,k) = 4; 
        % end
    elseif contains(exposure{k,1}, 'GE')
        gr_code(1,k)=3; 
    end
end
scatter(ones(size(gr_code(gr_code==0))).*(1+(rand(size(gr_code(gr_code==0)))-0.5)/10), ...
    sum(DFT_pks_exp(gr_code==0,1:4),2)-sum(DFT_pks_baseline(gr_code==0,1:4),2), 20, ...
    'filled', 'MarkerEdgeColor',rd, 'MarkerFaceColor',rd)

scatter(2.*ones(size(gr_code(gr_code==1))).*(1+(rand(size(gr_code(gr_code==1)))-0.5)/10), ...
    sum(DFT_pks_exp(gr_code==1,1:4),2)-sum(DFT_pks_baseline(gr_code==1,1:4),2), 20,...
    'filled', 'MarkerEdgeColor',blu, 'MarkerFaceColor',blu)

scatter(3.*ones(size(gr_code(gr_code==2))).*(1+(rand(size(gr_code(gr_code==2)))-0.5)/10), ...
    sum(DFT_pks_exp(gr_code==2,1:4),2)-sum(DFT_pks_baseline(gr_code==2,1:4),2),20, ...
    'filled', 'MarkerEdgeColor',yel, 'MarkerFaceColor',yel)

scatter(4.*ones(size(gr_code(gr_code==3))).*(1+(rand(size(gr_code(gr_code==3)))-0.5)/10), ...
    sum(DFT_pks_exp(gr_code==3,1:4),2)-sum(DFT_pks_baseline(gr_code==3,1:4),2),20, ...
    'filled', 'MarkerEdgeColor',gre, 'MarkerFaceColor',gre)

% scatter(2.*ones(size(gr_code(gr_code==5))).*(1+(rand(size(gr_code(gr_code==5)))-0.5)/10), sum(DFT_pks_exp(gr_code==5,1:4),2)-sum(DFT_pks_baseline(gr_code==5,1:4),2), 15,  'filled', 'MarkerEdgeColor',blu, 'MarkerFaceColor',blu, 'Marker', 'x', 'HandleVisibility','off')
% scatter(3.*ones(size(gr_code(gr_code==4))).*(1+(rand(size(gr_code(gr_code==4)))-0.5)/10),sum(DFT_pks_exp(gr_code==4,1:4),2)-sum(DFT_pks_baseline(gr_code==4,1:4),2),15, 'filled', 'MarkerEdgeColor',yel, 'MarkerFaceColor',yel, 'Marker', 'x', 'HandleVisibility','off')


set(gca, 'XTickLabel', {'TTS', 'CA', 'PTS', 'GE'})
%xlabel('Exposure')
ylabel({'Change in Amplitude (\muV)';'\SigmaEFR_{1:4} Exp - Baseline'})
ylim([-.8, .3])
yticks([-.75:.25:.25])
set(gca, "FontSize", 12)
title('RAM EFR', 'FontWeight', 'bold', 'FontSize', 16)

%% Mean Plots


figure;
hold on;

CA_inds = strcmp(exposure,'CA_2wksPost'); 
ca_mean = [mean(DFT_pks_baseline(CA_inds,:),1,"omitmissing")',mean(DFT_pks_exp(CA_inds,:),1,"omitmissing")']; %col1 pre col2 post
ca_std= [std(DFT_pks_baseline(CA_inds,:),[],1,"omitmissing")',std(DFT_pks_exp(CA_inds,:),[],1,"omitmissing")'];

TTS_inds = strcmp(exposure,'TTS_2wksPost');
tts_mean = [mean(DFT_pks_baseline(TTS_inds,:),1,"omitmissing")',mean(DFT_pks_exp(TTS_inds,:),1,"omitmissing")']; 
tts_std= [std(DFT_pks_baseline(TTS_inds,:),[],1,"omitmissing")',std(DFT_pks_exp(TTS_inds,:),[],1,"omitmissing")'];

PTS_inds = strcmp(exposure,'PTS_2wksPost');
pts_mean = [mean(DFT_pks_baseline(PTS_inds,:),1,"omitmissing")',mean(DFT_pks_exp(PTS_inds,:),1,"omitmissing")']; 
pts_std= [std(DFT_pks_baseline(PTS_inds,:),[],1,"omitmissing")',std(DFT_pks_exp(PTS_inds,:),[],1,"omitmissing")'];


GE_inds = strcmp(exposure,'GE_1wkPost');
ge_mean = [mean(DFT_pks_baseline(GE_inds,:),1,"omitmissing")',mean(DFT_pks_exp(GE_inds,:),1,"omitmissing")']; 
ge_std= [std(DFT_pks_baseline(GE_inds,:),[],1,"omitmissing")',std(DFT_pks_exp(GE_inds,:),[],1,"omitmissing")'];

freq = LOCS/1000; 
%Plot

subplot(2,2,1);
hold on
errorbar(freq,tts_mean(:,1),tts_std(:,1),'o-','color',blck,'LineWidth',2.5)
errorbar(freq,tts_mean(:,2),tts_std(:,2),'o-','color',rd,'LineWidth',2.5)
hold off
xticks([0:.5:2.5]);
xlim([.0, 2])
yticks(0:.1:.5);
ylim([0,.5]);
xlabel('Frequency (kHz)');
ylabel('Amplitude (\muV)');
title('Synaptopathy','color',rd)
grid on


subplot(2,2,3);
hold on
errorbar(freq,ca_mean(:,1),ca_std(:,1),'o-','color',blck,'LineWidth',2.5)
errorbar(freq,ca_mean(:,2),ca_std(:,2),'o-','color',blu,'LineWidth',2.5)
hold off
xticks([0:.5:2.5]);
xlim([.0, 2])
yticks(0:.1:.5);
ylim([0,.5]);
xlabel('Frequency (kHz)');
ylabel('Amplitude (\muV)');
title('IHC Dys.','color',blu)
grid on

subplot(2,2,2);
hold on
errorbar(freq,pts_mean(:,1),pts_std(:,1),'o-','color',blck,'LineWidth',2.5)
errorbar(freq,pts_mean(:,2),pts_std(:,2),'o-','color',yel,'LineWidth',2.5)
hold off
xticks([0:.5:2.5]);
xlim([.0, 2])
yticks(0:.1:.5);
ylim([0,.5]);
xlabel('Frequency (kHz)');
ylabel('Amplitude (\muV)');
title('Complex SNHL','color',yel)
grid on

subplot(2,2,4);
hold on
errorbar(freq,ge_mean(:,1),ge_std(:,1),'o-','color',blck,'LineWidth',2.5)
errorbar(freq,ge_mean(:,2),ge_std(:,2),'o-','color',gre,'LineWidth',2.5)
hold off
xticks([0:.5:2.5]);
xlim([.0, 2])
yticks(0:.1:.5);
ylim([0,.5]);
xlabel('Frequency (kHz)');
ylabel('Amplitude (\muV)');
title('OHC Dys.','color',gre)
grid on

%set(gcf,'Position',[675 240 1012 725])

cd(cwd);

%% Save Table
Subject = [chin_ID; chin_ID];  
Group = [exposure; exposure];
stat = repmat([zeros(size(chin_ID)); 2*ones(size(chin_ID))], 1,1); 
RAM_baseline = sum(DFT_pks_baseline(:,1:4),2);
RAM_exp = sum(DFT_pks_exp(:,1:4),2); 
RAM_amplitude_sum = [RAM_baseline; RAM_exp]; 
RAM_amp_1 = [DFT_pks_baseline(:,1); DFT_pks_exp(:,1)]
RAM_amp_2 = [DFT_pks_baseline(:,2); DFT_pks_exp(:,2)]
RAM_amp_3 = [DFT_pks_baseline(:,3); DFT_pks_exp(:,3)]
RAM_amp_4 = [DFT_pks_baseline(:,4); DFT_pks_exp(:,4)]


ram_table = table(Subject, Group, stat, RAM_amplitude_sum, RAM_amp_1, ...
    RAM_amp_2, RAM_amp_3, RAM_amp_4); 

cd('D:\THESIS\Pitch_Diagnostics_Data\Processed')
writetable(ram_table, 'Data_RAM_2.csv')
cd(cwd)

subj = repmat(Subject, 2,1); 
stat = [zeros(size(Subject)); ones(size(Subject))]; 
group = repmat(Group, 2, 1); 
RAM = [RAM_baseline; RAM_exp]; 

ram_table = table(subj, group, stat, RAM); 

cd('C:\Users\saman\Desktop\Pitch_Diagnostics_Code\Profiling\Processed')
writetable(ram_table, 'Data_RAM_lme.csv')
cd(cwd)