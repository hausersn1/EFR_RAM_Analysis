
% Get all humans
cd(prefix)
all_humans = dir(sprintf('*/S*/Processed/S*.mat'));
numFiles = numel(all_humans); 

% initialize data matricies
subject=cell(numFiles, 1); 
group= cell(numFiles,1);
PLV = zeros(numFiles,2801);
f= zeros(numFiles,2801);
DFT = zeros(numFiles,2801);
LOCS = zeros(numFiles,4);
PLV_pks = zeros(numFiles,4);
DFT_pks = zeros(numFiles,4);
t = zeros(numFiles,5601);
t_env = zeros(numFiles,5601);

% Load data
for z = 1:numFiles
    subj = all_humans(z).name(1:4);
    cd(all_humans(z).folder)
    condition = extractBetween(all_humans(z).name, '_223_', '.mat'); 

    % Load the file
    load(all_humans(z).name, 'data')

    % Save data needed for plotting
    subject{z,1} = subj; 
    group{z,1} = condition;
    PLV(z,:) = data.result.PLV_env';
    f(z,:) = data.result.f;
    DFT(z,:) = data.result.DFT_corrected'; 
    LOCS(z,:) = data.result.PLV_LOCS; 
    PLV_pks(z,:) = data.result.PLV_PKS; 
    DFT_pks(z,:) = data.result.DFT_PKS; 
    t(z,:) = data.result.t; 
    t_env(z,:) = data.result.T_env'; 

end


%% Plot Data

blck = [0.25, 0.25, 0.25];
rd = [194 106 119]./255; %TTS
blu = [148 203 236]./255; %CA
yel = [220 205 125]./255; %PTS
gre = [93 168 153]./255; %GE

i_rd = [194 106 119 ]./255; %TTS
i_blu = [148 203 236 ]./255; %CA
i_yel = [220 205 125 ]./255; %PTS
i_gre = [93 168 153]./255; %GE

i_cols = [ i_rd; i_blu; i_yel; i_gre];
cols = [blck; rd; blu; yel; gre];
groups = {'YNH', 'MANH', 'HL'};
subp = [1 2 3]';

figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 12, 8])

for j = 1:numFiles
        grp = strcmp(groups, group{j,1});
        if sum(grp) > 0
            subplot(1,3,grp * subp)
            hold on;
            
            % Plot DFT
            plot(f(j,:)./1e3, DFT(j,:), '-','Color',i_cols(grp,:),'linewidth',1)
            alpha(.3)
            plot(LOCS(j,:)./1e3, DFT_pks(j,:), '*','Color',i_cols(grp,:),'linewidth',1)
            alpha(.3)

            %text(9,baseline(j,5), all_chins(j).name, 'Units', 'Data', 'Color', i_blck(1,1:3))
            %text(9,post(j,5), all_chins(j).name, 'Units', 'Data', 'Color',grp * i_cols)
        end
end

% YNH
subplot(1,3,1)
ylim([-0, .5])
xlim([.15,2])
xticks([.25, .5, 1, 2])
yticks([.02:0.02:1].*5)
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Young Normal Hearing', 'FontSize', 18, 'Color', rd)
set(gca, 'XScale', 'log', 'FontSize', 14)


% MANH
subplot(1,3,2)
ylim([-0, .1])
xlim([.15,2])
xticks([.25, .5, 1, 2])
yticks([.02:0.02:1])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Middle Aged NH', 'FontSize', 18, 'Color', blu)
set(gca, 'XScale', 'log', 'FontSize', 14)


% HL
subplot(1,3,3)
ylim([-0, .1])
xlim([.15,2])
xticks([.25, .5, 1, 2])
yticks([.02:0.02:1])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Hearing Loss', 'FontSize', 18, 'Color', yel)
set(gca, 'XScale', 'log', 'FontSize', 14)


%% Boxplot of Amplitudes

% Colors for humans
y_col = [046 037 133]./255; 
m_col = [159 074 150]./255; 
h_col = [051 117 056]./255; 

% Set up for boxplot
for k = 1:length(group)
    if strcmp(group{k,1}, 'YNH')
        gr_code(1,k)=0; 
    elseif strcmp(group{k,1}, 'MANH')
        gr_code(1,k)=1; 
    elseif strcmp(group{k,1}, 'HL')
        gr_code(1,k)=2; 
    end
end

%% actually plot
figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 3.3, 2.5])
boxplot(sum(DFT_pks,2)', gr_code, 'PlotStyle', 'Traditional', 'ColorGroup', gr_code, 'Colors',[y_col;m_col; h_col])
hold on; 
scatter(ones(size(gr_code(gr_code==0))).*(1+(rand(size(gr_code(gr_code==0)))-0.5)/10),sum(DFT_pks(gr_code==0,:),2), 10, 'filled', 'MarkerEdgeColor',y_col, 'MarkerFaceColor',y_col)
scatter(2.*ones(size(gr_code(gr_code==1))).*(1+(rand(size(gr_code(gr_code==1)))-0.5)/10),sum(DFT_pks(gr_code==1,:),2),10, 'filled', 'MarkerEdgeColor',m_col, 'MarkerFaceColor',m_col)
scatter(3.*ones(size(gr_code(gr_code==2))).*(1+(rand(size(gr_code(gr_code==2)))-0.5)/10),sum(DFT_pks(gr_code==2,:),2),10, 'filled', 'MarkerEdgeColor',h_col, 'MarkerFaceColor',h_col)
set(gca, 'XTickLabel', {'YNH', 'MANH', 'HL'})
xlabel('Group')
ylabel('\SigmaEFR_{1:4} Amplitude (\muV)')
ylim([0,.6])
set(gca, "FontSize", 11)
title('RAM EFR', 'FontWeight', 'bold', 'FontSize', 14)

%%
% Set up for boxplot
for k = 1:length(group)
    if strcmp(group{k,1}, 'YNH_tiptrodes')
        gr_code(1,k)=0; 
    elseif strcmp(group{k,1}, 'MANH_tiptrodes')
        gr_code(1,k)=1; 
    elseif strcmp(group{k,1}, 'HL_tiptrodes')
        gr_code(1,k)=2; 
    else
        gr_code(1,k) = nan;
    end
end

%% actually plot
figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 3.3, 2.5])
boxplot(sum(DFT_pks,2)', gr_code, 'PlotStyle', 'Traditional', 'ColorGroup', gr_code, 'Colors',[y_col;m_col; h_col])
hold on; 
scatter(ones(size(gr_code(gr_code==0))).*(1+(rand(size(gr_code(gr_code==0)))-0.5)/10),sum(DFT_pks(gr_code==0,:),2), 10, 'filled', 'MarkerEdgeColor',y_col, 'MarkerFaceColor',y_col)
scatter(2.*ones(size(gr_code(gr_code==1))).*(1+(rand(size(gr_code(gr_code==1)))-0.5)/10),sum(DFT_pks(gr_code==1,:),2),10, 'filled', 'MarkerEdgeColor',m_col, 'MarkerFaceColor',m_col)
scatter(3.*ones(size(gr_code(gr_code==2))).*(1+(rand(size(gr_code(gr_code==2)))-0.5)/10),sum(DFT_pks(gr_code==2,:),2),10, 'filled', 'MarkerEdgeColor',h_col, 'MarkerFaceColor',h_col)
set(gca, 'XTickLabel', {'YNH', 'MANH', 'HL'})
xlabel('Group')
ylabel('\SigmaEFR_{1:4} Amplitude (\muV)')
ylim([0,.6])
set(gca, "FontSize", 11)
title('RAM EFR', 'FontWeight', 'bold', 'FontSize', 14)

%%
%% For LME
Subject = subject; 
Group = group; 
EFR_mag = sum(DFT_pks(:,1:4),2);
RAM_amp_1 = [DFT_pks(:,1)]; 
RAM_amp_2 = [DFT_pks(:,2)]; 
RAM_amp_3 = [DFT_pks(:,3)]; 
RAM_amp_4 = [DFT_pks(:,4)]; 


ram_table = table(Subject, Group, EFR_mag, RAM_amp_1, ...
    RAM_amp_2, RAM_amp_3, RAM_amp_4); 

cd('C:\Users\saman\Desktop\Pitch_Diagnostics_Code\Profiling\Processed')
writetable(ram_table, 'Data_RAM_human.csv')
cd(cwd)
