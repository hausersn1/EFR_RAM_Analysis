

for z = 1:length(conds)
    condition = conds{z};
    suffix = [condition, filesep,subj,filesep,'Processed'];
    datapath = [prefix,suffix];
    cd(datapath)
    
    file = dir("*_RAM_223*_*.mat"); 
    % Import Data
    % Load the file
    if ~isempty(file)
        load(file.name, 'data')
    else
        return
    end

    % Save data needed for plotting
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

%% Plot

blck = [0.25, 0.25, 0.25];
rd = [194 106 119]./255; %TTS
blu = [148 203 236]./255; %CA

i_rd = [194 106 119 ]./255; %TTS
i_blu = [148 203 236 ]./255; %CA
i_yel = [220 205 125 ]./255; %PTS
i_gre = [93 168 153]./255; %GE

cols = [blck; rd; blu];
groups = conds;
subp = [1 2]';

figure;
hold on;
set(gcf, 'Units', 'inches', 'Position', [.25, .25, 12, 8])
subplot(1,3,1:2)
hold on; 
for j = 1:length(conds)
    % Plot DFT
    plot(f(j,:)./1e3, DFT(j,:), '-','Color',cols(j,:),'linewidth',1.5)
    alpha(.3)
    plot(LOCS(j,1:8)./1e3, DFT_pks(j,1:8), '*','Color',cols(j,:),'linewidth',1.5, 'HandleVisibility','off')
    alpha(.3)
end
legend(conds)

% YNH
subplot(1,3,1:2)
maxy = max(DFT,[],"all"); 
ylim([-0.05, maxy+.05])
xlim([.15,2])
xticks([.25:0.25:2])
yticks([-.05:0.05:.5])
ylabel('Amplitude (\muV)','FontWeight','bold')
xlabel('Frequency(kHz)','FontWeight','bold')
title('Pre vs Post DFT', 'FontSize', 18)
set(gca, 'XScale', 'linear', 'FontSize', 14)


%% Boxplot of Amplitudes

% Set up for boxplot
for k = 1:length(group)
    if contains(group{k,1}, 'Baseline')
        gr_code(1,k)=0; 
    elseif contains(group{k,1}, 'Post')
        gr_code(1,k)=1; 
    end
end

%% actually plot
subplot(1,3,3);
hold on;
X = categorical(groups); 
numOfBars = size(X,1)
colors = [blck; rd]; 
for y = 1:numOfBars
    hold on;
    bar(X(y,1), sum(DFT_pks(y,1:4),2), 'FaceColor', colors(y,:))
end

set(gca, 'XTickLabel', {'Pre', 'Post'})
xlabel('Group')
ylabel('\SigmaEFR_{1:4} Amplitude (\muV)')

ylim([0,0.05+max(sum(DFT_pks(1,1:4),2), sum(DFT_pks(1,1:4),2))])
set(gca, "FontSize", 14,'FontWeight','bold')

cd(cwd)