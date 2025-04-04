%% Run RAM Analysis (after preprocessing bdf file)
% Chin version
clear;
set(0,'defaultfigurerenderer','opengl')

%% Enter information here:
subj = 'Q431';                    % e.g., 'Q440'
condition = 'Baseline';           % e.g, 'HL', 'YNH', 'MANH'
location = 3;                       % 0 == mac, 1 == Desktop, 2 == SNAPlab, 3 == Laptop

export = 1;     % Save data?
yes_plot = 1;   % plot a figure?
mode = 0;       % 0 - Process single chins,
                % 1 - Process pre and post for single chin,
                % 2 - batch process all chins,
                % 3 - Only Plot Processed Group Data
                % 4 - Only Plot Pre-Post

%% Set directories
if location == 1
    comp = 'F:/';
elseif location == 0
    comp = '/Volumes/SNH/';
elseif location == 3
    comp = 'D:/';
end

% Set data directory: 
directories = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
    filesep, 'Chin', filesep];

prefix = [comp, directories];
cwd = pwd;

%% Run analysis

switch mode
    case 0
        disp("Processing Single Chin...");
        suffix = [ condition,filesep,subj, filesep, 'Raw'];
        datapath = [prefix,suffix];
        processChin;

    case 1
        disp("Processing Pre and Post for a Single Chin");
        cd(prefix)
        chins_cond = dir(sprintf('*/%s/Raw/*-*', subj));
        for folder_i = 1:size(chins_cond,1)
            cond_cell(folder_i,1) = extractBetween(chins_cond(folder_i).folder,'EFR_RAM\Chin\',sprintf('%s%s%s%s', filesep, subj, filesep, 'Raw'));
        end
        [conds,~,ind] = unique(cond_cell);
        cd(cwd);

        for c = 1:size(conds,1)
            condition = conds{c};
            suffix = [condition, filesep,subj,filesep,'Raw'];
            datapath = [prefix,suffix];
            processChin;
        end

        preVpost_RAM;
        cd(cwd)

    case 2
        disp("Batch Processing every chin, pre and post...This might take a while!")
        cd(prefix)
        all_chins = dir(sprintf('Baseline/Q*'));

        for z = 1:numel(all_chins)
            subj = all_chins(z).name;
            cd(prefix)
            chins_cond = dir(sprintf('*/%s/Raw/*-*', subj));
            if isempty(chins_cond)
                continue
            end
            clear cond_cell;
            for folder_i = 1:size(chins_cond,1)
                cond_cell(folder_i,1) = extractBetween(chins_cond(folder_i).folder,...
                    'EFR_RAM\Chin\',sprintf('%s%s%s%s', filesep, subj, filesep, 'Raw'));
            end

            [conds,~,ind] = unique(cond_cell);
            cd(cwd);

            for c = 1:size(conds,1)
                condition = conds{c};
                suffix = [condition, filesep,subj,filesep,'Raw'];
                datapath = [prefix,suffix];
                processChin;
            end
        end



    case 3
        disp("Plotting Group Data...")
        all_chins_RAM; %Need to fix
        cd(cwd)

    case 4
        disp("Plotting Pre and Post for a Single Chin");
        cd(prefix)
        chins_cond = dir(sprintf('*/%s/Processed/*-*', subj));
        if isempty(chins_cond)
            disp("No processed data for this chin yet...try mode 2")
            cd(cwd)
            return
        end

        for folder_i = 1:size(chins_cond,1)
            cond_cell(folder_i,1) = extractBetween(chins_cond(folder_i).folder,'EFR_RAM\Chin\',sprintf('%s%s%s%s', filesep, subj, filesep, 'Processed'));
        end
        [conds,~,ind] = unique(cond_cell);
        cd(cwd);
        preVpost_RAM;
        cd(cwd)
end