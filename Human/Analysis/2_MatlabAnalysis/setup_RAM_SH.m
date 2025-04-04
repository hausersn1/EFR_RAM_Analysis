%% Run RAM Analysis (after preprocessing bdf file)
% Human version
clear;
set(0,'defaultfigurerenderer','opengl')

%% Enter information here:
subj = 'S115';                    % e.g., 'Q440'
condition = 'MANH';                 % e.g, 'HL', 'YNH', 'MANH'
uname = 'samhauser';
location = 3;                       % 0 == mac, 1 == Desktop, 2 == SNAPlab, 3 == Laptop

export = 1;     % Save data?
yes_plot = 1;   % plot a figure? 
mode = 0;       % 0 - Process single human,
                % 2 - batch process all humans,
                % 3 - Plot Processed Group Data

%% Set directories
if location == 1
    comp = 'F:/';
elseif location == 0
    comp = '/Volumes/SNH/';
elseif location == 3
    comp = 'D:/';
end

directories = ['THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'EFR_RAM', ...
    filesep, 'Human', filesep];

prefix = [comp, directories];
cwd = pwd;

%% Run analysis

% Load data table for saving
tableDir = [comp, 'THESIS', filesep, 'Pitch_Diagnostics_Data', filesep, 'compiledData', filesep, 'Human']; 
cd(tableDir)
M = readtable("Data_RAM.csv","ReadVariableNames",true); 
cd(cwd)

switch mode
    case 0
        disp("Processing Single Human...");
        suffix = [ condition,filesep,subj, filesep, 'Preprocessed'];
        datapath = [prefix,suffix];
        processSubject;

        % Save new datatable
        cd(tableDir)
        writetable(M, "Data_RAM.csv")
        cd(cwd)
        disp("------------------Done!------------------")

    case 2
        disp("Batch Processing every human...This might take a while!")
        cd(prefix)
        all_humans = dir(sprintf('*/S*'));

        for z = 1:numel(all_humans)
            subj = all_humans(z).name;
            cd(prefix)
            human_cond = dir(sprintf('*/%s/Preprocessed/*_RAM_223_EFR_preProcessed.mat', subj));

            if isempty(human_cond)
                continue
            else
                for folder_i = 1:size(human_cond,1)
                    cond_cell(folder_i,1) = extractBetween(human_cond(folder_i).folder,...
                        '\Human\',sprintf('%s%s%s%s', filesep, subj, filesep, 'Preprocessed'));
                end

                [conds,~,ind] = unique(cond_cell);
                cd(cwd);

                for c = 1:size(conds,1)
                    condition = conds{c};
                    suffix = [condition, filesep,subj,filesep,'Preprocessed'];
                    datapath = [prefix,suffix];
                    processSubject;
                end
            end
            fprintf("done with subject: %s\n", subj)
        end

        % Save new datatable
        cd(tableDir)
        writetable(M, "Data_RAM.csv")
        cd(cwd)
        disp("------------------Done!------------------")


    case 3
        disp("Plotting Group Data...")
        %all_humans_RAM; 
        cd(cwd)
end