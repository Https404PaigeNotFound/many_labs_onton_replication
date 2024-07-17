%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% binopdf requires Statistics and Machine Learning Toolbox.

%%

% Memory load and other task comparisons

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');


%%

% Load ERSP results
load(fullfile(pathToEEGLAB, '4_memory_load', 'memory_load_results.mat'), 'ersp_results');

%%

% Define significance threshold
p_threshold = 0.01;

% Extract unique participant identifiers
unique_participants = fieldnames(ersp_results);

% Define memory loads
memory_loads = 0:5;

% Initialise variables for group ERSP statistics
group_ersp_memorise = cell(1, length(memory_loads));
group_ersp_ignore = cell(1, length(memory_loads));
group_ersp_diff = cell(1, length(memory_loads));
sig_ersp_memorise = cell(1, length(memory_loads));
sig_ersp_ignore = cell(1, length(memory_loads));
sig_ersp_diff = cell(1, length(memory_loads));

%%

% Loop through memory loads to calculate group ERSPs
for load = memory_loads
    % Initialise matrices to store ERSP data
    ersp_memorise_data = [];
    ersp_ignore_data = [];
    
    % Loop through participants to collect ERSP data
    for p = 1:length(unique_participants)
        participant_id = unique_participants{p};
        
        % Get the ERSP data for the current memory load
        ersp_memorise = ersp_results.(participant_id).memorise{load + 1};
        ersp_ignore = ersp_results.(participant_id).ignore{load + 1};
        
        % Convert the ERSP data to numeric if necessary
        if ~isnumeric(ersp_memorise.icaact_ersp)
            ersp_memorise.icaact_ersp = cell2mat(ersp_memorise.icaact_ersp);
        end
        if ~isnumeric(ersp_ignore.icaact_ersp)
            ersp_ignore.icaact_ersp = cell2mat(ersp_ignore.icaact_ersp);
        end
        
        % Store the ERSP data
        ersp_memorise_data = cat(4, ersp_memorise_data, ersp_memorise.icaact_ersp);
        ersp_ignore_data = cat(4, ersp_ignore_data, ersp_ignore.icaact_ersp);
    end
    
    % Calculate mean ERSP for memorise and ignore conditions
    group_ersp_memorise{load + 1} = mean(ersp_memorise_data, 4);
    group_ersp_ignore{load + 1} = mean(ersp_ignore_data, 4);
    
    % Calculate ERSP difference
    group_ersp_diff{load + 1} = group_ersp_memorise{load + 1} - group_ersp_ignore{load + 1};
    
    % Binomial probability test for significance across subjects
    N = size(ersp_memorise_data, 4); % Number of participants
    sig_ersp_memorise{load + 1} = zeros(size(group_ersp_memorise{load + 1}));
    sig_ersp_ignore{load + 1} = zeros(size(group_ersp_ignore{load + 1}));
    sig_ersp_diff{load + 1} = zeros(size(group_ersp_diff{load + 1}));
    
    for t = 1:size(group_ersp_memorise{load + 1}, 3)
        for f = 1:size(group_ersp_memorise{load + 1}, 2)
            k_memorise = sum(ersp_memorise_data(f, t, :, :) > 0, 4); % Number of significant memorise results
            k_ignore = sum(ersp_ignore_data(f, t, :, :) > 0, 4); % Number of significant ignore results
            k_diff = sum((ersp_memorise_data(f, t, :, :) - ersp_ignore_data(f, t, :, :)) > 0, 4); % Number of significant difference results
            
            % Calculate binomial probabilities
            P_memorise = binopdf(k_memorise, N, p_threshold);
            P_ignore = binopdf(k_ignore, N, p_threshold);
            P_diff = binopdf(k_diff, N, p_threshold);
            
            % Identify significant time/frequency points
            if P_memorise < p_threshold
                sig_ersp_memorise{load + 1}(f, t) = 1;
            end
            if P_ignore < p_threshold
                sig_ersp_ignore{load + 1}(f, t) = 1;
            end
            if P_diff < p_threshold
                sig_ersp_diff{load + 1}(f, t) = 1;
            end
        end
    end
end

%%

% Save group ERSP results
save(fullfile(pathToEEGLAB, '5_group_ERSP', 'group_ersp_results.mat'), 'group_ersp_memorise', 'group_ersp_ignore', 'group_ersp_diff', 'sig_ersp_memorise', 'sig_ersp_ignore', 'sig_ersp_diff');

%%

% Define visualization parameters
time_range = [-1 2]; % Time range in seconds
freq_range = [1 30]; % Frequency range in Hz

% Visualise group ERSP results
for load = memory_loads
    for comp = 1:size(group_ersp_memorise{load + 1}, 1)
        figure;
        
        % Group Memorise ERSP
        subplot(3, 1, 1);
        imagesc(time_range, freq_range, squeeze(group_ersp_memorise{load + 1}(comp, :, :)));
        set(gca, 'YDir', 'normal');
        title(['Group Memorise Load ' num2str(load) ' Component ' num2str(comp)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;
        hold on;
        % Check if the matrix has enough dimensions for contour
        if size(squeeze(sig_ersp_memorise{load + 1}(comp, :, :)), 1) > 1 && size(squeeze(sig_ersp_memorise{load + 1}(comp, :, :)), 2) > 1
            contour(time_range, freq_range, squeeze(sig_ersp_memorise{load + 1}(comp, :, :)), [1 1], 'LineColor', 'k');
        end
        
        % Group Ignore ERSP
        subplot(3, 1, 2);
        imagesc(time_range, freq_range, squeeze(group_ersp_ignore{load + 1}(comp, :, :)));
        set(gca, 'YDir', 'normal');
        title(['Group Ignore Load ' num2str(load) ' Component ' num2str(comp)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;
        hold on;
        % Check if the matrix has enough dimensions for contour
        if size(squeeze(sig_ersp_ignore{load + 1}(comp, :, :)), 1) > 1 && size(squeeze(sig_ersp_ignore{load + 1}(comp, :, :)), 2) > 1
            contour(time_range, freq_range, squeeze(sig_ersp_ignore{load + 1}(comp, :, :)), [1 1], 'LineColor', 'k');
        end
        
        % Group Difference ERSP
        subplot(3, 1, 3);
        imagesc(time_range, freq_range, squeeze(group_ersp_diff{load + 1}(comp, :, :)));
        set(gca, 'YDir', 'normal');
        title(['Group Difference Load ' num2str(load) ' Component ' num2str(comp)]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;
        hold on;
        % Check if the matrix has enough dimensions for contour
        if size(squeeze(sig_ersp_diff{load + 1}(comp, :, :)), 1) > 1 && size(squeeze(sig_ersp_diff{load + 1}(comp, :, :)), 2) > 1
            contour(time_range, freq_range, squeeze(sig_ersp_diff{load + 1}(comp, :, :)), [1 1], 'LineColor', 'k');
        end
    end
end

disp('Group ERSP statistics and visualisation complete.');






