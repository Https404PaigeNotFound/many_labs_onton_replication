%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% Paths to each epoch
epoch_dirs = {'2_epoch_data/probe', '2_epoch_data/memorise', '2_epoch_data/ignore', '2_epoch_data/fixation', '2_epoch_data/maintenance'};
epoch_conditions = {'probe', 'memorise', 'ignore', 'fixation', 'maintenance'};
ersp_save_dirs = {'3_ERSP_data/probe', '3_ERSP_data/memorise', '3_ERSP_data/ignore', '3_ERSP_data/maintenance'};

%%

% Define the baseline period (in seconds)
baseline_window = [0, 4]; % First 4 seconds of the fixation period

% Define the window length for FFT and step size (in samples)
ersp_window = 256; % Window length in samples
ersp_step = 128; % Step size in samples (50% overlap). Source material does not spec, value chosen to balance temporal resolution and computational efficiency.

% Number of bootstrap resamples for statistical significance
num_bootstraps = 1000;

% Flag to indicate whether to visualize the results
visualize = true;

% Get list of files for each condition
file_lists = cell(length(epoch_conditions), 1);
for i = 1:length(epoch_conditions)
    epoch_files = dir(fullfile(pathToEEGLAB, epoch_dirs{i}, '*.set'));
    file_lists{i} = {epoch_files.name};
end

% Extract unique participant identifiers
unique_participants = unique(cellfun(@(x) regexp(x, '^(pilot_Manylabs_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})', 'match', 'once'), file_lists{1}, 'UniformOutput', false));

% Loop through each participant
for p = 1:length(unique_participants)
    participant_id = unique_participants{p};
    
    % Initialize a structure to store the loaded EEG datasets
    EEG_data = struct();
    
    % Load the fixation baseline for the participant
    fixation_file = fullfile(pathToEEGLAB, '2_epoch_data/fixation', [participant_id, '_fixation.set']);
    EEG_fixation = pop_loadset('filename', fixation_file);

    % Loop through each condition and load the corresponding files for the participant
    for c = 1:length(epoch_conditions)
        condition = epoch_conditions{c};
        condition_file = fullfile(pathToEEGLAB, epoch_dirs{c}, [participant_id, '_', condition, '.set']);
        EEG_data.(condition) = pop_loadset('filename', condition_file);
    end
    
    % Calculate ERSPs for each condition using the fixation baseline
    for c = 1:length(epoch_conditions)
        if strcmp(epoch_conditions{c}, 'fixation')
            continue; % Skip fixation as we use it as baseline
        end
        
        condition = epoch_conditions{c};
        EEG = EEG_data.(condition);
        EEG_fixation = compute_ica_activations(EEG_fixation); % Compute independent component activations with custom function
        EEG = compute_ica_activations(EEG); % Compute independent component activations with custom function
        EEG = calculate_ersps(EEG, EEG_fixation, baseline_window, ersp_window, ersp_step, num_bootstraps, visualize);
        
        % Define the save directory explicitly
        switch condition
            case 'probe'
                save_dir = fullfile(pathToEEGLAB, '3_ERSP_data/probe');
            case 'memorise'
                save_dir = fullfile(pathToEEGLAB, '3_ERSP_data/memorise');
            case 'ignore'
                save_dir = fullfile(pathToEEGLAB, '3_ERSP_data/ignore');
            case 'maintenance'
                save_dir = fullfile(pathToEEGLAB, '3_ERSP_data/maintenance');
        end
        
        % Check if the directory exists, if not, create it
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        % Save the ERSP results
        save_ersp_filepath = fullfile(save_dir, [participant_id, '_', condition, '_ersp.set']);
        pop_saveset(EEG, 'filename', save_ersp_filepath);
        
        disp(['ERSPs calculated and saved for ' participant_id ' in ' condition ' condition']);
    end
end

disp('All ERSP calculations complete.');
