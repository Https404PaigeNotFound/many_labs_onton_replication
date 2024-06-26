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
%ersp_save_dirs = {'3_ERSP_data/probe', '3_ERSP_data/memorise', '3_ERSP_data/ignore', '3_ERSP_data/maintenance'};

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
%%

% Define memory loads and letter positions % CHECK: The original instructions mention 0-5 letters
memorise_memory_loads = [3, 5, 7];
ignore_memory_loads = [4, 6, 8];  % Corresponding to memorise loads + 1
letter_positions = 0:7;  % 0 to 7 to cover all possible positions

% Function to extract event codes from EEG.event structure
function codes = extract_event_codes(EEG)
    codes = {EEG.event.type};
end

% Separate memorise and ignore data by memory load and letter position
[memorise_by_load_pos, ignore_by_load_pos] = deal(cell(length(memorise_memory_loads), length(letter_positions)));
for load_idx = 1:length(memorise_memory_loads)
    memorise_load = memorise_memory_loads(load_idx);
    ignore_load = ignore_memory_loads(load_idx);
    
    for pos_idx = 1:length(letter_positions)
        pos = letter_positions(pos_idx);
        
        % For Memorize events
        memorise_event = sprintf('s%d%d', memorise_load, pos);
        memorise_indices = find(strcmp(extract_event_codes(memorise), memorise_event));
        memorise_by_load_pos{load_idx, pos_idx} = memorise.data(:, :, memorise_indices);
        
        % For Ignore events
        ignore_event = sprintf('s%d%d', ignore_load, pos);
        ignore_indices = find(strcmp(extract_event_codes(ignore), ignore_event));
        ignore_by_load_pos{load_idx, pos_idx} = ignore.data(:, :, ignore_indices);
    end
end

% Compute ERSPs for each memory load and letter position
% ADD

% Perform regression analysis on ERSPs for each condition
for condition = 1:length(memory_load_conditions)
    for t = 1:length(times)
        for f = 1:length(freqs)
            y = cellfun(@(x) x(f, t), ersp_memorize(condition, :));
            X = (0:5)';
            b = regress(log(y), [ones(size(X)) X]);
            % Store regression coefficients or p-values as needed
        end
    end
end


% Mean ERSP Differences and Bootstrap Statistics

% Finally, longer data epochs time-locked to Maintenance period onsets and to Probe letter presentations were transformed to produce ERSP time/ frequency images using the same Fixation baseline spectrum.