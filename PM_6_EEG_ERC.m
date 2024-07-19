%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Event-Related Component Phase Coherence
Phase Coherence Calculation and Bootstrap Statistics
Action: Compute phase coherence between component pairs using a complex sinusoidal wavelet. Use bootstrap methods to estimate the significance of coherence increases.
Purpose: To examine fixed phase relationships between component processes.
Figures: Fig. 5B shows co-modulation of theta power between component pairs.
%}
%%

% Event-related component phase coherence

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');

%%

% Paths to each epoch
epoch_dirs = {'2_epoch_data/probe', '2_epoch_data/memorise', '2_epoch_data/ignore', '2_epoch_data/fixation', '2_epoch_data/maintenance'};
epoch_conditions = {'probe', 'memorise', 'ignore', 'fixation', 'maintenance'};

%%

% Get list of files for each condition
file_lists = cell(length(epoch_conditions), 1);
for i = 1:length(epoch_conditions)
    epoch_files = dir(fullfile(pathToEEGLAB, epoch_dirs{i}, '*.set'));
    file_lists{i} = {epoch_files.name};
end

%%

% Extract unique participant identifiers
unique_participants = unique(cellfun(@(x) regexp(x, '^(pilot_Manylabs_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})', 'match', 'once'), file_lists{1}, 'UniformOutput', false));

% Parameters for wavelet transformation
frequencies = linspace(3, 50, 100); % Frequency range from 3 Hz to 50 Hz
num_cycles = linspace(3, 6, length(frequencies)); % Number of cycles increases smoothly from 3 to 6

% Define conditions
condition_names = {'fixation', 'memorise', 'ignore', 'maintenance', 'probe'};
num_conditions = length(condition_names);

% Initialise variables for ERC calculation
erc_values = cell(num_conditions, 1);
significant_erc = cell(num_conditions, 1);
num_bootstraps = 1000;
p_threshold = 0.01;

% Define the wavelet parameters
wavelet_params = struct('frequencies', frequencies, 'num_cycles', num_cycles);

% Loop through each participant
for p_idx = 1:length(unique_participants)
    participant_id = unique_participants{p_idx};
    
    % Initialise a structure to store the loaded EEG datasets
    EEG_data = struct();
    
    % Load the corresponding files for the participant
    for c = 1:length(epoch_conditions)
        condition = epoch_conditions{c};
        condition_file = fullfile(pathToEEGLAB, epoch_dirs{c}, [participant_id, '_', condition, '.set']);
        EEG_data.(condition) = pop_loadset('filename', condition_file);
    end
    
    % Compute ICA activations for each condition
    component_activations = struct();
    for c_idx = 1:num_conditions
        condition = condition_names{c_idx};
        component_activations.(condition) = compute_ica_activations(EEG_data.(condition));
    end
    
    % Calculate ERC for each condition
    for cond_idx = 1:num_conditions
        condition = condition_names{cond_idx};
        
        % Get component activations for the condition
        component_data = component_activations.(condition);
        num_components = size(component_data, 1);
        num_trials = size(component_data, 2);
        num_timepoints = size(component_data, 3);
        
        % Initialise ERC values
        erc_matrix = zeros(num_components, num_components, length(frequencies), num_timepoints);
        
        % Loop through pairs of components
        for comp1 = 1:num_components
            for comp2 = comp1+1:num_components
                % Calculate phase coherence for each pair
                [coherence, ~] = compute_phase_coherence(squeeze(component_data(comp1, :, :)), squeeze(component_data(comp2, :, :)), wavelet_params);
                
                % Store ERC values
                erc_matrix(comp1, comp2, :, :) = coherence;
                erc_matrix(comp2, comp1, :, :) = coherence; % Symmetric matrix
            end
        end
        
        % Store ERC values for the condition
        erc_values{cond_idx} = erc_matrix;
        
        % Calculate significance using bootstrap method
        significant_erc{cond_idx} = calculate_erc_significance(erc_matrix, component_activations.fixation, wavelet_params, num_bootstraps, p_threshold);
        
        disp(['ERC calculated for ' condition ' for participant ' participant_id]);

    end
    
    % Save the ERC results for the participant
    save_filename = fullfile(pathToEEGLAB, '6_ERC', [participant_id, '_erc_results.mat']);
    save(save_filename, 'erc_values', 'significant_erc');
    
    % Visualise significant ERC results
    %visualise_erc_results(significant_erc, wavelet_params, condition_names, participant_id);
end

disp('ERC calculation and significance testing complete.');


