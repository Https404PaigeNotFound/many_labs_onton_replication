%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Maximally Independent Single-Trial Time/Frequency Templates
Template Identification
Action: Use PCA and ICA to identify time/frequency templates that capture independent patterns of spectral variation.
Purpose: To decompose single-trial ERSPs into independent time/frequency modes.
Figures: Fig. 3A shows the strongest time/frequency templates. Fig. 5A shows the distribution of single-trial weights for ERSP templates.
Results Section: Single-trial ERSP decomposition.
%}


% Maximally independent single-trial time/frequency templates

% NOTE: Look into using the exisiting ERSP files.... 

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contains custom functions
addpath('utils');

% Paths to each epoch
epoch_dirs = {'2_epoch_data/memorise', '2_epoch_data/fixation'};
epoch_conditions = {'memorise', 'fixation'};
% TO DO: extend out to the other relevant conditons
%epoch_dirs = {'2_epoch_data/probe', '2_epoch_data/memorise', '2_epoch_data/ignore', '2_epoch_data/fixation', '2_epoch_data/maintenance'};
%epoch_conditions = {'probe', 'memorise', 'ignore', 'fixation', 'maintenance'};

% Get list of files for each condition
file_lists = cell(length(epoch_conditions), 1);
for i = 1:length(epoch_conditions)
    epoch_files = dir(fullfile(pathToEEGLAB, epoch_dirs{i}, '*.set'));
    file_lists{i} = {epoch_files.name};
end

% Extract unique participant identifiers
unique_participants = unique(cellfun(@(x) regexp(x, '^(pilot_Manylabs_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})', 'match', 'once'), file_lists{1}, 'UniformOutput', false));

% Parameters
num_subjects = length(unique_participants);
num_templates = 12; % Number of templates to extract
num_freqs = 50; % 1-Hz intervals from 1 to 50 Hz
num_times = 200;
%num_times = 56; % Time range from -1 to 2 seconds with 0.05s intervals
freq_range = [1 50]; % Frequency range from 1 to 50 Hz
time_range = [-1000 2000]; % Time range from -1000 to 2000 milliseconds
max_components = 7; % Limit the number of components to no more than 7 for now

% Initialise a matrix to store the concatenated single-trial ERSPs
all_ersp = [];

% Function to identify anterior components
function anterior_components = get_anterior_components(EEG)
    % Define anterior channels (example: Fz, AFz, Fp1, Fp2, etc.)
    anterior_labels = {'Fz', 'AFz', 'Fp1', 'Fp2', 'F1', 'F2'};
    anterior_components = find(ismember({EEG.chanlocs.labels}, anterior_labels));
end

% Loop through each participant to compute the single-trial ERSPs and concatenate
for p_idx = 1:num_subjects
    participant_id = unique_participants{p_idx};
    
    % Load memorise and fixation data
    memorise_file = fullfile(pathToEEGLAB, epoch_dirs{1}, [participant_id, '_memorise.set']);
    fixation_file = fullfile(pathToEEGLAB, epoch_dirs{2}, [participant_id, '_fixation.set']);
    
    EEG_memorise = pop_loadset('filename', memorise_file);
    EEG_fixation = pop_loadset('filename', fixation_file);
    
    %% FOR TESTING
    % Set the number of trials you want to use
    num_trials = 10;
    
    % Limit memorise data to first 10 trials
    if EEG_memorise.trials > num_trials
        EEG_memorise = pop_select(EEG_memorise, 'trial', 1:num_trials);
        disp(['Memorise data limited to first ' num2str(num_trials) ' trials']);
    else
        disp(['Using all available memorise trials: ' num2str(EEG_memorise.trials)]);
    end
    
    % Limit fixation data to first 10 trials
    if EEG_fixation.trials > num_trials
        EEG_fixation = pop_select(EEG_fixation, 'trial', 1:num_trials);
        disp(['Fixation data limited to first ' num2str(num_trials) ' trials']);
    else
        disp(['Using all available fixation trials: ' num2str(EEG_fixation.trials)]);
    end

    %%


    % Compute ICA activations for fixation condition
    EEG_fixation = compute_ica_activations(EEG_fixation);

    % Compute ICA activations for memorise condition
    EEG_memorise = compute_ica_activations(EEG_memorise);
    
    % Get anterior components
    anterior_components = get_anterior_components(EEG_memorise);
    
    % If no anterior components found, skip this participant
    if isempty(anterior_components)
        continue;
    end

    % Limit the number of components to max_components
    anterior_components = anterior_components(1:min(max_components, length(anterior_components)));
    component_activations = EEG_memorise.icaact(anterior_components, :, :);
    num_components = length(anterior_components); % Number of anterior components


    %% FOR TESTING

    EEG_fixation.icaact = EEG_fixation.icaact(anterior_components, :, :);
    EEG_memorise.icaact = EEG_memorise.icaact(anterior_components, :, :);




    %%

   % Compute mean fixation log power baseline spectrum
    fixation_ersp = calculate_single_trial_ersps(EEG_fixation);
    % Check if fixation_ersp is empty
    if isempty(fixation_ersp)
        error('fixation_ersp is empty, cannot proceed with analysis.');
    end
    
    mean_fixation_log_power = mean(log(fixation_ersp), 3);
    disp('Mean fixation log power baseline spectrum computed for fixation');
    
    % Compute single-trial ERSPs for memorise condition
    memorise_ersp = calculate_single_trial_ersps(EEG_memorise);
    
    % Check if memorise_ersp is empty
    if isempty(memorise_ersp)
        error('memorise_ersp is empty, cannot proceed with analysis.');
    end
    
    % Debug: Display sizes of fixation_ersp and memorise_ersp
    disp('Size of fixation_ersp:');
    disp(size(fixation_ersp));
    disp('Size of memorise_ersp:');
    disp(size(memorise_ersp));
    disp('Size of mean_fixation_log_power before repmat:');
    disp(size(mean_fixation_log_power));
    
    % Ensure mean_fixation_log_power has the same dimensions as memorise_ersp
    mean_fixation_log_power = repmat(mean_fixation_log_power, [1, 1, size(memorise_ersp, 3), 1]); % Adjusting the repmat dimensions
    disp('Size of mean_fixation_log_power after repmat:');
    disp(size(mean_fixation_log_power));
    
    single_trial_ersp = log(memorise_ersp) - mean_fixation_log_power;
    disp('Mean fixation log power baseline spectrum computed for memorise');


    % Ensure proper dimension handling
    if ndims(single_trial_ersp) == 4
        ersp_batch = reshape(single_trial_ersp, num_components, num_freqs * num_times, []);
        all_ersp = [all_ersp, ersp_batch(:)]; % Preallocate and concatenate data

    else
        error('Unexpected number of dimensions in single_trial_ersp');
    end
end

% Perform PCA to reduce dimensionality
[coeff, score, ~] = pca(all_ersp', 'NumComponents', num_templates);
reduced_ersp_matrix = score';

% Perform ICA on the PCA-reduced data
[weights, sphere] = runica(reduced_ersp_matrix, 'pca', num_templates);
ica_templates = pinv(weights * sphere);

% Save the results
save('single_trial_ersp_templates.mat', 'ica_templates', 'weights', 'sphere');

disp('Maximally independent single-trial time/frequency templates calculated.');
