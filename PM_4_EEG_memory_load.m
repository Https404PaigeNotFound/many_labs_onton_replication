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

% Path to preprocessed data 
preprocessed_dir = '1_preprocessed_data';

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


%%

% Get list of preprocessed files
preprocessed_files = dir(fullfile(pathToEEGLAB, preprocessed_dir, '*.set'));






%%

% Loop over all participants 
for subj = 1:length(preprocessed_files)
    % Load subject-specific epoch data
    EEG = pop_loadset('filename', preprocessed_files(subj).name, 'filepath', fullfile(preprocessed_files(subj).folder));
    %[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);
    disp('Data loaded')

   
end

disp('End of loop')


%%

% Re-labeling triggers happens in the preprocessng now 


% Compute for a whole (combine 3, 5 & 7). 
% seperate out the difficulty 3, 5, & 7. 



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