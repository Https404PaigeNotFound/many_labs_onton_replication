%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Co-Modulation of Component Activities
Clustering and Correlation Analysis
Action: 
    1. Cluster components into the frontal theta cluster based on spatial and frequency characteristics.
    2. Compute correlations of single-trial theta ERSP template weights between component pairs within the frontal theta cluster.
Purpose: 
    To identify and analyse co-modulation of theta activity across components within the frontal theta cluster.
Figures: 
    Fig. 5B shows co-modulation of theta power between component pairs.
Results Section: Component phase coherence and co-modulation.
%}

%% Co-Modulation of Component Activities with Clustering

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% Paths to each ERSP template file
ersp_template_filepath = fullfile(pathToEEGLAB, '4_memory_load', 'memory_load_results.mat');

% Load ERSP results
load(ersp_template_filepath, 'ersp_results');

% Parameters for correlation analysis
num_bootstraps = 200;
correlation_threshold = 0.05; % Set your significance threshold here

% Initialise variables for storing results
correlation_results = struct();
significance_results = struct();

% Loop through each participant in the ERSP results
unique_participants = fieldnames(ersp_results);

for p = 1:length(unique_participants)
    participant_id = unique_participants{p};
    
    % Extract theta ERSP template weights for this participant
    theta_weights_memorise = ersp_results.(participant_id).memorise_theta_weights; % Assuming this field contains the theta weights
    theta_weights_ignore = ersp_results.(participant_id).ignore_theta_weights;
    
    % Combine theta weights from both conditions
    theta_weights = [theta_weights_memorise; theta_weights_ignore];
    
    % Load the corresponding EEG dataset to access component information
    eeg_filename = ersp_results.(participant_id).eeg_filename; % Assuming this field contains the EEG filename
    EEG = pop_loadset('filename', eeg_filename, 'filepath', ersp_results.(participant_id).eeg_filepath); % Adjust fields as necessary
    
    % Perform Component Clustering
    % Criteria: Frontal topography and significant theta activity
    % Define frontal electrodes (adjust based on your montage)
    frontal_electrodes = {'Fz', 'FCz', 'Fp1', 'Fp2', 'AFz', 'AF3', 'AF4'};
    
    % Identify components with maximal activity at frontal electrodes
    is_frontal = false(1, EEG.nbchan);
    for fe = 1:length(frontal_electrodes)
        electrode_idx = find(strcmp({EEG.chanlocs.labels}, frontal_electrodes{fe}));
        if ~isempty(electrode_idx)
            is_frontal(electrode_idx) = true;
        end
    end
    
    % Compute Global Field Power (GFP) or other spatial criteria
    % Here, we'll use the component's topography correlation with frontal electrodes
    frontal_components = [];
    for comp = 1:EEG.nbchan
        topography = EEG.icawinv(:, comp);
        topography_frontal = topography(is_frontal);
        if any(abs(topography_frontal) > 0.1) % Threshold can be adjusted
            frontal_components = [frontal_components, comp];
        end
    end
    
    % Compute Frequency Criteria using previously calculated ERSP
    % Assuming ERSP is computed per component and frequency
    % Here, we'll define components as frontal theta if their mean theta power is above a threshold
    theta_power = mean(theta_weights(:, ersp_results.(participant_id).theta_freq_indices), 2); % Adjust 'theta_freq_indices'
    theta_threshold = prctile(theta_power, 75); % Top 25% as an example
    frontal_theta_components = frontal_components(theta_power(frontal_components) >= theta_threshold);
    
    % If no components meet both criteria, skip
    if isempty(frontal_theta_components)
        warning(['No frontal theta components found for participant ' participant_id]);
        continue;
    end
    
    % Extract theta weights for frontal theta components
    frontal_theta_weights = theta_weights(frontal_theta_components, :);
    
    % Initialise matrices for storing correlation results
    num_components = size(frontal_theta_weights, 1);
    correlation_matrix = zeros(num_components);
    significance_matrix = zeros(num_components);
    
    % Compute correlations between each pair of frontal theta components
    for i = 1:num_components
        for j = i+1:num_components
            % Compute correlation between component i and component j
            correlation_matrix(i, j) = corr(frontal_theta_weights(i, :)', frontal_theta_weights(j, :)');
            correlation_matrix(j, i) = correlation_matrix(i, j); % Symmetric matrix
            
            % Bootstrap for significance testing
            bootstrap_dist = zeros(1, num_bootstraps);
            for b = 1:num_bootstraps
                % Shuffle one of the components' weights
                shuffled_weights = frontal_theta_weights(i, randperm(size(frontal_theta_weights, 2)));
                bootstrap_dist(b) = corr(shuffled_weights', frontal_theta_weights(j, :)');
            end
            
            % Determine significance based on the bootstrap distribution
            p_value = sum(abs(bootstrap_dist) >= abs(correlation_matrix(i, j))) / num_bootstraps;
            significance_matrix(i, j) = p_value < correlation_threshold;
            significance_matrix(j, i) = significance_matrix(i, j); % Symmetric matrix
        end
    end
    
    % Store results in the structure
    correlation_results.(participant_id) = correlation_matrix;
    significance_results.(participant_id) = significance_matrix;
end

% Save the correlation and significance results
save(fullfile(pathToEEGLAB, '5_correlations', 'correlation_results.mat'), 'correlation_results', 'significance_results');

% Visualise the correlation matrices
figure;
for p = 1:length(unique_participants)
    participant_id = unique_participants{p};
    
    if ~isfield(correlation_results, participant_id)
        continue; % Skip participants without frontal theta components
    end
    
    subplot(ceil(sqrt(length(unique_participants))), ceil(sqrt(length(unique_participants))), p);
    imagesc(correlation_results.(participant_id));
    title(['Participant ' participant_id ' Correlation Matrix']);
    xlabel('Component');
    ylabel('Component');
    colorbar;
end

sgtitle('Co-Modulation of Frontal Theta Components Across Participants');

disp('Component clustering, correlation analysis, and visualisation complete.');
