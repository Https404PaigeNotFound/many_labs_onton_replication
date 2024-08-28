%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Co-Modulation of Component Activities
Correlation Analysis
Action: Compute correlations of single-trial theta ERSP template weights between component pairs.
Purpose: To analyse co-modulation of theta activity across components.
Figures: Fig. 5B shows co-modulation of theta power between component pairs.
Results Section: Component phase coherence and co-modulation.
%}

%% Co-Modulation of Component Activities

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

% Initialize variables for storing results
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
    
    % Initialize matrices for storing correlation results
    num_components = size(theta_weights, 1);
    correlation_matrix = zeros(num_components);
    significance_matrix = zeros(num_components);
    
    % Compute correlations between each pair of components
    for i = 1:num_components
        for j = i+1:num_components
            % Compute correlation between component i and component j
            correlation_matrix(i, j) = corr(theta_weights(i, :)', theta_weights(j, :)');
            correlation_matrix(j, i) = correlation_matrix(i, j); % Symmetric matrix
            
            % Bootstrap for significance testing
            bootstrap_dist = zeros(1, num_bootstraps);
            for b = 1:num_bootstraps
                % Shuffle one of the components' weights
                shuffled_weights = theta_weights(i, randperm(size(theta_weights, 2)));
                bootstrap_dist(b) = corr(shuffled_weights', theta_weights(j, :)');
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

% Visualize the correlation matrices
figure;
for p = 1:length(unique_participants)
    participant_id = unique_participants{p};
    
    subplot(ceil(sqrt(length(unique_participants))), ceil(sqrt(length(unique_participants))), p);
    imagesc(correlation_results.(participant_id));
    title(['Participant ' participant_id ' Correlation Matrix']);
    xlabel('Component');
    ylabel('Component');
    colorbar;
end

disp('Correlation analysis and visualization complete.');
