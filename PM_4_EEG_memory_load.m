%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Memory Load and Task Comparisons
Memory Load Separation and ERSP Regression
Action: Separate Memorise and Ignore letters by memory load and compute ERSPs for each condition. Regress ERSPs on memory load at each time/frequency point.
Purpose: To study how spectral power evolves with increasing memory load.
Figures: Fig. 6 shows memory load and task comparisons in ERSPs.
Results Section: Mean spectral power change.
%}
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

% Define the baseline period (in seconds)
baseline_window = [0, 4]; % First 4 seconds of the fixation period

% Define the window length for FFT and step size (in samples)
ersp_window = 256; % Window length in samples
ersp_step = 128; % Step size in samples (50% overlap). Source material does not spec, value chosen to balance temporal resolution and computational efficiency.

% Number of bootstrap resamples for statistical significance
num_bootstraps = 1000;

% Flag to indicate whether to visualise the results
visualise = true;

% Define memory loads
memory_loads = 0:5;

%%

% Paths to each epoch
epoch_dirs = {'2_epoch_data/probe', '2_epoch_data/memorise', '2_epoch_data/ignore', '2_epoch_data/fixation', '2_epoch_data/maintenance'};
epoch_conditions = {'probe', 'memorise', 'ignore', 'fixation', 'maintenance'};
ersp_save_dirs = {'4_ERSP_data/probe', '4_ERSP_data/memorise', '4_ERSP_data/ignore', '4_ERSP_data/maintenance'};

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

% Initialise a structure to store ERSP results
ersp_results = struct();

%%

% Loop through each participant
for p = 1:length(unique_participants)
    participant_id = unique_participants{p};
    
    % Initialise a structure to store the loaded EEG datasets
    EEG_data = struct();
    
    % Load the fixation baseline for the participant
    fixation_file = fullfile(pathToEEGLAB, epoch_dirs{4}, [participant_id, '_fixation.set']);
    EEG_fixation = pop_loadset('filename', fixation_file);
    EEG_fixation = compute_ica_activations(EEG_fixation);
    
    % Load the corresponding files for the participant
    for c = 1:length(epoch_conditions)
        condition = epoch_conditions{c};
        condition_file = fullfile(pathToEEGLAB, epoch_dirs{c}, [participant_id, '_', condition, '.set']);
        EEG_data.(condition) = pop_loadset('filename', condition_file);
    end
    
    % Initialise arrays to store ERSPs for Memorise and Ignore conditions
    ERSP_memorise = cell(1, length(memory_loads));
    ERSP_ignore = cell(1, length(memory_loads));
    
    % Calculate ERSPs for Memorise and Ignore conditions
    for load = memory_loads
        % Memorise condition
        EEG_memorise = filter_epochs_by_load(EEG_data.memorise, load); % Custom function to filter epochs by memory load
        EEG_memorise = compute_ica_activations(EEG_memorise);
        ERSP_memorise{load + 1} = calculate_ersps(EEG_memorise, EEG_fixation, baseline_window, ersp_window, ersp_step, num_bootstraps, visualise);
        
        % Ignore condition
        EEG_ignore = filter_epochs_by_load(EEG_data.ignore, load); % Custom function to filter epochs by memory load
        EEG_ignore = compute_ica_activations(EEG_ignore);
        ERSP_ignore{load + 1} = calculate_ersps(EEG_ignore, EEG_fixation, baseline_window, ersp_window, ersp_step, num_bootstraps, visualise);
    end
    
    % Sanitise participant ids remove - _
    san_participant_id = strrep(participant_id, '-', '_');


    % Store ERSP results for the participant
    ersp_results.(san_participant_id).memorise = ERSP_memorise;
    ersp_results.(san_participant_id).ignore = ERSP_ignore;
    
    % Save the ERSP results for Memorise and Ignore conditions
    save_ersp_filepath = fullfile(pathToEEGLAB, '4_memory_load', [participant_id, '_memory_load_ersp.set']);
    pop_saveset(EEG_data.memorise, 'filename', save_ersp_filepath);
    
    disp(['Memory load ERSPs calculated and saved for ' participant_id]);
end

% Save the ERSP results
save(fullfile(pathToEEGLAB, '4_memory_load', 'memory_load_results.mat'), 'ersp_results');

% Check if there are multiple participants
if length(unique_participants) > 1
    % Perform statistical analysis across participants
    disp('Performing statistical analysis across participants...');
    
    % Initialise arrays to store grand average ERSPs
    grand_avg_ersp_memorise = cell(1, length(memory_loads));
    grand_avg_ersp_ignore = cell(1, length(memory_loads));
    
    % Calculate grand average ERSPs
    for load = memory_loads
        % Concatenate ERSPs across participants
        concatenated_ersp_memorise = cat(4, ersp_results.(unique_participants{1}).memorise{load + 1}.icaact_ersp{:});
        concatenated_ersp_ignore = cat(4, ersp_results.(unique_participants{1}).ignore{load + 1}.icaact_ersp{:});
        
        for p = 2:length(unique_participants)
            concatenated_ersp_memorise = cat(4, concatenated_ersp_memorise, ersp_results.(unique_participants{p}).memorise{load + 1}.icaact_ersp{:});
            concatenated_ersp_ignore = cat(4, concatenated_ersp_ignore, ersp_results.(unique_participants{p}).ignore{load + 1}.icaact_ersp{:});
        end
        
        % Calculate mean ERSPs across participants
        grand_avg_ersp_memorise{load + 1} = mean(concatenated_ersp_memorise, 4);
        grand_avg_ersp_ignore{load + 1} = mean(concatenated_ersp_ignore, 4);
    end
    
    % Perform regression analysis on grand average ERSPs
    for load = memory_loads
        for comp = 1:size(grand_avg_ersp_memorise{load + 1}, 1)
            for f = 1:size(grand_avg_ersp_memorise{load + 1}, 2)
                for t = 1:size(grand_avg_ersp_memorise{load + 1}, 3)
                    y_memorise = squeeze(grand_avg_ersp_memorise{load + 1}(comp, f, t, :));
                    y_ignore = squeeze(grand_avg_ersp_ignore{load + 1}(comp, f, t, :));
                    X = (0:5)';
                    b_memorise = regress(log(y_memorise), [ones(size(X)) X]);
                    b_ignore = regress(log(y_ignore), [ones(size(X)) X]);
                    % Store regression coefficients or p-values as needed
                end
            end
        end
    end
    
    % Perform statistical comparison between memorise and ignore conditions
    disp('Performing statistical comparison between memorise and ignore conditions...');
    
    % Mean ERSP Differences and Bootstrap Statistics
    for load = memory_loads
        for comp = 1:size(grand_avg_ersp_memorise{load + 1}, 1)
            mean_diff_ersp = grand_avg_ersp_memorise{load + 1}(comp, :, :) - grand_avg_ersp_ignore{load + 1}(comp, :, :);
            
            % Bootstrap significance testing
            bootstrap_distributions = zeros(num_bootstraps, size(mean_diff_ersp, 2), size(mean_diff_ersp, 3));
            for b = 1:num_bootstraps
                resample_indices_memorise = randi(size(grand_avg_ersp_memorise{load + 1}, 4), [1 size(grand_avg_ersp_memorise{load + 1}, 4)]);
                resample_indices_ignore = randi(size(grand_avg_ersp_ignore{load + 1}, 4), [1 size(grand_avg_ersp_ignore{load + 1}, 4)]);
                resampled_diff_ersp = mean(grand_avg_ersp_memorise{load + 1}(:, :, :, resample_indices_memorise), 4) - mean(grand_avg_ersp_ignore{load + 1}(:, :, :, resample_indices_ignore), 4);
                bootstrap_distributions(b, :, :) = resampled_diff_ersp;
            end
            
            bootstrap_means = mean(bootstrap_distributions, 1);
            bootstrap_stds = std(bootstrap_distributions, 0, 1);
            significance_threshold = squeeze(bootstrap_means + 2.58 * bootstrap_stds); % 99% confidence interval
            
            % Identify significant points of difference
            sig_diff_mask = mean_diff_ersp > significance_threshold;
            ersp_results.mean_diff_ersp{load + 1} = mean_diff_ersp;
            ersp_results.mean_diff_ersp_sig{load + 1} = sig_diff_mask;
        end
    end
else
    % Visualise single participant data
    disp('Only one participant found, visualising data without statistical analysis...');
    participant_id = unique_participants{1};
    san_participant_id = strrep(participant_id, '-', '_'); % Sanitise participant ID
    for load = memory_loads
        for comp = 1:length(ersp_results.(san_participant_id).memorise)
            figure;
            subplot(2, 1, 1);
            imagesc(ersp_results.(san_participant_id).memorise{load + 1}.times, ersp_results.(san_participant_id).memorise{load + 1}.freqs, ersp_results.(san_participant_id).memorise{load + 1}.icaact_ersp{comp});
            set(gca, 'YDir', 'normal');
            title(['Memorise Load ' num2str(load) ' Component ' num2str(comp)]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
    
            subplot(2, 1, 2);
            imagesc(ersp_results.(san_participant_id).ignore{load + 1}.times, ersp_results.(san_participant_id).ignore{load + 1}.freqs, ersp_results.(san_participant_id).ignore{load + 1}.icaact_ersp{comp});
            set(gca, 'YDir', 'normal');
            title(['Ignore Load ' num2str(load) ' Component ' num2str(comp)]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
        end
    end
end

disp('All memory load ERSP calculations and visualisations complete.');