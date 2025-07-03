
% Theta Power and Memory Load Analysis

%{
Input: Memorise epochs, fixation epochs
Output: Theta Power vs Memory Load results
Summary: Identify fmθ Component, Validate Near ACC, Compute ERSPs, and Perform ANOVA
%}

% Set variables
clear; clc;

% Set directories
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory
memoriseEpochFolder = fullfile(pathToEEGLAB, 'analysis_output/a_preprocessed_data/2_epoch_data/memorise'); 
fixationEpochFolder = fullfile(pathToEEGLAB, 'analysis_output/a_preprocessed_data/2_epoch_data/fixation'); 
outputFolder = fullfile(pathToEEGLAB, 'analysis_output/b_theta_power_memory'); % Output folder for preprocessed data
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% Initialise EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Analysis Parameters
thetaBand = [5 7]; % Theta frequency range (5–7 Hz)
memoryLoads = 0:5; % Accumulated memory loads L0 to L5
freqs = [2 30];
rv_threshold = 0.15; % 15% residual variance
mni_constraints = [-10 10; 10 50; 10 50]; % Midline ACC constraints
acc_centroid = [-5, 20, 25];

% Load Preprocessed Memorise and Fixation Epochs
memoriseFiles = dir(fullfile(memoriseEpochFolder, '*_memorise.set'));
fixationFiles = dir(fullfile(fixationEpochFolder, '*_fixation.set'));
if isempty(memoriseFiles) || isempty(fixationFiles)
    error('Required epoch files not found in specified folders.');
end

% Ensure the number of files match
if length(memoriseFiles) ~= length(fixationFiles)
    error('Mismatch between number of memorise and fixation files.');
end

% Preallocate variables for CSVs and ANOVA
csvData = [];
participantIDs = {};
dipole_info = [];
fm_theta_idx_di = [];
fm_theta_rv_di = [];
posX = [];
posY = [];
posZ = [];



% Load Memorise and Fixation Epochs
memoriseFiles = dir(fullfile(memoriseEpochFolder, '*_memorise.set'));
fixationFiles = dir(fullfile(fixationEpochFolder, '*_fixation.set'));

% Match files by participant ID
matchedFiles = struct(); % To store matched pairs
for i = 1:length(memoriseFiles)
    % Extract participant ID from memorise file name
    [~, memoriseName, ~] = fileparts(memoriseFiles(i).name);
    participantID = extractBefore(memoriseName, '_memorise');
    
    % Find corresponding fixation file
    matchingFixation = contains({fixationFiles.name}, participantID);
    if sum(matchingFixation) == 1
        matchedFiles(i).memoriseFile = memoriseFiles(i).name;
        matchedFiles(i).fixationFile = fixationFiles(matchingFixation).name;
    else
        error('No unique fixation file found for participant: %s', participantID);
    end
end

%% Loop Through Participant Files
for i = 1:length(matchedFiles)
    fprintf('Processing Participant %d: %s -> %s\n', i, matchedFiles(i).memoriseFile, matchedFiles(i).fixationFile);
    participantIDs{end+1} = extractBefore(matchedFiles(i).memoriseFile, '_memorise');

    % Load memorise and fixation files
    EEG_memorise = pop_loadset('filename', matchedFiles(i).memoriseFile, 'filepath', memoriseEpochFolder);
    EEG_fixation = pop_loadset('filename', matchedFiles(i).fixationFile, 'filepath', fixationEpochFolder);

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);
    EEG_fixation = ensureICAActivations(EEG_fixation);

    %% Step A1: Identify fmθ Components Using Spectral Theta Peak and DIPFIT
    fprintf('Identifying fmθ components via PSD and DIPFIT...\n');

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);

    % Run DIPFIT
    EEG_memorise = performDipfit(EEG_memorise);

    % Compute power spectrum for each IC
    num_ICs = size(EEG_memorise.icaact, 1);
    theta_scores = zeros(1, num_ICs);

    % Get power spectra (mean over all epochs)
    for ic = 1:num_ICs
        data_ic = squeeze(EEG_memorise.icaact(ic, :, :)); % time x trials
        data_ic = data_ic(:); % collapse epochs
        [pxx, f] = pwelch(data_ic, [], [], [], EEG_memorise.srate);
        theta_band = f >= thetaBand(1) & f <= thetaBand(2);
        theta_scores(ic) = mean(pxx(theta_band));
    end

    % Rank ICs by theta power
    [~, top_theta_idx] = max(theta_scores);

    % Validate dipoles for top ICs
    [valid_idx, distances] = validate_dipoles(EEG_memorise, rv_threshold, mni_constraints, acc_centroid);

    % Keep ICs with valid dipoles
    fm_theta_candidates = find(valid_idx);
    fprintf('Valid fmθ IC candidates: %s\n', mat2str(fm_theta_candidates));

    % Select best fmθ IC
    if isempty(fm_theta_candidates)
        % fallback: closest dipole
        [~, closest_idx] = min(distances);
        fm_theta_idx = closest_idx;
        fprintf('No valid ICs. Using fallback closest IC: %d\n', fm_theta_idx);
    else
        [~, best_idx] = max(theta_scores(fm_theta_candidates));
        fm_theta_idx = fm_theta_candidates(best_idx);
        fprintf('Selected fmθ IC: %d\n', fm_theta_idx);
    end


    %% Step A2: Compute Single-Trial ERSPs for fmθ and Run PCA + ICA
    fprintf('Computing ERSPs for selected fmθ component: IC %d...\n', fm_theta_idx);

    % Compute ERSPs (freqs x times x trials)
    ERSP_tensor = memorise_ERSP_baseline_norm(EEG_fixation, EEG_memorise, fm_theta_idx);

    % Reshape ERSPs: trials × (freqs × times)
    [num_freqs, num_times, num_trials] = size(ERSP_tensor);
    ERSP_matrix = reshape(permute(ERSP_tensor, [3, 1, 2]), num_trials, []);  % [trials × (freqs × times)]

    fprintf('Reshaped ERSP matrix size for PCA: [%d trials x %d features]\n', size(ERSP_matrix));

    % Run PCA
    [~, pca_scores, latent] = pca(ERSP_matrix);
    explained = cumsum(latent) / sum(latent) * 100;
    num_components = find(explained >= 85, 1);
    fprintf('PCA retained %d components to explain 85%% variance.\n', num_components);

    % Truncate PCA scores
    pca_scores_reduced = pca_scores(:, 1:num_components)';

    % Run ICA
    fprintf('Running ICA on PCA-reduced ERSPs...\n');
    [weights, sphere] = runica(pca_scores_reduced, 'extended', 1, 'stop', 1e-7);

    % Compute trial weights (contributions) for each template
    template_weights = weights * pca_scores_reduced;  % [components × trials]

    % Transpose to get [trials × components]
    template_weights = template_weights';
    fprintf('Computed trial weights for %d templates.\n', size(template_weights, 2));


    % Get independent spectral templates
    ica_templates = pinv(weights * sphere); % Each column is a template: [features × components]
    ica_templates = ica_templates'; % [components × features]

    % Optionally reshape templates back to [freqs × times]
    templates_reshaped = reshape(ica_templates, num_components, num_freqs, num_times); % [components × freqs × times]

    fprintf('Log Spectral ICA completed. Extracted %d time-frequency templates.\n', num_components);

    %% Step A3: Cluster Components Across Participants
    % (This part will be implemented after processing all participants.)
    fprintf('Clustering components will be performed after all participants are processed.\n');

    %% Step B1: Extract Component Activity
    fprintf('Extracting component activity for fmθ...\n');
    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);

    % Extract the time series activity for the identified fmθ component
    fm_theta_activity = EEG_memorise.icaact(fm_theta_idx, :, :); % Dimensions: [1 x time points x epochs]

    % Reshape the activity into [time points x epochs]
    fm_theta_activity = squeeze(fm_theta_activity); % Dimensions: [time points x epochs]
    fprintf('fmθ component activity extracted with dimensions: [%d x %d]\n', size(fm_theta_activity, 1), size(fm_theta_activity, 2));

    %% Step B2: Separate Epochs by Memory Load
    fprintf('Separating epochs by memory load...\n');
% Preallocate containers for memory load subsets and theta power means
accumulatedLoads = 0:5; % Accumulated loads (L0 to L5)
erspMeans = cell(1, length(accumulatedLoads)); % One cell per load

% Sort memorise epochs using eventAccLoad
for e = 1:length(EEG_memorise.epoch)
    if isfield(EEG_memorise.epoch(e), 'eventAccLoad')
        accLoadStr = EEG_memorise.epoch(e).eventAccLoad{1};
        loadNum = sscanf(accLoadStr, 'L%d');
        if ~isempty(loadNum) && loadNum <= 5
            memorise_by_load{loadNum+1} = [memorise_by_load{loadNum+1}, e];
        else
            warning('Unexpected AccLoad format in epoch %d: %s', e, accLoadStr);
        end
    else
        warning('Missing eventAccLoad in epoch %d', e);
    end
end
fprintf('Memorise epochs split into memory load conditions based on AccLoad.\n');

    fprintf('Assigning fixation epochs cyclically by memory load...\n');
    
    % Preallocate containers for fixation load subsets
    fixation_by_load = cell(1, length(memoryLoads));
    for e = 1:length(EEG_fixation.epoch)
        load_idx = mod(e - 1, length(memoryLoads)) + 1; % Cyclic assignment
        fixation_by_load{load_idx} = [fixation_by_load{load_idx}, e];
    end
    fprintf('Fixation epochs assigned cyclically to memory load conditions.\n');

 
    %% Step B3: Calculate Theta Power
    fprintf('Calculating theta power...\n');
    
    % Function to compute theta power for a set of epochs
    compute_power = @(activity, epochs, theta_band, srate) arrayfun(@(e) ...
        mean(mean(abs(fft(activity(:, e), [], 1)).^2, 1)), epochs);
    %{
    % Calculate theta power for memorise epochs
    theta_power_memorise = cellfun(@(epochs) ...
        compute_power(fm_theta_activity, epochs, thetaBand, EEG_memorise.srate), ...
        memorise_by_load, 'UniformOutput', false);
    %}
    %-------------
    % Calculate theta power for memorise epochs
    theta_power_memorise = cellfun(@(epochs) ...
    computeThetaPower(fm_theta_activity, epochs, thetaBand, EEG_memorise.srate), ...
    memorise_by_load, 'UniformOutput', false);

    % Calculate theta power for fixation epochs
    theta_power_fixation = cellfun(@(epochs) ...
        compute_power(fm_theta_activity, epochs, thetaBand, EEG_fixation.srate), ...
        fixation_by_load, 'UniformOutput', false);

    % Baseline-normalise theta power
    %theta_power_normalised = cellfun(@(memorise, fixation) ...
    %mean(memorise) - mean(fixation), ...
    %theta_power_memorise, theta_power_fixation);
    %-------------

    % Aggregate theta power
% Aggregate theta power for all loads (L0 to L5)
for load = 0:length(accumulatedLoads)
    mean_memorise_power = mean(theta_power_memorise{load+1});
    mean_fixation_power = mean(theta_power_fixation{load+1});
    erspMeans{load+1} = mean_memorise_power - mean_fixation_power;
    fprintf('Load %d: Mean memorise theta power = %.4f, Mean fixation theta power = %.4f, Normalised = %.4f\n', ...
        load, mean_memorise_power, mean_fixation_power, erspMeans{load+1});
end
    % Preallocate for mean theta power
    mean_theta_power = zeros(1, length(memoryLoads));
    
    % Loop through memory load conditions
    for i = 1:length(memoryLoads)
        % Aggregate mean theta power for memorise and fixation epochs
        mean_memorise_power = mean(theta_power_memorise{i});
        mean_fixation_power = mean(theta_power_fixation{i});
    
        % Baseline-normalise memorise power
        mean_theta_power(i) = mean_memorise_power - mean_fixation_power;
    
        % Debug output
        fprintf('Load %d: Mean memorise theta power = %.4f, Mean fixation theta power = %.4f, Normalised = %.4f\n', ...
            memoryLoads(i), mean_memorise_power, mean_fixation_power, mean_theta_power(i));
    end

    % Store dipole information
    fm_theta_idx_di = [fm_theta_idx_di; fm_theta_idx];
    fm_theta_rv_di = [fm_theta_rv_di; EEG_memorise.dipfit.model(fm_theta_idx).rv;];
    posX = [posX, EEG_memorise.dipfit.model(fm_theta_idx).posxyz(1)]; % Extract X-coordinates
    posY = [posY, EEG_memorise.dipfit.model(fm_theta_idx).posxyz(2)]; % Extract Y-coordinates
    posZ = [posZ, EEG_memorise.dipfit.model(fm_theta_idx).posxyz(3)]; % Extract Z-coordinates
        




    %% Save Participant-Level Data
    csvData = [csvData; mean_theta_power];

disp("Exiting the big loop!")
end

%% Export Results to CSV
% Create a table for CSV export
% Save the csv 
csvFileName = fullfile(outputFolder, 'theta_power_memory_load.csv');
writetable(csvTable, csvFileName);
disp(['CSV file saved: ' csvFileName]);

% Print sizes to console
fprintf('Size of fm_theta_idx_di: [%d x %d]\n', size(fm_theta_idx_di, 1), size(fm_theta_idx_di, 2));
fprintf('Size of fm_theta_rv_di: [%d x %d]\n', size(fm_theta_rv_di, 1), size(fm_theta_rv_di, 2));
fprintf('Size of posX: [%d x %d]\n', size(posX, 1), size(posX, 2));
fprintf('Size of posY: [%d x %d]\n', size(posY, 1), size(posY, 2));
fprintf('Size of posZ: [%d x %d]\n', size(posZ, 1), size(posZ, 2));


% Create the table
csvTableDipole = table(participantIDs', fm_theta_idx_di', posX, posY, posZ, fm_theta_rv_di', ...
% Save the table as a CSV file
csvFileName_dt = fullfile(outputFolder, 'dipole_info.csv'); % Define file path
writetable(csvTableDipole, csvFileName_dt);
% Display success message
disp(['Dipole information saved as CSV: ' csvFileName_dt]);



%% Step C1: Statistical Analysis
% Organise data for repeated measures
loadNames = {'Load0', 'Load1', 'Load2', 'Load3', 'Load4', 'Load5'};
participantID = (1:size(data, 1))'; % Create participant IDs

% Create table
dataTable.ParticipantID = participantID;

% Define within-subject factors

% Fit repeated measures model
rm = fitrm(dataTable, 'Load3-Load7~1', 'WithinDesign', withinDesign);

% Run repeated measures ANOVA
ranovaResults = ranova(rm);

% Display results
disp('Repeated Measures ANOVA Results:');
disp(ranovaResults);

% Post-hoc tests
%posthoc_results = multcompare(stats, 'CType', 'bonferroni');

%% Step C2: Visualisation
% ---- Plot theta power by memory load (from time-domain activity of fmθ) ----
fprintf('Plotting memory-load theta power summary...\n');
figure;
plot(memoryLoads, mean_theta_power, '-o', 'LineWidth', 2);
xlabel('Memory Load (Letters)');
ylabel('Baseline-Normalised Theta Power');
title('Theta Power Across Memory Loads');
grid on;
saveas(gcf, fullfile(outputFolder, 'theta_power_memory_loads.png'));

% ---- Plot ERSP patterns by template weights ----
fprintf('Visualising ERSP patterns for top/bottom 20%% of trials by template activation...\n')
% Compute template weights (activation per trial per template)
template_weights = weights * pca_scores_reduced; % [components x trials]

% Create output folder if needed
templateOutputFolder = fullfile(outputFolder, 'templates_top_bottom_trials');
if ~exist(templateOutputFolder, 'dir')
    mkdir(templateOutputFolder);
end

% Loop through each ICA template/component
for component_id = 1:num_components
    % Identify top and bottom 20% trials based on template weight
    comp_weights = template_weights(component_id, :);
    top_20_percent = comp_weights > prctile(comp_weights, 80);
    bottom_20_percent = comp_weights < prctile(comp_weights, 20);

    % Average ERSPs over selected trials
    ersp_top = mean(ERSP_tensor(:, :, top_20_percent), 3);
    ersp_bottom = mean(ERSP_tensor(:, :, bottom_20_percent), 3);

    % Plot side-by-side
    figure;
    subplot(1,2,1);
    imagesc(times, freqs, ersp_top); axis xy;
    title(sprintf('Template %d – Top 20%% Trials', component_id));
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); colorbar;

    subplot(1,2,2);
    imagesc(times, freqs, ersp_bottom); axis xy;
    title(sprintf('Template %d – Bottom 20%% Trials', component_id));
    xlabel('Time (ms)'); ylabel('Frequency (Hz)'); colorbar;

    % Save figure
    filename = fullfile(templateOutputFolder, sprintf('template%d_top_vs_bottom.png', component_id));
    saveas(gcf, filename);
    close(gcf);

    fprintf('Saved top/bottom comparison for Template %d → %s\n', component_id, filename);
end


%% Utility Functions
function EEG = ensureICAActivations(EEG)
    % Compute ICA activations if missing
    if isempty(EEG.icaact)
        fprintf('Computing ICA activations...\n');
        % Reshape data into 2D (channels x samples)
        reshapedData = reshape(EEG.data, size(EEG.data, 1), []);
        % Compute ICA activations
        icaActTemp = EEG.icaweights * EEG.icasphere * reshapedData;
        % Reshape activations back to 3D (components x time points x epochs)
        EEG.icaact = reshape(icaActTemp, size(EEG.icaweights, 1), EEG.pnts, EEG.trials);
        disp('ICA activations computed.');
    end
end

function thetaPower = computeThetaPower(activity_data, epochs, thetaBand, srate)
    % computeThetaPower Computes theta power for given epochs and theta band
    %
    % Inputs:
    %   activity_data - Time x Epochs matrix (e.g., 1500 x 609)
    %   epochs - Epoch indices to compute power for
    %   thetaBand - [min_theta max_theta], frequency range for theta
    %   srate - Sampling rate in Hz
    %
    % Output:
    %   thetaPower - Mean theta power across epochs and time
    
    % Extract epochs from activity data
    selected_data = activity_data(:, epochs); % Dimensions: Time x Selected Epochs
    
    % Compute power spectral density (PSD) using FFT
    nfft = size(selected_data, 1); % Number of FFT points (equal to time samples)
    freqs = linspace(0, srate / 2, floor(nfft / 2) + 1); % Frequency vector
    theta_range = find(freqs >= thetaBand(1) & freqs <= thetaBand(2)); % Indices for theta band
    
    % Calculate FFT for each epoch
    thetaPower_epochs = zeros(1, size(selected_data, 2)); % Preallocate theta power
    for epoch = 1:size(selected_data, 2)
        % FFT for single epoch
        fft_result = fft(selected_data(:, epoch), nfft);
        psd = abs(fft_result(1:floor(nfft / 2) + 1)).^2; % One-sided PSD
        
        % Compute mean theta power
        thetaPower_epochs(epoch) = mean(psd(theta_range));
    end
    
    % Compute mean theta power across all selected epochs
    thetaPower = mean(thetaPower_epochs);
end




function EEG = performDipfit(EEG)
    % Configure DIPFIT settings
    EEG = pop_dipfit_settings(EEG, 'hdmfile', 'standard_BEM/standard_vol.mat', ...
        'mrifile', 'standard_BEM/standard_mri.mat', ...
        'chanfile', 'standard_BEM/elec/standard_1005.elc', ...
        'coordformat', 'MNI');
    
    % Suppress FieldTrip feedback during DIPFIT
    cfg = [];
    cfg.feedback = 'no'; % Suppress feedback
    cfg.verbose = 'no';  % Suppress verbose output
    
    % Perform the fitting
    evalc('EEG = pop_multifit(EEG, 1:size(EEG.icaweights, 1), ''threshold'', 15, ''plotopt'', {''normlen'' ''on''});');
end



function ERSP_tensor = memorise_ERSP_baseline_norm(EEG_fixation, EEG_memorise, fm_theta_idx)
    % memorise_ERSP_baseline_norm: Computes ERSPs for a selected IC (e.g., fmθ)
    %
    % Inputs:
    %   EEG_fixation - EEGLAB structure (fixation)
    %   EEG_memorise - EEGLAB structure (memorise)
    %   fm_theta_idx - Index of selected fmθ component
    %
    % Output:
    %   ERSP_tensor - ERSP data for single trials: [freqs × times × trials]

    % Extract ICA activations for selected IC
    activity_data = squeeze(EEG_memorise.icaact(fm_theta_idx, :, :)); % [time × trials]
    baseline_data = squeeze(EEG_fixation.icaact(fm_theta_idx, :, :)); % [time × trials]
    
    % Flatten baseline to compute single baseline spectrum
    baseline_flat = baseline_data(:);
    baseline_mean = mean(abs(baseline_flat).^2);

    % Initialize ERSP storage
    num_trials = size(activity_data, 2);

    % Use newtimef on first trial to get dims
    [test_ersp, ~, ~, times, freqs] = newtimef( ...
        activity_data(:, 1), EEG_memorise.pnts, ...
        [EEG_memorise.xmin, EEG_memorise.xmax]*1000, ...
        EEG_memorise.srate, [3], ...
        'baseline', baseline_mean, ...
        'plotersp', 'off', 'plotitc', 'off');

    num_freqs = size(test_ersp, 1);
    num_times = size(test_ersp, 2);
    ERSP_tensor = zeros(num_freqs, num_times, num_trials);

    % Loop over trials
    for trial = 1:num_trials
        [ersp, ~, ~, ~, ~] = newtimef( ...
            activity_data(:, trial), EEG_memorise.pnts, ...
            [EEG_memorise.xmin, EEG_memorise.xmax]*1000, ...
            EEG_memorise.srate, [3], ...
            'baseline', baseline_mean, ...
            'plotersp', 'off', 'plotitc', 'off');
        ERSP_tensor(:, :, trial) = ersp;
    end
end
    


function theta_template_idx = find_theta_template(ica_templates, freqs, thetaBand, visualise)
    % find_theta_template Identifies the dominant theta template (5-7 Hz) from ICA templates
    %
    % Inputs:
    %   ica_templates - Independent components (ICs), size: [num_ICs x features]
    %   freqs - [min_freq max_freq], frequency range of the decomposition
    %   thetaBand - [min_theta max_theta], frequency range of the theta band
    %   visualise - Boolean, whether to visualise the dominant theta template
    %
    % Output:
    %   theta_template_idx - Index of the IC with dominant theta activity (5-7 Hz)
    
    % Compute the dimensions of the frequency and time axes
    num_freqs = length(freqs); % Use provided frequency vector length
    num_times = size(ica_templates, 2) / num_freqs; % Derive time dimension

    % Check validity of reshaping
    if round(num_times) ~= num_times
        error('Frequency-time features do not align with the frequency vector length.');
    end

    % Compute the frequency vector
    freq_vector = linspace(freqs(1), freqs(2), num_freqs);
    
    % Identify indices corresponding to the theta band
    theta_range = find(freq_vector >= thetaBand(1) & freq_vector <= thetaBand(2));
    
    % Initialise theta scores for each IC
    num_ICs = size(ica_templates, 1);
    theta_scores = zeros(1, num_ICs);

    % Loop through each IC
    for ic = 1:num_ICs
        % Reshape IC template back to frequency-time representation
        ic_template = reshape(ica_templates(ic, :), num_freqs, num_times);
        
        % Extract power in the theta band
        theta_power = mean(ic_template(theta_range, :), 1); % Mean over theta frequencies
        theta_scores(ic) = mean(theta_power); % Average over time
    end

    % Identify the IC with the maximum theta score
    [~, theta_template_idx] = max(theta_scores);

    % Visualise the dominant theta template if requested
    if visualise
        dominant_theta_template = reshape(ica_templates(theta_template_idx, :), num_freqs, num_freqs);
        times = linspace(0, 1, num_freqs); % Example time vector (adjust as needed)
        
        figure;
        imagesc(times, freq_vector, dominant_theta_template);
        axis xy;
        colorbar;
        title('Dominant Theta Template (5–7 Hz)');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
    end
end



function theta_scores = score_theta_alignment(ica_templates, theta_template_idx, freqs, thetaBand)
    % score_theta_alignment Scores ICs based on alignment with the theta template
    %
    % Inputs:
    %   ica_templates - IC templates (num_ICs x features)
    %   theta_template_idx - Index of the dominant theta template
    %   freqs - Frequency vector (based on original ERSP dimensions)
    %   thetaBand - [min_theta max_theta], frequency range of the theta band

    % Recreate frequency vector based on the original ERSP dimensions
    num_features = size(ica_templates, 2);
    freq_vector = linspace(freqs(1), freqs(2), num_features);

    % Identify indices corresponding to the theta band
    theta_indices = find(freq_vector >= thetaBand(1) & freq_vector <= thetaBand(2));
    fprintf('Debug: Theta band indices: %s\n', mat2str(theta_indices));

    % Extract and normalise theta template
    theta_template = ica_templates(theta_template_idx, :);
    theta_template_norm = theta_template / norm(theta_template);

    % Calculate theta scores
    num_ICs = size(ica_templates, 1);
    theta_scores = zeros(1, num_ICs);
    for ic = 1:num_ICs
        % Extract IC template
        ic_template = ica_templates(ic, :);
        ic_template_norm = ic_template / norm(ic_template);

        % Compute theta alignment as dot product (cosine similarity)
        theta_scores(ic) = dot(theta_template_norm(theta_indices), ic_template_norm(theta_indices));
    end
end


function [valid_idx, distances] = validate_dipoles(EEG, rv_threshold, mni_constraints, acc_centroid)
    % Validate dipoles based on residual variance and MNI constraints, and calculate distances from ACC centroid
    % Inputs:
    %   EEG - EEG structure with dipfit model
    %   rv_threshold - Residual variance threshold (e.g., 0.15)
    %   mni_constraints - [X_min X_max; Y_min Y_max; Z_min Z_max] matrix defining MNI constraints
    %   acc_centroid - [x, y, z], ACC centroid coordinates (e.g., [-5, 20, 25])
    % Outputs:
    %   valid_idx - Logical array indicating valid ICs (1 for valid, 0 for invalid)
    %   distances - Array of distances from ACC centroid for all ICs

    num_ICs = length(EEG.dipfit.model);
    valid_idx = false(1, num_ICs); % Initialise all ICs as invalid
    distances = inf(1, num_ICs); % Initialise distances array with Inf

    for ic = 1:num_ICs
        dipole = EEG.dipfit.model(ic);

        % Skip invalid dipoles
        if isempty(dipole.posxyz) || isempty(dipole.rv)
            distances(ic) = Inf;
            continue;
        end

        % Compute Euclidean distance from ACC centroid
        posxyz = dipole.posxyz;
        distances(ic) = norm(posxyz - acc_centroid);

        % Check if the dipole meets residual variance threshold
        if dipole.rv > rv_threshold
            continue; % Skip if residual variance is too high
        end

        % Check MNI coordinate constraints
        if posxyz(1) >= mni_constraints(1, 1) && posxyz(1) <= mni_constraints(1, 2) && ... % X-axis
           posxyz(2) >= mni_constraints(2, 1) && posxyz(2) <= mni_constraints(2, 2) && ... % Y-axis
           posxyz(3) >= mni_constraints(3, 1) && posxyz(3) <= mni_constraints(3, 2)       % Z-axis
            valid_idx(ic) = true; % Mark as valid
        end
    end

    % Debugging output for distances and valid indices
    fprintf('Debug: Distances from ACC centroid: %s\n', mat2str(distances));
    fprintf('Debug: Valid indices: %s\n', mat2str(find(valid_idx)));

    % Handle mismatched IC and dipole counts
    if length(valid_idx) ~= length(distances)
        warning('Mismatch: Number of ICs (%d) does not match number of dipoles (%d). Adjusting to smaller count.', ...
                length(valid_idx), length(distances));
        valid_idx = valid_idx(1:min(length(valid_idx), length(distances))); % Adjust valid_idx length
        distances = distances(1:min(length(valid_idx), length(distances))); % Adjust distances length
    end
end
