
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
memoryLoads = [3, 5, 7]; % Memory load conditions
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
erspL3_means = [];
erspL5_means = [];
erspL7_means = [];
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

    %% Step A1: Compute ERSPs for Memorise Epochs (Baseline-Normalised)
    fprintf('Computing ERSPs for memorise epochs...\n');
    ERSP_baseline_norm = memorise_ERSP_baseline_norm(EEG_fixation, EEG_memorise);
    fprintf('Baseline normalisation complete.\n');


    %% Step A2: ERSP Decomposition and Template Matching
    fprintf('Performing ERSP decomposition...\n');

    % Dimensions of ERSP_baseline_norm: [channels x frequencies x time points]
    [num_channels, num_frequencies, num_time_points] = size(ERSP_baseline_norm);
    fprintf('Debug: ERSP dimensions - Channels: %d, Frequencies: %d, Time Points: %d\n', ...
        num_channels, num_frequencies, num_time_points);

    % Number of ICs
    num_ICs = size(EEG_memorise.icaact, 1); % [num ICs, timepoint/epoch, num epochs]

    % Reshape ERSP to isolate the time dimension for PCA
    reshaped_for_time_PCA = reshape(ERSP_baseline_norm, num_frequencies * num_time_points, num_ICs); % [frequencies*time points x ICs]
    %reshaped_for_time_PCA = reshape(reshaped_for_time_PCA, num_frequencies * num_time_points, num_channels); % [frequencies*time points x channels]
    
    % Debug reshaping
    fprintf('Debug: Reshaped ERSP for PCA across time - Dimensions: [%d x %d]\n', ...
        size(reshaped_for_time_PCA, 1), size(reshaped_for_time_PCA, 2));
    
    % Perform PCA across time for each frequency
    fprintf('Performing PCA across time for each frequency...\n');
    [pca_coeff_time, pca_scores_time, latent_time] = pca(reshaped_for_time_PCA);
    
    % Retain only the top components that explain the most variance
    explained_variance_threshold = 85; % Retain components explaining 85% variance  % TO DO TEST WITH 99
    cumulative_variance = cumsum(latent_time) / sum(latent_time) * 100;
    num_retained_time_components = find(cumulative_variance >= explained_variance_threshold, 1);
    
    % Ensure at least one component is retained
    if isempty(num_retained_time_components) || num_retained_time_components == 0
        num_retained_time_components = size(pca_scores_time, 2); % Retain all components
    end
     
    fprintf('Debug: Numer of num_retained_time_components: %d\n', num_retained_time_components);

    % Adjust to ensure divisibility by num_frequencies
    num_retained_time_components = max(1, ...
        floor(num_retained_time_components / num_frequencies) * num_frequencies);
    fprintf('Debug: Retaining %d components (adjusted for divisibility by %d frequencies).\n', ...
        num_retained_time_components, num_frequencies);

    fprintf('Debug: Numer of num_retained_time_components after ensuring divisibility by num_frequencies: %d\n', num_retained_time_components);
    
    % Truncate PCA scores and coefficients
    pca_scores_time = pca_scores_time(:, 1:num_retained_time_components);
    pca_coeff_time = pca_coeff_time(:, 1:num_retained_time_components);
    
    % Reshape PCA scores back into [frequencies x components x time points]
    reshaped_scores_time = reshape(pca_scores_time, num_frequencies, num_retained_time_components, num_time_points);
    
    % Combine frequency and time dimensions for ICA
    fprintf('Combining frequency and time dimensions for ICA...\n');
    reshaped_for_ICA = reshape(reshaped_scores_time, num_frequencies * num_retained_time_components, num_time_points);
    
    fprintf('Debug: Dimensions of PCA-reduced input to ICA - [%d x %d]\n', ...
        size(reshaped_for_ICA, 1), size(reshaped_for_ICA, 2));
    
    % Apply ICA on the reduced PCA dimensions
    fprintf('Performing ICA on PCA-reduced data...\n');
    [ica_weights, ica_templates] = runica(reshaped_for_ICA);
    
    % Validate ICA output dimensions
    fprintf('Debug: ICA templates dimensions before reshaping: [%d x %d]\n', size(ica_templates, 1), size(ica_templates, 2));
    
    % Automatically adjust for mismatched dimensions
    expected_elements = num_frequencies * num_retained_time_components * num_time_points;
    actual_elements = numel(ica_templates);
    
    if actual_elements ~= expected_elements
        warning('Mismatch in elements: Expected %d, but got %d. Reshaping based on actual dimensions.', ...
            expected_elements, actual_elements);
        % Dynamically adjust reshaping parameters
        adjusted_num_time_points = floor(actual_elements / (num_frequencies * num_retained_time_components));
        fprintf('Debug: Adjusted time points for reshaping: %d\n', adjusted_num_time_points);
        num_time_points = adjusted_num_time_points; % Adjust time points dynamically
    end
    
    % Reshape ICA templates
    ica_templates_reshaped = reshape(ica_templates, num_frequencies, num_retained_time_components, num_time_points);
    
    % Debug reshaped dimensions
    fprintf('Debug: ICA templates reshaped to [%d frequencies x %d components x %d time points].\n', ...
        size(ica_templates_reshaped, 1), size(ica_templates_reshaped, 2), size(ica_templates_reshaped, 3));

    % Scoring Theta Components
    theta_scores = zeros(1, size(ica_templates_reshaped, 3));
    freq_vector = linspace(freqs(1), freqs(2), num_frequencies);
    
    for ic = 1:size(ica_templates_reshaped, 3)
        % Extract ICA template from reshaped data
        ic_template = ica_templates_reshaped(:, :, ic); % [frequencies x time_points]
    
        % Identify theta band indices
        theta_indices = find(freq_vector >= thetaBand(1) & freq_vector <= thetaBand(2));
    
        % Compute theta power
        theta_power = mean(ic_template(theta_indices, :), 1); % Mean across time
        theta_scores(ic) = mean(theta_power); % Average over time
    
        % Debug theta computation
        fprintf('Debug: Theta indices: %s\n', mat2str(theta_indices));
        fprintf('Debug: Theta power for IC %d: %.4f\n', ic, theta_scores(ic));
    end
    
    % Identify dominant theta template
    [~, theta_template_idx] = max(theta_scores);
    fprintf('Selected Theta Template Index: %d\n', theta_template_idx);
    
    % Spatial Validation           
    fprintf('Performing spatial validation with DIPFIT...\n');
    EEG_memorise = performDipfit(EEG_memorise);

    fprintf('Validating IC-to-dipole mapping...\n');
    num_ICs = size(ica_weights, 1);
    num_dipoles = length(EEG_memorise.dipfit.model);
    
   if num_ICs ~= num_dipoles
        warning('Mismatch: Number of ICs (%d) does not match number of dipoles (%d). Adjusting to the smaller count.', ...
            num_ICs, num_dipoles);
        num_dipoles = min(num_ICs, num_dipoles);
   end

    % Validate dipoles based on spatial and residual variance criteria
    [valid_idx, distances] = validate_dipoles(EEG_memorise, rv_threshold, mni_constraints, acc_centroid);

    % Debug output
    fprintf('Valid dipole indices: %s\n', mat2str(find(valid_idx)));
    fprintf('Distances from ACC centroid: %s\n', mat2str(distances));
    
    % Retain spatially valid ICs:
    valid_ICs = find(valid_idx);

    % Adjust valid_idx length to match theta_scores
    if length(valid_idx) ~= length(theta_scores)
        warning('Adjusting valid_idx length to match theta_scores...');
        valid_idx = valid_idx(1:min(length(valid_idx), length(theta_scores))); % Trim to match theta_scores
    end
    fprintf('Debug: Adjusted valid indices: %s\n', mat2str(find(valid_idx)));

    
    % Check if any dipoles are valid
    if any(valid_idx)
        % Retain theta scores for spatially valid components
        valid_theta_scores = theta_scores(valid_idx);
        valid_indices = find(valid_idx);
    
        % Select the dipole with the maximum absoulte theta score
        [~, max_valid_idx] = max(abs(valid_theta_scores));
        fm_theta_idx = valid_indices(max_valid_idx);
        fprintf('Selected valid dipole (Index: %d) with max theta score: %.4f.\n', ...
            fm_theta_idx, theta_scores(fm_theta_idx));
    else
        % No valid dipoles, use closest based on distance
        fprintf('No valid dipoles found. Using fallback based on distance to ACC centroid.\n');
    
        % Find top 5 closest dipoles
        [~, sorted_indices] = sort(distances);
        top_5_closest = sorted_indices(1:min(5, length(sorted_indices)));
    
        % Retrieve theta scores for the closest dipoles
        closest_theta_scores = theta_scores(top_5_closest);
    
        % Select the dipole with the maximum absoulte theta score among the closest
        [~, max_closest_idx] = max(abs(closest_theta_scores));
        fm_theta_idx = top_5_closest(max_closest_idx);
        fprintf('Selected fallback dipole (Index: %d) with max theta score: %.4f (Distance: %.2f mm).\n', ...
            fm_theta_idx, theta_scores(fm_theta_idx), distances(fm_theta_idx));
    end
    
    % Final selected dipole index
    fprintf('Final selected dipole index: %d\n', fm_theta_idx);


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
    fprintf('Separating memorise epochs by memory load...\n');

    % Preallocate containers for memory load subsets
    memorise_by_load = cell(1, length(memoryLoads)); % memoryLoads = [3, 5, 7]
    
    % Loop through epochs to sort by memory load
    for e = 1:length(EEG_memorise.epoch)
        % Extract the first event type string
        load_type_str = EEG_memorise.epoch(e).eventtype{1}; % Get the first cell element
    
        % Extract the numeric prefix from the event type string
        if startsWith(load_type_str, 's')
            load_type = str2double(load_type_str(2)); % Extract the first digit after 's'
        else
            warning('Event type %s does not start with "s". Skipping epoch %d.', load_type_str, e);
            continue;
        end
    
        % Match the extracted load type to memoryLoads
        load_idx = find(memoryLoads == load_type, 1);
        if ~isempty(load_idx)
            memorise_by_load{load_idx} = [memorise_by_load{load_idx}, e];
        else
            warning('Unknown memory load type in epoch %d: %s', e, load_type_str);
        end
    end
    fprintf('Memorise epochs split into memory load conditions.\n');
    
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
    fprintf('Aggregating and baseline-normalising theta power...\n');
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
    csvData = [csvData; mean_theta_power]; % Append to CSV data
    erspL3_means = [erspL3_means; mean_theta_power(1)];
    erspL5_means = [erspL5_means; mean_theta_power(2)];
    erspL7_means = [erspL7_means; mean_theta_power(3)];

disp("Exiting the big loop!")
end

%{
%% Step A3: Cluster Components Across Participants
fprintf('Clustering fmθ components across participants...\n');

% Loop the dipole_info data?? 

% Perform spatial clustering around the ACC
fprintf('Performing spatial clustering around ACC...\n');

% Define cluster radius
cluster_radius = 100; % Radius in mm for clustering

% Compute distances from ACC centre
distances = pdist2(all_dipoles, acc_centroid);
cluster_idx = distances < cluster_radius; % Identify dipoles within the radius

% Retain components within the ACC cluster
final_dipoles = all_dipoles(cluster_idx, :);
final_participants = selected_fm_theta_idx(cluster_idx);

% Save or visualise the final fmθ cluster
fprintf('Final cluster includes %d components from %d participants.\n', ...
    size(final_dipoles, 1), length(final_participants));
%}

%% Export Results to CSV
% Create a table for CSV export
csvTable = table(participantIDs', erspL3_means', erspL5_means', erspL7_means', ...
    'VariableNames', {'ParticipantID', 'ERSP_L3', 'ERSP_L5', 'ERSP_L7'});
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
    'VariableNames', {'ParticipantID', 'fm_theta_idx_di', 'PosX', 'PosY', 'PosZ', 'fm_theta_rv_di'});
% Save the table as a CSV file
csvFileName_dt = fullfile(outputFolder, 'dipole_info.csv'); % Define file path
writetable(csvTableDipole, csvFileName_dt);
% Display success message
disp(['Dipole information saved as CSV: ' csvFileName_dt]);



%% Step C1: Statistical Analysis
% Organise data for repeated measures
data = [erspL3_means, erspL5_means, erspL7_means]; % Combine means for all memory loads
loadNames = {'Load3', 'Load5', 'Load7'};
participantID = (1:size(data, 1))'; % Create participant IDs

% Create table
dataTable = array2table(data, 'VariableNames', loadNames);
dataTable.ParticipantID = participantID;

% Define within-subject factors
withinDesign = table(loadNames', 'VariableNames', {'MemoryLoad'});

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
fprintf('Visualising results...\n');
figure;
errorbar(memoryLoads, mean([erspL3_means, erspL5_means, erspL7_means]), ...
    std([erspL3_means, erspL5_means, erspL7_means]) / sqrt(length(participantIDs)), 'o-');
xlabel('Memory Load (Letters)');
ylabel('Baseline-Normalised Theta Power');
title('Theta Power Across Memory Loads');
saveas(gcf, fullfile(outputFolder, 'theta_power_memory_loads.png'));



%% TO DO
%Check when ammending the length if it needs to be based on the indexs
% rather than the length


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


function candidateICs = findComponentsNearACC(EEG) % Finds dipoles in ACC using MNI Coordinate range X: -10 to +10 (midline), Y: 10 to 50 (anterior-posterior), Z: 10 to 50 (superior-inferior)
    candidateICs = [];
    for ic = 1:length(EEG.dipfit.model)
        dipole = EEG.dipfit.model(ic);
        % Check for residual variance and ACC location
        if dipole.rv < 0.2 && abs(dipole.posxyz(1)) <= 10 ... % Midline constraint
                && dipole.posxyz(2) >= 10 && dipole.posxyz(2) <= 50 ... % Anterior-posterior constraint
                && dipole.posxyz(3) >= 10 && dipole.posxyz(3) <= 50 % inferior - superior constraint
            candidateICs = [candidateICs, ic];
        end
    end
end



function ERSP_baseline_norm = memorise_ERSP_baseline_norm(EEG_fixation, EEG_memorise)
%{
This function computes the Event-Related Spectral Perturbation (ERSP) 
for all ICs in the "memorise" EEG dataset (EEG_memorise) and 
normalises it using the baseline power derived from the "fixation" 
EEG dataset (EEG_fixation). It uses the newtimef function to perform 
time-frequency decomposition on each IC, averages the ERSP across 
epochs, and outputs the normalized ERSP for all ICs as a 
3D array (channels × frequencies × time points).
%}
    % Extract data
    activity_data = EEG_memorise.icaact; % IC x Time x Epochs
    baseline_data = EEG_fixation.icaact; % IC x Time x Epochs
    
    % Convert baseline to 2D (collapse across epochs)
    %baseline_data_2d = reshape(baseline_data, size(baseline_data, 1), []); % Channels x (time*epochs)
    baseline_data_2d = reshape(baseline_data, size(baseline_data, 1), []); % ICs x (Time*Epochs)

    % Compute mean power across baseline epochs for normalisation
    %baseline_mean = mean(abs(baseline_data_2d).^2, 2); % Channels x frequencies (later applied per freq)
    baseline_mean = mean(abs(baseline_data_2d).^2, 2); % ICs x Frequencies
    
     % Initialise outputs
    %num_channels = size(activity_data, 1);
    
    % Use newtimef once on the first IC to get output dimensions
    [test_ersp, ~, ~, times, freqs] = newtimef( ...
        squeeze(activity_data(1, :, :)), ...
        EEG_memorise.pnts, ...
        [EEG_memorise.xmin, EEG_memorise.xmax] * 1000, ...
        EEG_memorise.srate, ...
        [3], ... #cycles
        'padratio', 2, ...        % Zero-padding to improve frequency resolution
        'winsize', 512, ...       % Fixed window length
        'baseline', baseline_mean(1), ... % Baseline normalisation
        'plotersp', 'off', ...
        'plotitc', 'off' ...
    );
    
    % Get actual output dimensions from newtimef
    num_freqs = size(test_ersp, 1);
    num_times = size(test_ersp, 2);
    %{
    % Preallocate output array
    ersp_all_channels = zeros(num_channels, num_freqs, num_times);
    %}

    % Number of ICs
    num_ICs = size(activity_data, 1); % Rows correspond to ICs
    
    % Preallocate output array
    ersp_all_ICs = zeros(num_ICs, num_freqs, num_times);

    % Process all channels
    %for chan = 1:num_channels

    % Process all ICs
    for ic = 1:num_ICs
        % Run newtimef for each IC
        [ersp, ~, ~, ~, ~] = newtimef( ...
            squeeze(activity_data(ic, :, :)), ...
            EEG_memorise.pnts, ...
            [EEG_memorise.xmin, EEG_memorise.xmax] * 1000, ...
            EEG_memorise.srate, ...
            [3], ... #cycles
            'padratio', 2, ...        % Zero-padding to improve frequency resolution
            'winsize', 512, ...       % Fixed window length
            'baseline', baseline_mean(ic), ... % Baseline normalisation
            'plotersp', 'off', ...
            'plotitc', 'off' ...
        );

        % Store ERSP for current channel
        %ersp_all_channels(chan, :, :) = mean(ersp, 3); % Average across epochs

        % Store ERSP for current IC
        ersp_all_ICs(ic, :, :) = mean(ersp, 3); % Average across epochs
    end

    % Return baseline-normalised ERSP
    %ERSP_baseline_norm = ersp_all_channels;
    ERSP_baseline_norm = ersp_all_ICs;
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
