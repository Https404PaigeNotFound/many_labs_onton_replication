


%WE ARE NOW DOING A ANOVA IV:2(condition)x3(memory load) DV: Baseline-normalised alpha power from the alpha component. 




% alpha Power and Memory Load Analysis

%{
Input: Memorise epochs, ignore epochs, fixation epochs
Output: Condition (Memorise, Ignore), Memory load (3,5,7) and alpha Power 
Summary: Identify parietal alpha Component, Validate Near PPC, Compute ERSPs, and Perform ANOVA
%}

% Set variables
clear; clc;

% Set directories
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory
memoriseEpochFolder = fullfile(pathToEEGLAB, 'analysis_output/a_preprocessed_data/2_epoch_data/memorise'); 
ignoreEpochFolder = fullfile(pathToEEGLAB, 'analysis_output/a_preprocessed_data/2_epoch_data/ignore'); 
fixationEpochFolder = fullfile(pathToEEGLAB, 'analysis_output/a_preprocessed_data/2_epoch_data/fixation'); 
outputFolder = fullfile(pathToEEGLAB, 'analysis_output/c_alpha_power_memory_load_memorise_ignore'); % Output folder for preprocessed data
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
alphaBand = [8 12]; % Alpha frequency range (8–12 Hz)
freqs = [2 30];
rv_threshold = 0.15; % 15% residual variance
mni_constraints = [-30 50; -70 -50; 30 50]; % posterior parietal cortex  constraints
ppc_centroid = [-48, -66, 34];

% Load Preprocessed Memorise and Fixation Epochs
memoriseFiles = dir(fullfile(memoriseEpochFolder, '*_memorise.set'));
fixationFiles = dir(fullfile(fixationEpochFolder, '*_fixation.set'));
ignoreFiles = dir(fullfile(ignoreEpochFolder, '*_ignore.set'));
if isempty(memoriseFiles) || isempty(fixationFiles) || isempty(ignoreFiles)
    error('Required epoch files not found in specified folders.');
end

% Ensure the number of files match
if length(memoriseFiles) ~= length(fixationFiles) || length(memoriseFiles) ~= length(ignoreFiles)
    error('Mismatch between number of memorise and fixation/ignore files.');
end

% Preallocate for dynamic loads L0 to L5
accumulatedLoads = 0:5;
memorise_erspMeans = cell(1, length(accumulatedLoads));
ignore_erspMeans = cell(1, length(accumulatedLoads));

fm_alpha_idx_di = [];
fm_alpha_rv_di = [];
posX = [];
posY = [];
posZ = [];




% Match files by participant ID
matchedFiles = struct(); % To store matched pairs
for i = 1:length(memoriseFiles)
    % Extract participant ID from memorise file name
    [~, memoriseName, ~] = fileparts(memoriseFiles(i).name);
    participantID = extractBefore(memoriseName, '_memorise');
    
    % Find corresponding files
    matchingFixation = contains({fixationFiles.name}, participantID);
    if sum(matchingFixation) == 1
        matchedFiles(i).memoriseFile = memoriseFiles(i).name;
        matchedFiles(i).fixationFile = fixationFiles(matchingFixation).name;
        matchedFiles(i).ignoreFile = ignoreFiles(matchingIgnore).name;
    else
        error('No unique fixation or ignore file found for participant: %s', participantID);
    end
end

%% Loop Through Participant Files
for i = 1:length(matchedFiles)
    fprintf('Processing Participant %d: %s -> %s %s\n', i, matchedFiles(i).memoriseFile, matchedFiles(i).fixationFile, matchedFiles(i).ignoreFile);
    participantIDs{end+1} = extractBefore(matchedFiles(i).memoriseFile, '_memorise');

    % Load memorise and fixation files
    EEG_memorise = pop_loadset('filename', matchedFiles(i).memoriseFile, 'filepath', memoriseEpochFolder);
    EEG_fixation = pop_loadset('filename', matchedFiles(i).fixationFile, 'filepath', fixationEpochFolder);
    EEG_ignore = pop_loadset('filename', matchedFiles(i).ignoreFile, 'filepath', ignoreEpochFolder);

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);
    EEG_fixation = ensureICAActivations(EEG_fixation);
    EEG_ignore = ensureICAActivations(EEG_ignore);

    %% Step A1: Identify fmα Components Using Alpha Peak and DIPFIT (on "Ignore")
    fprintf('Identifying fmα components via PSD and DIPFIT on ignore epochs...\n');

    EEG_ignore = ensureICAActivations(EEG_ignore);
    EEG_ignore = performDipfit(EEG_ignore);

    num_ICs = size(EEG_ignore.icaact, 1);
    alpha_scores = zeros(1, num_ICs);

    for ic = 1:num_ICs
        data_ic = squeeze(EEG_ignore.icaact(ic, :, :));
        data_ic = data_ic(:);
        [pxx, f] = pwelch(data_ic, [], [], [], EEG_ignore.srate);
        alpha_band = f >= alphaBand(1) & f <= alphaBand(2);
        alpha_scores(ic) = mean(pxx(alpha_band));
    end

    [valid_idx, distances] = validate_dipoles(EEG_ignore, rv_threshold, mni_constraints, ppc_centroid);
    fm_alpha_candidates = find(valid_idx);

    if isempty(fm_alpha_candidates)
        [~, fallback_idx] = min(distances);
        fm_alpha_idx = fallback_idx;
        fprintf('No valid PPC dipoles. Using fallback IC %d.\n', fm_alpha_idx);
    else
        [~, best_idx] = max(alpha_scores(fm_alpha_candidates));
        fm_alpha_idx = fm_alpha_candidates(best_idx);
        fprintf('Selected fmα IC: %d\n', fm_alpha_idx);
    end

    %% Step A2: Compute ERSPs for fmα (on "memorise" data) and Apply PCA + ICA
    fprintf('Computing ERSPs for fmα component in memorise epochs...\n');

    ERSP_tensor = condition_ERSP_baseline_norm(EEG_fixation, EEG_memorise, fm_alpha_idx);

    [num_freqs, num_times, num_trials] = size(ERSP_tensor);
    ERSP_matrix = reshape(permute(ERSP_tensor, [3, 1, 2]), num_trials, []);

    [~, pca_scores, latent] = pca(ERSP_matrix);
    explained = cumsum(latent) / sum(latent) * 100;
    num_components = find(explained >= 85, 1);
    pca_scores_reduced = pca_scores(:, 1:num_components)';

    [weights, sphere] = runica(pca_scores_reduced, 'extended', 1, 'stop', 1e-7);
    ica_templates = pinv(weights * sphere)';
    templates_reshaped = reshape(ica_templates, num_components, num_freqs, num_times);

    fprintf('Extracted %d ICA templates.\n', num_components);


    %% Step A3: Cluster Components Across Participants
    % (This part will be implemented after processing all participants.)
    fprintf('Clustering components will be performed after all participants are processed.\n');

    %% Step B1: Extract Component Activity
    fprintf('Extracting component activity for fmα...\n');
    % Ensure ICA activations
    EEG_ignore = ensureICAActivations(EEG_ignore);

    % Extract the time series activity for the identified fmα component
    fm_alpha_activity = EEG_ignore.icaact(fm_alpha_idx, :, :); % Dimensions: [1 x time points x epochs]

    % Reshape the activity into [time points x epochs]
    fm_alpha_activity = squeeze(fm_alpha_activity); % Dimensions: [time points x epochs]memoryLoads
    fprintf('fmα component activity extracted with dimensions: [%d x %d]\n', size(fm_alpha_activity, 1), size(fm_alpha_activity, 2));

    %% Step B2: Separate Epochs by Memory Load (memorise)

    fprintf('Separating memorise epochs by memory load...\n');
% Sort memorise epochs using eventAccLoad
memorise_by_load = cell(1, length(accumulatedLoads));
for e = 1:length(EEG_memorise.epoch)
    if isfield(EEG_memorise.epoch(e), 'eventAccLoad')
        accLoadStr = EEG_memorise.epoch(e).eventAccLoad{1};
        loadNum = sscanf(accLoadStr, 'L%d');
        if ~isempty(loadNum) && loadNum <= 5
            memorise_by_load{loadNum+1} = [memorise_by_load{loadNum+1}, e];
        else
            warning('Unexpected AccLoad format in memorise epoch %d: %s', e, accLoadStr);
        end
    else
        warning('Missing eventAccLoad in memorise epoch %d', e);
    end
end
% Preallocate containers for accumulated memory load subsets (L0 to L5)
accumulatedLoads = 0:5;
memorise_erspMeans = cell(1, length(accumulatedLoads));
ignore_erspMeans = cell(1, length(accumulatedLoads));

    fprintf('Assigning fixation epochs cyclically by memory load...\n');
    
    %% Preallocate containers for fixation load subsets
    fixation_by_load = cell(1, length(memoryLoads));
    for e = 1:length(EEG_fixation.epoch)
        load_idx = mod(e - 1, length(memoryLoads)) + 1; % Cyclic assignment
        fixation_by_load{load_idx} = [fixation_by_load{load_idx}, e];
    end
    fprintf('Fixation epochs assigned cyclically to memory load conditions.\n');

    %% Step B2: Separate Epochs by Memory Load (ignore)

    fprintf('Separating ignore epochs by memory load...\n');
% Sort ignore epochs using eventAccLoad
ignore_by_load = cell(1, length(accumulatedLoads));
for e = 1:length(EEG_ignore.epoch)
    if isfield(EEG_ignore.epoch(e), 'eventAccLoad')
        accLoadStr = EEG_ignore.epoch(e).eventAccLoad{1};
        loadNum = sscanf(accLoadStr, 'L%d');
        if ~isempty(loadNum) && loadNum <= 5
            ignore_by_load{loadNum+1} = [ignore_by_load{loadNum+1}, e];
        else
            warning('Unexpected AccLoad format in ignore epoch %d: %s', e, accLoadStr);
        end
    else
        warning('Missing eventAccLoad in ignore epoch %d', e);
    end
end
% Preallocate containers for accumulated memory load subsets (L0 to L5)
accumulatedLoads = 0:5;
memorise_erspMeans = cell(1, length(accumulatedLoads));
ignore_erspMeans = cell(1, length(accumulatedLoads));

    fprintf('Assigning fixation epochs cyclically by memory load...\n');
    
 
    %% Step B3: Calculate alpha Power (Memorise)
    fprintf('Calculating alpha power for memorise...\n');
    
    % Function to compute alpha power for a set of epochs
    compute_power = @(activity, epochs, alpha_band, srate) arrayfun(@(e) ...
        mean(mean(abs(fft(activity(:, e), [], 1)).^2, 1)), epochs);

    % Calculate alpha power for memorise epochs
    alpha_power_memorise = cellfun(@(epochs) ...
    computealphaPower(fm_alpha_activity, epochs, alphaBand, EEG_memorise.srate), ...
    memorise_by_load, 'UniformOutput', false);

    % Calculate alpha power for fixation epochs
    alpha_power_fixation = cellfun(@(epochs) ...
        compute_power(fm_alpha_activity, epochs, alphaBand, EEG_fixation.srate), ...
        fixation_by_load, 'UniformOutput', false);

    % Aggregate alpha power
    fprintf('Aggregating and baseline-normalising alpha power...\n');
    % Preallocate for mean alpha power
    mean_alpha_power = zeros(1, length(memoryLoads));
    
    % Loop through memory load conditions
    for i = 1:length(memoryLoads)
        % Aggregate mean alpha power for memorise and fixation epochs
        mean_memorise_power = mean(alpha_power_memorise{i});
        mean_fixation_power = mean(alpha_power_fixation{i});
    
        % Baseline-normalise memorise power
        mean_alpha_power(i) = mean_memorise_power - mean_fixation_power;
    
        % Debug output
        fprintf('Load %d: MEMORISE Mean memorise alpha power = %.4f, Mean fixation alpha power = %.4f, Normalised = %.4f\n', ...
            memoryLoads(i), mean_memorise_power, mean_fixation_power, mean_alpha_power(i));
    end

    % Store dipole information
    fm_alpha_idx_di = [fm_alpha_idx_di; fm_alpha_idx];
    fm_alpha_rv_di = [fm_alpha_rv_di; EEG_memorise.dipfit.model(fm_alpha_idx).rv;];
    posX = [posX, EEG_memorise.dipfit.model(fm_alpha_idx).posxyz(1)]; % Extract X-coordinates
    posY = [posY, EEG_memorise.dipfit.model(fm_alpha_idx).posxyz(2)]; % Extract Y-coordinates
    posZ = [posZ, EEG_memorise.dipfit.model(fm_alpha_idx).posxyz(3)]; % Extract Z-coordinates
        

    %% Step B3: Calculate alpha Power (ignore)
    fprintf('Calculating alpha power for ignore...\n');
    

    % Calculate alpha power for memorise epochs
    alpha_power_ignore = cellfun(@(epochs) ...
    computealphaPower(fm_alpha_activity, epochs, alphaBand, EEG_ignore.srate), ...
    ignore_by_load, 'UniformOutput', false);


    % Aggregate alpha power
    fprintf('Aggregating and baseline-normalising alpha power...\n');
    % Preallocate for mean alpha power
    ignore_mean_alpha_power = zeros(1, length(ignoreMemoryLoads));
    
    % Loop through memory load conditions
    for i = 1:length(ignoreMemoryLoads)
        % Aggregate mean alpha power for memorise and fixation epochs
        mean_ignore_power = mean(alpha_power_ignore{i});
        mean_fixation_power = mean(alpha_power_fixation{i});
    
        % Baseline-normalise memorise power
        ignore_mean_alpha_power(i) = mean_ignore_power - mean_fixation_power;
    
        % Debug output
        fprintf('Load %d: IGNORE Mean memorise alpha power = %.4f, Mean fixation alpha power = %.4f, Normalised = %.4f\n', ...
            ignoreMemoryLoads(i), mean_ignore_power, mean_fixation_power, ignore_mean_alpha_power(i));
    end


    %% Save Participant-Level Data
% Store participant-level data dynamically
for load = accumulatedLoads
    memorise_erspMeans{load+1} = [memorise_erspMeans{load+1}; mean_alpha_power(load+1)];
    ignore_erspMeans{load+1} = [ignore_erspMeans{load+1}; ignore_mean_alpha_power(load+1)];
end

end

%% Export Results to CSV
% Create a table for CSV export
csvTable = table(participantIDs');
for load = accumulatedLoads
    csvTable.(['Memorise_ERSP_L' num2str(load)]) = memorise_erspMeans{load+1};
    csvTable.(['Ignore_ERSP_L' num2str(load)]) = ignore_erspMeans{load+1};
end

disp(['CSV file saved: ' csvFileName]);

% Print sizes to console
fprintf('Size of fm_alpha_idx_di: [%d x %d]\n', size(fm_alpha_idx_di, 1), size(fm_alpha_idx_di, 2));
fprintf('Size of fm_alpha_rv_di: [%d x %d]\n', size(fm_alpha_rv_di, 1), size(fm_alpha_rv_di, 2));
fprintf('Size of posX: [%d x %d]\n', size(posX, 1), size(posX, 2));
fprintf('Size of posY: [%d x %d]\n', size(posY, 1), size(posY, 2));
fprintf('Size of posZ: [%d x %d]\n', size(posZ, 1), size(posZ, 2));


% Create the table
csvTableDipole = table(participantIDs', fm_alpha_idx_di', posX, posY, posZ, fm_alpha_rv_di', ...
    'VariableNames', {'ParticipantID', 'fm_alpha_idx_di', 'PosX', 'PosY', 'PosZ', 'fm_alpha_rv_di'});
% Save the table as a CSV file
csvFileName_dt = fullfile(outputFolder, 'dipole_info.csv'); % Define file path
writetable(csvTableDipole, csvFileName_dt);
% Display success message
disp(['Dipole information saved as CSV: ' csvFileName_dt]);



%% Step C1: Statistical Analysis
% Create table dynamically for ANOVA
memorise_data = cell2mat(cellfun(@(x) mean(x), memorise_erspMeans, 'UniformOutput', false));
ignore_data = cell2mat(cellfun(@(x) mean(x), ignore_erspMeans, 'UniformOutput', false));
all_data = [memorise_data, ignore_data];

loadNames = [arrayfun(@(x) sprintf('Memorise_Load%d', x), accumulatedLoads, 'UniformOutput', false), ...
             arrayfun(@(x) sprintf('Ignore_Load%d', x), accumulatedLoads, 'UniformOutput', false)];

dataTable = array2table(all_data, 'VariableNames', loadNames);
dataTable.ParticipantID = (1:size(all_data, 1))';

% Dynamic within-subject factors
condition = categorical([repmat({'Memorise'}, 1, length(accumulatedLoads)), ...
                        repmat({'Ignore'}, 1, length(accumulatedLoads))])';
memoryLoad = categorical(repmat(arrayfun(@(x) sprintf('Load%d', x), accumulatedLoads, 'UniformOutput', false), 1, 2))';
withinDesign = table(condition, memoryLoad, 'VariableNames', {'Condition', 'MemoryLoad'});

% Fit repeated measures model and run ANOVA
rm = fitrm(dataTable, sprintf('%s~1', strjoin(loadNames, '-')), 'WithinDesign', withinDesign);
ranovaResults = ranova(rm);
disp('Repeated Measures ANOVA Results:');
disp(ranovaResults);

fprintf('Visualising results...\n');
fprintf('Visualising results dynamically for all loads...\n');
memoryLoads = 0:5;

% Mean and standard error already computed in previous step

figure;
errorbar(memoryLoads, mean_memorise, std_memorise, 'o-', 'DisplayName', 'Memorise');
hold on;
errorbar(memoryLoads, mean_ignore, std_ignore, 'o-', 'DisplayName', 'Ignore');
xlabel('Memory Load (Letters)');
ylabel('Baseline-Normalised Alpha Power');
title('Alpha Power Across Memory Loads by Condition');
legend;
saveas(gcf, fullfile(outputFolder, 'alpha_power_memory_loads_conditions.png'));



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

function alphaPower = computealphaPower(activity_data, epochs, alphaBand, srate)
    % computealphaPower Computes alpha power for given epochs and alpha band
    %
    % Inputs:
    %   activity_data - Time x Epochs matrix (e.g., 1500 x 609)
    %   epochs - Epoch indices to compute power for
    %   alphaBand - [min_alpha max_alpha], frequency range for alpha
    %   srate - Sampling rate in Hz
    %
    % Output:
    %   alphaPower - Mean alpha power across epochs and time
    
    % Extract epochs from activity data
    selected_data = activity_data(:, epochs); % Dimensions: Time x Selected Epochs
    
    % Compute power spectral density (PSD) using FFT
    nfft = size(selected_data, 1); % Number of FFT points (equal to time samples)
    freqs = linspace(0, srate / 2, floor(nfft / 2) + 1); % Frequency vector
    alpha_range = find(freqs >= alphaBand(1) & freqs <= alphaBand(2)); % Indices for alpha band
    
    % Calculate FFT for each epoch
    alphaPower_epochs = zeros(1, size(selected_data, 2)); % Preallocate alpha power
    for epoch = 1:size(selected_data, 2)
        % FFT for single epoch
        fft_result = fft(selected_data(:, epoch), nfft);
        psd = abs(fft_result(1:floor(nfft / 2) + 1)).^2; % One-sided PSD
        
        % Compute mean alpha power
        alphaPower_epochs(epoch) = mean(psd(alpha_range));
    end
    
    % Compute mean alpha power across all selected epochs
    alphaPower = mean(alphaPower_epochs);
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



function ERSP_tensor = condition_ERSP_baseline_norm(EEG_fixation, EEG_memorise, fm_alpha_idx)
    % Compute ERSPs for selected fmα component across trials
    
        activity_data = squeeze(EEG_memorise.icaact(fm_alpha_idx, :, :));
        baseline_data = squeeze(EEG_fixation.icaact(fm_alpha_idx, :, :));
    
        baseline_flat = baseline_data(:);
        baseline_mean = mean(abs(baseline_flat).^2);
    
        num_trials = size(activity_data, 2);
    
        [test_ersp, ~, ~, times, freqs] = newtimef( ...
            activity_data(:, 1), EEG_memorise.pnts, ...
            [EEG_memorise.xmin, EEG_memorise.xmax]*1000, ...
            EEG_memorise.srate, [3], ...
            'baseline', baseline_mean, ...
            'plotersp', 'off', 'plotitc', 'off');
    
        num_freqs = size(test_ersp, 1);
        num_times = size(test_ersp, 2);
        ERSP_tensor = zeros(num_freqs, num_times, num_trials);
    
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
    


function alpha_template_idx = find_alpha_template(ica_templates, freqs, alphaBand, visualise)
    % find_alpha_template Identifies the dominant alpha template (8-12 Hz) from ICA templates
    %
    % Inputs:
    %   ica_templates - Independent components (ICs), size: [num_ICs x features]
    %   freqs - [min_freq max_freq], frequency range of the decomposition
    %   alphaBand - [min_alpha max_alpha], frequency range of the alpha band
    %   visualise - Boolean, whether to visualise the dominant alpha template
    %
    % Output:
    %   alpha_template_idx - Index of the IC with dominant alpha activity (8-12 Hz)
    
    % Compute the dimensions of the frequency and time axes
    num_freqs = length(freqs); % Use provided frequency vector length
    num_times = size(ica_templates, 2) / num_freqs; % Derive time dimension

    % Check validity of reshaping
    if round(num_times) ~= num_times
        error('Frequency-time features do not align with the frequency vector length.');
    end

    % Compute the frequency vector
    freq_vector = linspace(freqs(1), freqs(2), num_freqs);
    
    % Identify indices corresponding to the alpha band
    alpha_range = find(freq_vector >= alphaBand(1) & freq_vector <= alphaBand(2));
    
    % Initialise alpha scores for each IC
    num_ICs = size(ica_templates, 1);
    alpha_scores = zeros(1, num_ICs);

    % Loop through each IC
    for ic = 1:num_ICs
        % Reshape IC template back to frequency-time representation
        ic_template = reshape(ica_templates(ic, :), num_freqs, num_times);
        
        % Extract power in the alpha band
        alpha_power = mean(ic_template(alpha_range, :), 1); % Mean over alpha frequencies
        alpha_scores(ic) = mean(alpha_power); % Average over time
    end

    % Identify the IC with the maximum alpha score
    [~, alpha_template_idx] = max(alpha_scores);

    % Visualise the dominant alpha template if requested
    if visualise
        dominant_alpha_template = reshape(ica_templates(alpha_template_idx, :), num_freqs, num_freqs);
        times = linspace(0, 1, num_freqs); % Example time vector (adjust as needed)
        
        figure;
        imagesc(times, freq_vector, dominant_alpha_template);
        axis xy;
        colorbar;
        title('Dominant alpha Template (5–7 Hz)');
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
    end
end



function alpha_scores = score_alpha_alignment(ica_templates, alpha_template_idx, freqs, alphaBand)
    % score_alpha_alignment Scores ICs based on alignment with the alpha template
    %
    % Inputs:
    %   ica_templates - IC templates (num_ICs x features)
    %   alpha_template_idx - Index of the dominant alpha template
    %   freqs - Frequency vector (based on original ERSP dimensions)
    %   alphaBand - [min_alpha max_alpha], frequency range of the alpha band

    % Recreate frequency vector based on the original ERSP dimensions
    num_features = size(ica_templates, 2);
    freq_vector = linspace(freqs(1), freqs(2), num_features);

    % Identify indices corresponding to the alpha band
    alpha_indices = find(freq_vector >= alphaBand(1) & freq_vector <= alphaBand(2));
    fprintf('Debug: alpha band indices: %s\n', mat2str(alpha_indices));

    % Extract and normalise alpha template
    alpha_template = ica_templates(alpha_template_idx, :);
    alpha_template_norm = alpha_template / norm(alpha_template);

    % Calculate alpha scores
    num_ICs = size(ica_templates, 1);
    alpha_scores = zeros(1, num_ICs);
    for ic = 1:num_ICs
        % Extract IC template
        ic_template = ica_templates(ic, :);
        ic_template_norm = ic_template / norm(ic_template);

        % Compute alpha alignment as dot product (cosine similarity)
        alpha_scores(ic) = dot(alpha_template_norm(alpha_indices), ic_template_norm(alpha_indices));
    end
end


function [valid_idx, distances] = validate_dipoles(EEG, rv_threshold, mni_constraints, ppc_centroid)
    % Validate dipoles based on residual variance and MNI constraints, and calculate distances from PPC centroid
    % Inputs:
    %   EEG - EEG structure with dipfit model
    %   rv_threshold - Residual variance threshold (e.g., 0.15)
    %   mni_constraints - [X_min X_max; Y_min Y_max; Z_min Z_max] matrix defining MNI constraints
    %   ppc_centroid - [x, y, z], PPC centroid coordinates (e.g., [-5, 20, 25])
    % Outputs:
    %   valid_idx - Logical array indicating valid ICs (1 for valid, 0 for invalid)
    %   distances - Array of distances from PPC centroid for all ICs

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

        % Compute Euclidean distance from PPC centroid
        posxyz = dipole.posxyz;
        distances(ic) = norm(posxyz - ppc_centroid);

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
    fprintf('Debug: Distances from PPC centroid: %s\n', mat2str(distances));
    fprintf('Debug: Valid indices: %s\n', mat2str(find(valid_idx)));

    % Handle mismatched IC and dipole counts
    if length(valid_idx) ~= length(distances)
        warning('Mismatch: Number of ICs (%d) does not match number of dipoles (%d). Adjusting to smaller count.', ...
                length(valid_idx), length(distances));
        valid_idx = valid_idx(1:min(length(valid_idx), length(distances))); % Adjust valid_idx length
        distances = distances(1:min(length(valid_idx), length(distances))); % Adjust distances length
    end
end
