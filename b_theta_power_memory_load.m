
% Theta Power and Memory Load Analysis

%{
Input: Memorise epochs, fixation epochs
Output: Theta Power vs Memory Load results
Summary: Identify fmθ Component, Validate Near ACC, Compute ERSPs, and Perform ANOVA
%}

% Set variables
clear;

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

% Preallocate variables for CSV and ANOVA
csvData = [];
participantIDs = {};
erspL3_means = [];
erspL5_means = [];
erspL7_means = [];


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
    participantIDs{end+1} = matchedFiles(i).memoriseFile;

    % Load memorise and fixation files
    EEG_memorise = pop_loadset('filename', matchedFiles(i).memoriseFile, 'filepath', memoriseEpochFolder);
    EEG_fixation = pop_loadset('filename', matchedFiles(i).fixationFile, 'filepath', fixationEpochFolder);

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);
    EEG_fixation = ensureICAActivations(EEG_fixation);

    if isempty(EEG_memorise.icaact)
        disp("Something has gone wrong")
    end

    %% Dipole Localization: Identify All ICs Near ACC
    EEG_memorise = performDipfit(EEG_memorise);

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);
    EEG_fixation = ensureICAActivations(EEG_fixation);

    if isempty(EEG_memorise.icaact)
        disp("Something has gone very wrong")
    end

    candidateICs = findComponentsNearACC(EEG_memorise);

    if isempty(candidateICs)  
        if ~isfield(EEG_memorise.dipfit, 'model') || isempty(EEG_memorise.dipfit.model) || ~isstruct(EEG_memorise.dipfit.model)
            warning('DIPFIT model is empty or not initialized. Skipping participant.');
            continue;
        end
        
        % Compute distances of all ICs to the ACC (MNI [0, 20, 40])
        accCoords = [0, 20, 40];
        numICs = length(EEG_memorise.dipfit.model); % Safe to access after the above check
        distances = nan(numICs, 1); % Preallocate with NaN

        % Loop through all ICs to compute distances and theta power
        for ic = 1:numICs

            dipole = EEG_memorise.dipfit.model(ic);

            % Skip invalid dipoles
            if isempty(dipole.posxyz) || isempty(dipole.rv)
                distances(ic) = Inf;
                continue;
            end
    
            % Compute Euclidean distance to ACC
            distances(ic) = norm(dipole.posxyz - accCoords);
    
            % Safeguard against invalid indexing
            if ic > size(EEG_memorise.icaact, 1)
                warning('Skipping IC%d: Index exceeds available ICs in icaact.', ic);
                distances(ic) = Inf;
                continue;
            end
    
        end

        % Find the 10 closest ICs to the ACC
        [sortedDistances, sortedIndices] = sort(distances);
        top10ICs = sortedIndices(1:min(10, length(sortedIndices))); % Get up to 10 closest ICs
        fprintf('10 Closest ICs to ACC:\n');
        for i = 1:length(top10ICs)
            icIdx = top10ICs(i);
            fprintf('IC%d: Distance = %.2f mm\n', ...
                    icIdx, sortedDistances(i));
        end

        % Assign the closest IC as the fallback
        candidateICs = top10ICs;
        disp("The top 10 closest ICs are now considered the candiate ICs")

        % Compute theta power for the candidate ICs
        thetaPower = computeThetaPower(EEG_memorise, EEG_fixation, thetaBand, candidateICs);
        
        fprintf('Number of candidate ICs: %d\n', length(candidateICs));
        fprintf('Length of thetaPower array: %d\n', length(thetaPower));


        % Add theta power to the candidateICs for plotting
        for i = 1:length(candidateICs)
            icIdx = candidateICs(i);
            fprintf('IC%d: Distance = %.2f mm, Theta Power = %.2f\n', ...
                    icIdx, sortedDistances(i), thetaPower(i));
        end

        % Plot scalp topographies of the 10 closest ICs
        figure;
        for i = 1:length(candidateICs)
            icIdx = candidateICs(i);
            subplot(2, 5, i); % Arrange in a 2x5 grid
            pop_topoplot(EEG_memorise, 0, icIdx, sprintf('IC%d', icIdx), 0, 'electrodes', 'on');
            title(sprintf('IC%d\nDist: %.2f mm\nTheta: %.2f', icIdx, sortedDistances(i), thetaPower(i)));
        end
        sgtitle('Top 10 Closest ICs to ACC');
    
        % Plot the ACC location and top 10 ICs in MNI space
        figure;
        scatter3(accCoords(1), accCoords(2), accCoords(3), 100, 'r', 'filled'); hold on;
        for i = 1:length(candidateICs)
            icIdx = candidateICs(i);
            scatter3(EEG_memorise.dipfit.model(icIdx).posxyz(1), ...
                     EEG_memorise.dipfit.model(icIdx).posxyz(2), ...
                     EEG_memorise.dipfit.model(icIdx).posxyz(3), 100, 'b', 'filled');
            text(EEG_memorise.dipfit.model(icIdx).posxyz(1), ...
                 EEG_memorise.dipfit.model(icIdx).posxyz(2), ...
                 EEG_memorise.dipfit.model(icIdx).posxyz(3), sprintf('IC%d', icIdx));
        end
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
        legend({'ACC', 'Top 10 ICs'}, 'Location', 'best');
        grid on; view(3);
        title('ACC and Top 10 ICs in MNI Space');
    
    end

       
    %% Compute Theta Power and Select Best fmθ Component
    thetaPower = computeThetaPower(EEG_memorise, EEG_fixation, thetaBand, candidateICs);
    [~, fmThetaIC] = max(thetaPower);
    fmThetaIC = candidateICs(fmThetaIC); % Select IC with max theta power among candidates
    fprintf('Selected fmθ Component: IC%d\n', fmThetaIC);

    %% Visual Inspection of Scalp Topography
    figure;
    pop_topoplot(EEG_memorise, 0, fmThetaIC, sprintf('Topography of IC%d (Theta Candidate)', fmThetaIC), 0, 'electrodes', 'on');
    title(sprintf('Inspect IC%d: Does it align with frontal midline activity?', fmThetaIC));
    disp('Visually inspect the scalp topography for frontal midline activity.');

    %% Compute ERSP for Each Memory Load

    % Define memory load codes
    memoryLoadCodes = {'s3', 's5', 's7'}; % Corresponding to loads 3, 5, and 7
    
    % Initialize structures to hold data
    fmThetaData = struct('L3', [], 'L5', [], 'L7', []);
    baselineERSP = struct('L3', [], 'L5', [], 'L7', []);
    
    % Extract trial events and ensure compatibility
    trialEvents = {EEG_memorise.epoch.eventtype};
    if iscell(trialEvents)
        trialEvents = cellfun(@(x) x{1}, trialEvents, 'UniformOutput', false);
    elseif isstring(trialEvents)
        trialEvents = cellstr(trialEvents);
    elseif isnumeric(trialEvents)
        trialEvents = arrayfun(@num2str, trialEvents, 'UniformOutput', false);
    else
        error('Unsupported data type in trialEvents.');
    end
    
    % Separate fmThetaData and compute baseline ERSP for each memory load
    for loadIdx = 1:length(memoryLoads)
        loadCode = memoryLoadCodes{loadIdx};
        loadTrials = contains(trialEvents, loadCode);
        
        if any(loadTrials)
            % Extract task data for current memory load
            fmThetaData.(sprintf('L%d', memoryLoads(loadIdx))) = EEG_memorise.icaact(fmThetaIC, :, loadTrials);
            
            % Compute baseline ERSP using all fixation trials
            baselineData = reshape(EEG_fixation.icaact(fmThetaIC, :, :), 1, []); % Combine all fixation data
            [baselineERSP.(sprintf('L%d', memoryLoads(loadIdx))), ~, ~, ~, ~] = ...
                newtimef(baselineData, EEG_fixation.pnts, ...
                         [EEG_fixation.xmin EEG_fixation.xmax]*1000, EEG_fixation.srate, [3 0.5], ...
                         'baseline', NaN, 'plotitc', 'off', 'plotersp', 'off', ...
                         'freqs', thetaBand, 'nfreqs', 10, 'freqscale', 'linear', 'padratio', 2);
        else
            warning('No trials found for memory load %d.', memoryLoads(loadIdx));
        end
    end
    
    % Function to compute ERSP for task data
    computeERSP = @(data) newtimef(data, EEG_memorise.pnts, ...
        [EEG_memorise.xmin EEG_memorise.xmax]*1000, EEG_memorise.srate, [3 0.5], ...
        'baseline', NaN, 'trialbase', 'full', 'freqs', thetaBand, 'nfreqs', 10, ...
        'freqscale', 'linear', 'basenorm', 'off', 'plotitc', 'off', 'plotersp', 'off', 'padratio', 2);
    
    % Compute ERSPs for each memory load
    [erspL3, ~, ~, timesL3, freqsL3] = computeERSP(fmThetaData.L3);
    %[erspL5, ~, ~, timesL5, freqsL5] = computeERSP(fmThetaData.L5);
    %[erspL7, ~, ~, timesL7, freqsL7] = computeERSP(fmThetaData.L7);
    
    % Normalize task ERSP by baseline ERSP
    erspL3 = erspL3 - baselineERSP.L3;
    %erspL5 = erspL5 - baselineERSP.L5;
    %erspL7 = erspL7 - baselineERSP.L7;
    
    %{
    % Store results in a structure for further analysis
    erspResults = struct('L3', struct('ersp', erspL3, 'times', timesL3, 'freqs', freqsL3), ...
                         'L5', struct('ersp', erspL5, 'times', timesL5, 'freqs', freqsL5), ...
                         'L7', struct('ersp', erspL7, 'times', timesL7, 'freqs', freqsL7));
    %}

    % Store results in a structure for further analysis
    erspResults = struct('L3', struct('ersp', erspL3, 'times', timesL3, 'freqs', freqsL3));

    % Save ERSP results
    save(fullfile(outputFolder, 'ersp_results.mat'), 'erspResults');
    disp('ERSP results saved.');

    % Define theta frequency range
    thetaRange = [4 7]; % in Hz

    % Find indices corresponding to theta frequencies
    thetaIndicesL3 = freqsL3 >= thetaRange(1) & freqsL3 <= thetaRange(2);
    %thetaIndicesL5 = freqsL5 >= thetaRange(1) & freqsL5 <= thetaRange(2);
    %thetaIndicesL7 = freqsL7 >= thetaRange(1) & freqsL7 <= thetaRange(2);

    % Compute mean ERSP within the theta band for each load
    meanERSP_L3 = mean(mean(erspL3(thetaIndicesL3, :), 1), 2);
    %meanERSP_L5 = mean(mean(erspL5(thetaIndicesL5, :), 1), 2);
    %meanERSP_L7 = mean(mean(erspL7(thetaIndicesL7, :), 1), 2);
    
    % Store the mean ERSP values
    erspL3_means = [erspL3_means; meanERSP_L3];
    %erspL5_means = [erspL5_means; meanERSP_L5];
    %erspL7_means = [erspL7_means; meanERSP_L7];


    %% Append data for CSV and ANOVA
    % Extract participant ID without file extension
    erspL3_means(end+1) = erspL3_means;
    %erspL5_means(end+1) = meanL5;
    %erspL7_means(end+1) = meanL7;
end

%% Export Results to CSV
% Create a table for CSV export
%{
csvTable = table(participantIDs', erspL3_means', erspL5_means', erspL7_means', ...
    'VariableNames', {'ParticipantID', 'ERSP_L3', 'ERSP_L5', 'ERSP_L7'});
%}
csvTable = table(participantIDs', erspL3_means', ...
    'VariableNames', {'ParticipantID', 'ERSP_L3'});
% Save the csv 
csvFileName = fullfile(outputFolder, 'theta_power_memory_load.csv');
writetable(csvTable, csvFileName);
disp(['CSV file saved: ' csvFileName]);

%% Perform Repeated-Measures ANOVA
loadLabels = arrayfun(@(x) sprintf('Load_%d', x), memoryLoads, 'UniformOutput', false);
thetaPowerTable = array2table(thetaPowerMatrix, 'VariableNames', loadLabels);

% ANOVA Model
rm = fitrm(thetaPowerTable, sprintf('%s-%s~1', loadLabels{1}, loadLabels{end}), 'WithinDesign', memoryLoads');
ranovaResults = ranova(rm);

% Display Results
disp('Repeated-Measures ANOVA Results:');
disp(ranovaResults);

%% Save ANOVA Results
save(fullfile(outputFolder, 'theta_power_memory_results.mat'), 'thetaPowerMatrix', 'ranovaResults');
disp('ANOVA results saved.');


%% Utility Functions
function EEG = ensureICAActivations(EEG)
    % Compute ICA activations if missing
    if isempty(EEG.icaact)
        reshapedData = reshape(EEG.data, size(EEG.data, 1), []);
        icaActTemp = EEG.icaweights * EEG.icasphere * reshapedData;
        EEG.icaact = reshape(icaActTemp, size(EEG.icaweights, 1), EEG.pnts, EEG.trials);
        disp('ICA activations computed.');
    end
end

function thetaPower = computeThetaPower(EEG_task, EEG_fixation, thetaBand, candidateICs)
    % Computes baseline-normalised theta power in dB
    thetaPower = zeros(length(candidateICs), 1);
    for idx = 1:length(candidateICs)
        ic = candidateICs(idx);

        % Task data
        icData_task = reshape(EEG_task.icaact(ic, :, :), 1, []);
        fftData_task = fft(icData_task, 512); % Zero-padded FFT
        powerSpectrum_task = abs(fftData_task(1:256)).^2; % Single-sided
        
        % Fixation (baseline) data
        icData_fixation = reshape(EEG_fixation.icaact(ic, :, :), 1, []);
        fftData_fixation = fft(icData_fixation, 512); % Zero-padded FFT
        powerSpectrum_fixation = abs(fftData_fixation(1:256)).^2; % Single-sided
        
        % Frequency vector
        freqRes = EEG_task.srate / 512;
        f = (0:255) * freqRes;
        
        % Theta band indices
        thetaIdx = f >= thetaBand(1) & f <= thetaBand(2);
        
        % Mean power in theta band for task and fixation
        meanThetaPower_task = mean(powerSpectrum_task(thetaIdx));
        meanThetaPower_fixation = mean(powerSpectrum_fixation(thetaIdx));
        
        % Convert to decibels (10 * log10)
        taskPower_dB = 10 * log10(meanThetaPower_task);
        fixationPower_dB = 10 * log10(meanThetaPower_fixation);
        
        % Baseline-normalized theta power (difference in dB)
        thetaPower(idx) = taskPower_dB - fixationPower_dB;
    end
end


%{
function EEG = performDipfit(EEG)
    EEG = pop_dipfit_settings(EEG, 'hdmfile', 'standard_BEM/standard_vol.mat', ...
        'mrifile', 'standard_BEM/standard_mri.mat', 'chanfile', ...
        'standard_BEM/elec/standard_1005.elc', 'coordformat', 'MNI');
    EEG = pop_multifit(EEG, 1:size(EEG.icaweights, 1), 'threshold', 100, 'dipplot', 'off');
end
%}

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
    evalc('EEG = pop_multifit(EEG, 1:size(EEG.icaweights, 1), ''threshold'', 100, ''dipplot'', ''off'');');
end


function candidateICs = findComponentsNearACC(EEG) % Finds dipoles in ACC using MNI Coordinate range X: -10 to +10 (midline), Y: 10 to 50 (anterior-posterior), Z: 10 to 50 (superior-inferior)
    candidateICs = [];
    for ic = 1:length(EEG.dipfit.model)
        dipole = EEG.dipfit.model(ic);
        % Check for residual variance and dACC location
        if dipole.rv < 0.2 && abs(dipole.posxyz(1)) <= 10 ... % Midline constraint
                && dipole.posxyz(2) >= 10 && dipole.posxyz(2) <= 50 ... % Anterior-posterior constraint
                && dipole.posxyz(3) >= 10 && dipole.posxyz(3) <= 50 % inferior - superior constraint
            candidateICs = [candidateICs, ic];
        end
    end
end