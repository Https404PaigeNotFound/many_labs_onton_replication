
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
timeWindow = [0 1.5]; % Task-relevant time window for analysis
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
thetaPowerMatrix = []; % Rows = Participants, Columns = Memory Loads

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

    % Load memorise and fixation files
    EEG_memorise = pop_loadset('filename', matchedFiles(i).memoriseFile, 'filepath', memoriseEpochFolder);
    EEG_fixation = pop_loadset('filename', matchedFiles(i).fixationFile, 'filepath', fixationEpochFolder);

    % Ensure ICA activations
    EEG_memorise = ensureICAActivations(EEG_memorise);
    EEG_fixation = ensureICAActivations(EEG_fixation);

    %% Dipole Localization: Identify All ICs Near ACC
    EEG_memorise = performDipfit(EEG_memorise);
    candidateICs = findComponentsNearACC(EEG_memorise);
    % If not ICs near the ACC
    %{
    if isempty(candidateICs)
        warning('No ICs near ACC found for participant. Skipping.');
        continue;
    end
    %}
    if isempty(candidateICs)  
        if ~isfield(EEG.dipfit, 'model') || isempty(EEG.dipfit.model) || ~isstruct(EEG.dipfit.model)
            warning('DIPFIT model is empty or not initialized. Skipping participant.');
            continue;
        end
        
        % Compute distances of all ICs to the ACC (MNI [0, 20, 40])
        accCoords = [0, 20, 40];
        numICs = length(EEG.dipfit.model); % Safe to access after the above check
        distances = nan(numICs, 1); % Preallocate with NaN


        for ic = 1:length(EEG.dipfit.model)
            dipole = EEG.dipfit.model(ic);
            if ~isempty(dipole.posxyz) && ~isempty(dipole.rv)
                distances(ic) = norm(dipole.posxyz - accCoords); % Euclidean distance
            else
                distances(ic) = Inf; % Handle empty or invalid positions
            end
        end
        
        % Find the IC closest to ACC
        [minDistance, nearestIC] = min(distances);
        fprintf('No ICs near ACC. Closest IC is IC%d with distance %.2f mm.\n', nearestIC, minDistance);
        fprintf('MNI Coordinates of closest IC: [%.2f, %.2f, %.2f]\n', ...
                EEG.dipfit.model(nearestIC).posxyz(1), ...
                EEG.dipfit.model(nearestIC).posxyz(2), ...
                EEG.dipfit.model(nearestIC).posxyz(3));
    
        % Plot the scalp topography of the nearest IC
        figure;
        subplot(1, 2, 1);
        pop_topoplot(EEG, 0, nearestIC, sprintf('Closest IC%d', nearestIC), 0, 'electrodes', 'on');
        title(sprintf('IC%d: Closest to ACC (Distance: %.2f mm)', nearestIC, minDistance));
    
        % Plot the ACC location as a reference
        subplot(1, 2, 2);
        scatter3(accCoords(1), accCoords(2), accCoords(3), 100, 'r', 'filled'); hold on;
        scatter3(EEG.dipfit.model(nearestIC).posxyz(1), ...
                 EEG.dipfit.model(nearestIC).posxyz(2), ...
                 EEG.dipfit.model(nearestIC).posxyz(3), 100, 'b', 'filled');
        xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
        legend({'ACC', sprintf('IC%d', nearestIC)}, 'Location', 'best');
        grid on; view(3);
        title('ACC and Closest IC in MNI Space');
        
        % Continue to next participant
        warning('No ICs near ACC. Nearest IC plotted for reference.');
        continue;
    end
    disp("IC near ACC found")

    %% Compute Theta Power and Select Best fmθ Component
    thetaPower = computeThetaPower(EEG_memorise, thetaBand, candidateICs);
    [~, fmThetaIC] = max(thetaPower);
    fmThetaIC = candidateICs(fmThetaIC); % Select IC with max theta power among candidates
    fprintf('Selected fmθ Component: IC%d\n', fmThetaIC);

    %% Visual Inspection of Scalp Topography
    figure;
    pop_topoplot(EEG_memorise, 0, fmThetaIC, sprintf('Topography of IC%d (Theta Candidate)', fmThetaIC), 0, 'electrodes', 'on');
    title(sprintf('Inspect IC%d: Does it align with frontal midline activity?', fmThetaIC));
    disp('Visually inspect the scalp topography for frontal midline activity.');

    %% Compute ERSP for Each Memory Load
    fmThetaData = EEG_memorise.icaact(fmThetaIC, :, :); 
    baselineData = EEG_fixation.icaact(fmThetaIC, :, :); 

    thetaByLoad = zeros(1, length(memoryLoads));
    trialEvents = {EEG_memorise.epoch.eventtype};

    for loadIdx = 1:length(memoryLoads)
        loadCode = ['s' num2str(memoryLoads(loadIdx))]; 
        loadTrials = contains(trialEvents, loadCode); 

        if any(loadTrials)
            fmThetaLoadData = fmThetaData(:, :, loadTrials);
            [ersp, ~, ~, times, freqs] = newtimef(fmThetaLoadData, EEG_memorise.pnts, ...
                [EEG_memorise.xmin EEG_memorise.xmax]*1000, EEG_memorise.srate, [3 0.5], ...
                'baseline', NaN, 'trialbase', 'full', 'freqs', thetaBand, 'nfreqs', 10, ...
                'freqscale', 'linear', 'baseline_data', baselineData, 'plotitc', 'off', 'plotersp', 'off', 'padratio', 2);
            timeIdx = times >= timeWindow(1)*1000 & times <= timeWindow(2)*1000;
            thetaByLoad(loadIdx) = mean(mean(ersp(freqs >= thetaBand(1) & freqs <= thetaBand(2), timeIdx), 1), 2);
        end
    end

    % Append data for CSV and ANOVA
    participantIDs{end+1} = matchedFiles(i).memoriseFile; % Participant ID
    csvData = [csvData; [thetaByLoad]]; % Append theta power
    thetaPowerMatrix = [thetaPowerMatrix; thetaByLoad];
end

%% Export Results to CSV
csvTable = array2table([csvData, participantIDs'], ...
    'VariableNames', {'ThetaPower_Load3', 'ThetaPower_Load5', 'ThetaPower_Load7', 'ParticipantID'});
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

function thetaPower = computeThetaPower(EEG, thetaBand, candidateICs)
    thetaPower = zeros(length(candidateICs), 1);
    for idx = 1:length(candidateICs)
        ic = candidateICs(idx);
        icData = reshape(EEG.icaact(ic, :, :), 1, []);
        fftData = fft(icData, 512); % Zero-padded FFT
        powerSpectrum = abs(fftData(1:256)).^2; % Single-sided
        freqRes = EEG.srate / 512;
        f = (0:255) * freqRes;
        thetaIdx = f >= thetaBand(1) & f <= thetaBand(2);
        thetaPower(idx) = mean(log(powerSpectrum(thetaIdx)));
    end
end

function EEG = performDipfit(EEG)
    EEG = pop_dipfit_settings(EEG, 'hdmfile', 'standard_BEM/standard_vol.mat', ...
        'mrifile', 'standard_BEM/standard_mri.mat', 'chanfile', ...
        'standard_BEM/elec/standard_1005.elc', 'coordformat', 'MNI');
    EEG = pop_multifit(EEG, 1:size(EEG.icaweights, 1), 'threshold', 100, 'dipplot', 'off');
end

function candidateICs = findComponentsNearACC(EEG) % Finds dipoles in Dorsal ACC using MNI Coordinate range X: -10 to +10 (midline), Y: 10 to 30 (anterior-posterior), Z: 30 to 50 (superior-inferior)
    candidateICs = [];
    for ic = 1:length(EEG.dipfit.model)
        dipole = EEG.dipfit.model(ic);
        % Check for residual variance and dACC location
        if dipole.rv < 0.2 && abs(dipole.posxyz(1)) <= 10 ... % Midline constraint
                && dipole.posxyz(2) >= 10 && dipole.posxyz(2) <= 30 ... % Anterior-posterior range
                && dipole.posxyz(3) >= 30 && dipole.posxyz(3) <= 50 % Dorsal Z constraint
            candidateICs = [candidateICs, ic];
        end
    end
end