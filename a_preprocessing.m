%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MATLAB Script for EEG Preprocessing using EEGLAB
% Author: Paige Metcalfe

%{
Input: Raw EEG data
Output: Cleaned data and epoch data files 
Summary of steps:
%}


% Set variables
clear;

% Set directories
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory
pathToEEGData = fullfile(pathToEEGLAB, 'brainvision_eeg'); % Assumes raw data is in the brainvision_eeg file
outputFolder = fullfile(pathToEEGLAB, 'a_preprocessed_data'); % Output folder for preprocessed data
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% Initialise EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Get a list of all .vhdr files in the directory
EEGfiles = dir(fullfile(pathToEEGData, '*.vhdr'));
% Raise error if no files found
if isempty(EEGfiles)
    error('No .vhdr files found in the specified directory: %s', pathToEEGData);
end
EEGfileNames = {EEGfiles.name}; % Extract the names of the files
disp(EEGfileNames); % Display the list of .vhdr files

% Define Trigger Events 
blockRestEvent = {'s19', 's29'}; % block_start % s19 = short, s29 = long 
interTrialRestEvent = {'s15'}; % rest
fixationCrossEvent = {'s10'}; % fixation_cross 
memoriseLetterEvents = {'s30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 's50', 's51', 's52', 's53', 's54', 's55', 's56', 's57', 's70', 's71', 's72', 's73', 's74', 's75', 's76', 's77'}; % trail_3_letters  
% The first number 3,5 or 7 indicates memory load and the following number is letter position (zero-indexed). 
ignoreLetterEvents = {'s40', 's41', 's42', 's43', 's44', 's45', 's46', 's47', 's60', 's61', 's62', 's63', 's64', 's65', 's66', 's67', 's80', 's81', 's82', 's83', 's84', 's85', 's86', 's87'}; % trail_3_letters  
% The first number 4,6 or 8 indicates (memory load+1) just so can differentiate between the memorise and ignore. The following number is letter position (zero-indexed). 
variableMaintenancePeriodEvent = {'s11'}; % memorise 
probeEvents = {'s38', 's39', 's58', 's59', 's78', 's79'}; % probe % s38 means 3 memory load and probe letter was in set and s39 is same memory load but probe letter NOT in set  
correctAnswerEvent = {'s88'}; % response_validation  
wrongAnswerEvent = {'s89'}; % response_validation 
blockEndEvent = {'s91', 's92'}; % blockEnd % s91 = short, s92 = long 
endEvent = {'s99'}; % End 

% Loop through files for preprocessing
for i = 1:length(EEGfileNames)
    % Load dataset 
    % ensure the brain vision extension is installed, can install
    % throught the GUI
    fprintf('Processing file: %s\n', EEGfileNames{i});
    try
        EEG = pop_loadbv(pathToEEGData, EEGfileNames{i});
    catch ME
        warning('Error loading file: %s\n%s', EEGfileNames{i}, ME.message);
        continue; % Skip to the next file
    end

    % Add channel locations
    chanlocs = EEG.chanlocs;
    disp('chanlocs');
    disp(chanlocs);
    % Add channel locations using the 10-10 system
    EEG = pop_chanedit(EEG, 'lookup', fullfile(pathToEEGLAB, 'channel_location_files', 'Standard-10-10-CUSTOM63.ced'));
    
    % Visualise the cleaned EEG data using EEGLAB GUI
    eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Not clean EEG data');
    
    % Reclassify events
    EEG = reclassify_events(EEG);

    % Bandpass filter 0.01 - 50 Hz 
    EEG = pop_eegfiltnew(EEG, 0.01, 50);
    disp('Bandpass filter comeplete');

    % Apply automated artifact rejection algorithms
    try
        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 4, 'ChannelCriterion', 0.85, 'LineNoiseCriterion', 4, 'Highpass', 'off', ...
            'BurstCriterion', 20, 'WindowCriterion', 0.25, 'BurstRejection', 'on', 'Distance', 'Euclidian', 'WindowCriterionTolerances', [-Inf 7]);
    catch ME
        warning('Artifact rejection failed for file: %s\n%s', EEGfileNames{i}, ME.message);
        continue; % Skip problematic datasets
    end
    disp('Automated artifact rejection complete');

    % Re-reference to the right mastoid 
    if any(strcmp({EEG.chanlocs.labels}, 'M2'))
        EEG = pop_reref(EEG, {'M2'});
        disp('Re-referencing to M2 complete.');
    else
        warning('M2 channel not found. Re-referencing skipped for file: %s', EEGfileNames{i});
    end

    disp('Re-referencing complete');

    % Interpolate removed channels based on the original channel locations
    EEG = pop_interp(EEG, EEG.chanlocs);
    disp('Interpolation of channels complete');

    % Remove incorrect trials 
    EEG = removeIncorrectTrials(EEG, fixationCrossEvent, correctAnswerEvent, wrongAnswerEvent);
    disp('Incorrect trails removed');

    % Apply extended infomax ICA with default extended-mode parameters and a stopping weight change of 1e-7
    EEG = pop_runica(EEG, 'extended', 1, 'stop', 1e-7);
    disp('ICA conducted');
    
    %{
    % TODO: Set up installation!
    % ensure DIPFIT plugin is installed
    % Apply DIPFIT settings to localise ICs anatomically.
    EEG = pop_dipfit_settings(EEG, 'hdmfile', 'standard_BEM/standard_vol.mat', ...
        'coordformat', 'MNI', 'mrifile', 'standard_BEM/standard_mri.mat', ...
        'chanfile', 'standard_BEM/elec/standard_1005.elc', ...
        'coord_transform', [0 0 0 0 0 0 1 1 1]);

    % Perform dipole fitting
    EEG = pop_multifit(EEG, 1:EEG.nbchan, 'threshold', 100, 'dipplot', 'off', 'plotopt', {'normlen' 'on'});

    % Remove components with >15% residual variance or outside brain
    dipoles = [EEG.dipfit.model.rv];
    inside = [EEG.dipfit.model.inside];
    reject = find(dipoles > 0.15 | ~inside);
    EEG = pop_subcomp(EEG, reject, 0);

    disp('Dipole fitting completed. Components with >15% residual variance or outside the brain removed.');

    % Optional: Plot dipoles for inspection
    %dipplot(EEG, 'components', 1:length(EEG.icachansind), 'view', [0 0], 'normlen', 'on');
    %disp('Dipole locations visualised.');

    %}

    % Automatically categorise components using ICLabel
    EEG = pop_iclabel(EEG, 'default');
    disp('ICA componets labled');
    
    % Automatically flag and remove artifact components
    EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
    disp('Automatically flag and remove artifact components');
    
    % Remove flagged components
    EEG = pop_subcomp(EEG, [], 0);
    disp('Remove flagged components');

    % Visualise the cleaned EEG data using EEGLAB GUI
    eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Prepreprocessed EEG Data After Artifact Rejection');

    % Check everything is saved
    disp(EEG.icaweights); % Display the ICA weights matrix
    disp(EEG.icasphere);  % Display the ICA sphere matrix
    disp(EEG.icaact);     % Display matrix of independent component activations
    disp(EEG.etc.ic_classification);  % Display the classification of components (e.g., ICLabel results)
    disp(EEG.dipfit);     % Display the dipole model information (if DIPFIT is applied)



    % Extract the base filename without extension from EEG.comments
    fileBase = extractBetween(EEG.comments, 'Original file: ', '.eeg');
    fileBase = fileBase{1}; % Convert the cell array to a string

    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '1_cleaned_continuous'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '1_cleaned_continuous'));
    end
    
    %Save preprocessed data with incorrect trails removed
    pop_saveset(EEG, 'filename', [fileBase '_cleaned_continuous.set'], 'filepath', fullfile(outputFolder, '1_cleaned_continuous'));
    disp('Cleaned continous data saved');

    % Segment data into epochs time-locked to the onsets of Memorise, Ignore, and Probe stimuli
    EEG_memorise = pop_epoch(EEG, memoriseLetterEvents, [-1 2], 'epochinfo', 'yes');
    EEG_ignore = pop_epoch(EEG, ignoreLetterEvents, [-1 2], 'epochinfo', 'yes');
    EEG_probe = pop_epoch(EEG, probeEvents, [-1 2], 'epochinfo', 'yes');
    % Segment longer epochs for fixation (5s) and memory maintenance (varies between 2s and 4s)
    EEG_fixation = pop_epoch(EEG, fixationCrossEvent, [-1 4], 'epochinfo', 'yes');
    % Custom function to extract variable-length maintenance epochs
    EEG_maintenance = extractVariableEpochs(EEG, variableMaintenancePeriodEvent, probeEvents);
    disp('Epochs created');

    % Save epochs

    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '2_epoch_data/probe'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '2_epoch_data/probe'));
    end
    % Save epoch
    pop_saveset(EEG_probe, 'filename', [fileBase '_probe.set'], 'filepath', fullfile(outputFolder, '2_epoch_data/probe'));
    
    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '2_epoch_data/memorise'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '2_epoch_data/memorise'));
    end
    % Save epoch
    pop_saveset(EEG_memorise, 'filename', [fileBase '_memorise.set'], 'filepath', fullfile(outputFolder, '2_epoch_data/memorise'));
    
    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '2_epoch_data/ignore'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '2_epoch_data/ignore'));
    end
    % Save epoch
    pop_saveset(EEG_ignore, 'filename', [fileBase '_ignore.set'], 'filepath', fullfile(outputFolder, '2_epoch_data/ignore'));
    
    
    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '2_epoch_data/fixation'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '2_epoch_data/fixation'));
    end
    % Save epoch
    pop_saveset(EEG_fixation, 'filename', [fileBase '_fixation.set'], 'filepath', fullfile(outputFolder, '2_epoch_data/fixation'));
    
    
    % Check if the directory exists
    if ~exist(fullfile(outputFolder, '2_epoch_data/maintenance'), 'dir')
        % Create the directory if it doesn't exist
        mkdir(fullfile(outputFolder, '2_epoch_data/maintenance'));
    end
    % Save epoch
    pop_saveset(EEG_maintenance, 'filename', [fileBase '_maintenance.set'], 'filepath', fullfile(outputFolder, '2_epoch_data/maintenance'));


    % Display loop progress
    disp(i)
    disp("of")
    disp(length(EEGfileNames))
    dips("=======================")

    % Clear variables for next iteration
    clear EEG EEG_memorise EEG_ignore EEG_probe EEG_fixation EEG_maintenance;

end

disp('Pre-processing complete and epoched datasets saved.');




