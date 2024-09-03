        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The paper specifies  referencing to the right mastoid
% and bandpass filter 0.01 - 100 Hz *online* at data acquistion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Preprocessing
Action: Digitally filter data to remove frequencies above 50 Hz. Separate data into non-overlapping epochs time-locked to each Memorise, Ignore, and Probe letter onset. Remove epochs with high-amplitude, high-frequency artifacts.
Purpose: To clean the data by removing artifacts, ensuring only relevant brain activity is analysed.
ICA Decomposition
Action: Apply extended infomax ICA to the preprocessed EEG data to separate the data into independent components.
Purpose: To isolate independent brain activity components from mixed EEG signals.
Figure: Results will produce equivalent dipole locations for the frontal components similar to Fig. 1B.
Results Section: Frontal midline theta cluster.
Component Selection
Component Categorisation and Dipole Modeling
Action: Assess and categorise components as brain activity or non-brain artifact by visual inspection of their scalp topographies, time courses, and activation spectra. Compute equivalent current dipole models for each brain activity component map using a four-shell spherical head model in the DIPFIT toolbox.
Purpose: To ensure that only brain-related components are included in further analysis and to localise the source of brain activity within the head model.
Figures: Fig. 1B shows equivalent dipole locations for the fmθ cluster components.
Results Section: Frontal midline theta cluster.
Data Epoch Selection
Epoch Creation
Action: Separate the continuous data into 3-s epochs starting 1 s before and ending 2 s after each stimulus onset.
Purpose: To focus on specific periods of interest around stimulus presentations.
Figure: Results will identify epochs with significant theta power contributions during Memorise-letter trials, similar to Fig. 1E.
%}

%%
% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory
pathToEEGData = fullfile(pathToEEGLAB, 'brainvision_eeg'); % Assumes raw data is in the brainvision_eeg file

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% Get a list of all .vhdr files in the directory
EEGfiles = dir(fullfile(pathToEEGData, '*.vhdr'));
EEGfileNames = {EEGfiles.name}; % Extract the names of the files
disp(EEGfileNames); % Display the list of .vhdr files

%%%% Load dataset 
% ensure the brain vision extension is installed, can install
% throught the GUI
EEG = pop_loadbv(pathToEEGData, EEGfileNames{1}); % Just using first file for testing

% Add channel locations
chanlocs = EEG.chanlocs;
disp('chanlocs');
disp(chanlocs);
% Add channel locations using the 10-10 system
EEG = pop_chanedit(EEG, 'lookup', fullfile(pathToEEGLAB, 'channel_location_files', 'Standard-10-10-CUSTOM63.ced'));

% Visualise the cleaned EEG data using EEGLAB GUI
eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Not clean EEG data');
disp('Data visulised');

%%
% Reclassify events
EEG = reclassify_events(EEG);

%% 

% Bandpass filter 0.01 - 100 Hz 
% This should be done online but here as well.
EEG = pop_eegfiltnew(EEG, 0.01, 100);
disp('Bandpass filter comeplete');

%%%% Digitally filter data to remove frequencies above 50 Hz (low-pass filter with a 50 Hz cutoff)
EEG = pop_eegfiltnew(EEG, [], 50);
disp('Lowpass filter complete');


%%

% Apply automated artifact rejection algorithms
% Removes the requirement to do them mannually - we'll add in a visual
% inspection later. 
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 4, 'ChannelCriterion', 0.85, 'LineNoiseCriterion', 4, 'Highpass', 'off', ...
    'BurstCriterion', 20, 'WindowCriterion', 0.25, 'BurstRejection', 'on', 'Distance', 'Euclidian', 'WindowCriterionTolerances', [-Inf 7]);
disp('Automated artifact rejection complete');


%%

% Referencing is performed after initial filtering and artifact rejection but before ICA. This approach ensures that the data is cleaner and less noisy, for accurate ICA decomposition.
% Re-reference to the right mastoid 
EEG = pop_reref(EEG, 'M2');
% Re-reference to the right mastoid (assume 'M2' is the right mastoid in your channel locations)
%EEG = pop_reref(EEG, find(strcmp({EEG.chanlocs.labels}, 'M2')));
disp('Re-referencing complete');


%%

% Interpolation is performed after artifact rejection but before ICA. This ensures that any channels removed during artifact rejection are interpolated, providing a complete dataset for ICA.
% Interpolate removed channels based on the original channel locations
EEG = pop_interp(EEG, EEG.chanlocs);
% Interpolate removed channels based on the original channel locations
%EEG = pop_interp(EEG, EEG.chanlocs, 'spherical');
disp('Interpolation of channels complete');


%%
%{
%% JUST FOR TESTING

% Find the first occurrence of 's88' and change it to 's89'
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type, 's88')
        EEG.event(i).type = 's89';
        break;
    end
end

%%
%}

%%

%{
%%%%%%%%% JUST TO CHECK IF REMOVE INCoorrect FUNCTION WORKS
% save event data to a csv file for inspection
% Check if EEG.event is not empty
if ~isempty(EEG.event)
    % Initialize variables to store event details
    eventType = {EEG.event.type}';
    eventLatency = num2cell([EEG.event.latency]');
    eventDuration = num2cell([EEG.event.duration]');
    
    % Check if all events have a 'duration' field; if not, adjust accordingly
    if isfield(EEG.event, 'duration')
        % Combine into a table (if 'duration' is available for all events)
        eventsTable = table(eventType, eventLatency, eventDuration, 'VariableNames', {'Type', 'Latency', 'Duration'});
    else
        % If 'duration' is not a universal field, create a table without it
        eventsTable = table(eventType, eventLatency, 'VariableNames', {'Type', 'Latency'});
    end
    
    % Define your output CSV file name
    outputCSVFileName = 'Before_incorrect_EEGEvents.csv';
    
    % Write the table to a CSV file
    writetable(eventsTable, outputCSVFileName);
    
    disp(['Events have been exported to ' outputCSVFileName]);
else
    disp('No events found in the EEG structure.');
end
%%%%%%%%%%%%%%%%
%}

%%

%%%%%%%%%%%%%% Define Trigger Events 
% welcome = {‘s0’} 
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

%%%%%%%%%%%%%% END OF Define Trigger Events 

% Remove incorrect trials 
EEG = removeIncorrectTrials(EEG, fixationCrossEvent, correctAnswerEvent, wrongAnswerEvent);
disp('Incorrect trails removed');


%%
%{
% save event data to a csv file for inspection
% Check if EEG.event is not empty
if ~isempty(EEG.event)
    % Initialize variables to store event details
    eventType = {EEG.event.type}';
    eventLatency = num2cell([EEG.event.latency]');
    eventDuration = num2cell([EEG.event.duration]');
    
    % Check if all events have a 'duration' field; if not, adjust accordingly
    if isfield(EEG.event, 'duration')
        % Combine into a table (if 'duration' is available for all events)
        eventsTable = table(eventType, eventLatency, eventDuration, 'VariableNames', {'Type', 'Latency', 'Duration'});
    else
        % If 'duration' is not a universal field, create a table without it
        eventsTable = table(eventType, eventLatency, 'VariableNames', {'Type', 'Latency'});
    end
    
    % Define your output CSV file name
    outputCSVFileName = 'After_incorrect_EEGEvents.csv';
    
    % Write the table to a CSV file
    writetable(eventsTable, outputCSVFileName);
    
    disp(['Events have been exported to ' outputCSVFileName]);
else
    disp('No events found in the EEG structure.');
end
%}

%%

% Apply extended infomax ICA with default extended-mode parameters and a stopping weight change of 1e-7
EEG = pop_runica(EEG, 'extended', 1, 'stop', 1e-7);
disp('ICA conducted');

% ensure the Fieldtrip-Lite extension is installed, can install
% throught the GUI
% Dipole fitting for component selection
EEG = pop_dipfit_settings(EEG, 'hdmfile', fullfile(pathToEEGLAB, 'standard_BEM/standard_vol.mat'), ...
                          'coordformat', 'MNI', 'mrifile', fullfile(pathToEEGLAB, 'standard_BEM/standard_mri.mat'), ...
                          'chanfile', fullfile(pathToEEGLAB, 'standard_BEM/elec/standard_1005.elc'), 'coord_transform', [0 0 0 0 0 0 1 1 1]);
EEG = pop_multifit(EEG, [], 'threshold', 15, 'plotopt', {'normlen', 'on'});
disp('Dipole fitting complete');

% Automatically categorise components using ICLabel
EEG = pop_iclabel(EEG, 'default');
disp('ICA componets labled');

% Automatically flag and remove artifact components
EEG = pop_icflag(EEG, [NaN NaN; 0.9 1; 0.9 1; NaN NaN; NaN NaN; NaN NaN; NaN NaN]);
disp('Automatically flag and remove artifact components');

% Remove flagged components
EEG = pop_subcomp(EEG, [], 0);
disp('Remove flagged components');

%%
% Extract all event types into a cell array
%allEventTypes = {EEG.event.type};
% Identify unique event types
%uniqueEventTypes = unique(allEventTypes);
% Print the unique event types
%disp('Unique event types:');
%disp(uniqueEventTypes);
%%
% Visualise the cleaned EEG data using EEGLAB GUI
eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Prepreprocessed EEG Data After Artifact Rejection');

%%

% Extract the base filename without extension from EEG.comments
fileBase = extractBetween(EEG.comments, 'Original file: ', '.eeg');
fileBase = fileBase{1}; % Convert the cell array to a string

%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '1_preprocessed_data'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '1_preprocessed_data'));
end

%Save preprocessed data with incorrect trails removed
pop_saveset(EEG, 'filename', [fileBase '_full_cnt_data_incorrect_removed.set'], 'filepath', fullfile(EEG.filepath, '1_preprocessed_data'));
disp('Data saved');


%%
% Segment data into epochs time-locked to the onsets of Memorise, Ignore, and Probe stimuli
% Segmented after preprocessing since in paper the epochs were to save
% computational resources.
%EEG_memorise = pop_epoch(EEG, {'s30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 's50', 's51', 's52', 's53', 's54', 's55', 's56', 's57', 's70', 's71', 's72', 's73', 's74', 's75', 's76', 's77'}, [-1 2], 'epochinfo', 'yes');
%EEG_ignore = pop_epoch(EEG, {'s40', 's41', 's42', 's43', 's44', 's45', 's46', 's47', 's60', 's61', 's62', 's63', 's64', 's65', 's66', 's67', 's80', 's81', 's82', 's83', 's84', 's85', 's86', 's87'}, [-1 2], 'epochinfo', 'yes');
%EEG_probe = pop_epoch(EEG, {'s38', 's39', 's58', 's59', 's78', 's79'}, [-1 2], 'epochinfo', 'yes');
EEG_memorise = pop_epoch(EEG, memoriseLetterEvents, [-1 2], 'epochinfo', 'yes');
EEG_ignore = pop_epoch(EEG, ignoreLetterEvents, [-1 2], 'epochinfo', 'yes');
EEG_probe = pop_epoch(EEG, probeEvents, [-1 2], 'epochinfo', 'yes');

% Segment longer epochs for fixation (5s) and memory maintenance (varies between 2s and 4s)
EEG_fixation = pop_epoch(EEG, fixationCrossEvent, [-1 4], 'epochinfo', 'yes');
% Custom function to extract variable-length maintenance epochs
EEG_maintenance = extractVariableEpochs(EEG, variableMaintenancePeriodEvent, probeEvents);
disp('Epochs created');

%%

%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '2_epoch_data/probe'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '2_epoch_data/probe'));
end
% Save epoch
pop_saveset(EEG_probe, 'filename', [fileBase '_probe.set'], 'filepath', fullfile(EEG.filepath, '2_epoch_data/probe'));

%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '2_epoch_data/memorise'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '2_epoch_data/memorise'));
end
% Save epoch
pop_saveset(EEG_memorise, 'filename', [fileBase '_memorise.set'], 'filepath', fullfile(EEG.filepath, '2_epoch_data/memorise'));

%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '2_epoch_data/ignore'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '2_epoch_data/ignore'));
end
% Save epoch
pop_saveset(EEG_ignore, 'filename', [fileBase '_ignore.set'], 'filepath', fullfile(EEG.filepath, '2_epoch_data/ignore'));

%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '2_epoch_data/fixation'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '2_epoch_data/fixation'));
end
% Save epoch
pop_saveset(EEG_fixation, 'filename', [fileBase '_fixation.set'], 'filepath', fullfile(EEG.filepath, '2_epoch_data/fixation'));


%%

% Check if the directory exists
if ~exist(fullfile(EEG.filepath, '2_epoch_data/maintenance'), 'dir')
    % Create the directory if it doesn't exist
    mkdir(fullfile(EEG.filepath, '2_epoch_data/maintenance'));
end
% Save epoch
pop_saveset(EEG_maintenance, 'filename', [fileBase '_maintenance.set'], 'filepath', fullfile(EEG.filepath, '2_epoch_data/maintenance'));


%%

disp('Epochs saved');
disp('Pre-processing complete and epoched datasets saved.');

%%
%{
% Compute and plot equivalent dipoles for ICA components

% Load ICA weights and sphere
load('ica_weights.mat', 'ica_weights', 'ica_sphere', 'ica_components');

% Load channel locations
chanlocs = load('channel_locations.mat');

% Compute equivalent dipoles
dipfit_settings = struct('model', 'single', 'rv_threshold', 0.15, 'vol', []);
[~, dipoles] = fit_dipoles(ica_weights, ica_sphere, chanlocs, dipfit_settings);

% Extract frontal components (example criteria, modify as needed)
frontal_dipoles = find_dipoles_in_region(dipoles, 'frontal');

% Plot dipoles
figure;
plot_dipoles(frontal_dipoles, 'color', 'pink');
title('Equivalent Dipole Locations for Frontal Components');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

% Save the figure
savefig('frontal_dipoles.fig');
%}