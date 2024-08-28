%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Spectral Analysis
Action: Analyse spectral responses during memory maintenance and after Probe letter onsets.
Purpose: To examine changes in theta power related to memory load and probe response.
Figures: Fig. 8A shows smoothed trial-by-trial time courses of log power at 6 Hz before and after Probe presentation. 
         Fig. 8B shows smoothed trial-by-trial time courses of log power at 3 Hz after Probe presentation.
Results Section: Spectral response to probe letters.
%}

%% Spectral Analysis of Memory Maintenance and Probe Responses

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contains custom functions
addpath('utils');

% Paths to memory maintenance and probe epochs
maintenance_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/maintenance');
probe_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/probe');

% Get a list of maintenance and probe .set files
Maintenance_epoch_files = dir(fullfile(maintenance_epoch_filepath, '*.set'));
Probe_epoch_files = dir(fullfile(probe_epoch_filepath, '*.set'));

% Extract the filenames
Maintenance_epoch_fileNames = {Maintenance_epoch_files.name};
Probe_epoch_fileNames = {Probe_epoch_files.name};

% Display the list of files
disp('Maintenance Files:');
disp(Maintenance_epoch_fileNames);
disp('Probe Files:');
disp(Probe_epoch_fileNames);

% Initialize variables for storing results
theta_power_maintenance = [];
theta_power_probe_3Hz = [];
theta_power_probe_6Hz = [];
subject_ids = {};

% Parameters for spectral analysis
theta_freq_range = [5 7]; % Theta frequency range for memory maintenance
probe_3Hz = 3; % 3 Hz for post-probe analysis
probe_6Hz = 6; % 6 Hz for pre- and post-probe analysis
time_window_maintenance = [-1000 2000]; % Time window for maintenance period (in ms)
time_window_probe = [-500 1000]; % Time window for probe response (in ms)
smoothing_window = 5; % Smoothing window size for time courses

% Loop through the subjects
for subj = 1:length(Maintenance_epoch_fileNames)
    % Load the maintenance and probe data
    EEG_maintenance = pop_loadset('filename', Maintenance_epoch_fileNames{subj}, 'filepath', maintenance_epoch_filepath);
    EEG_probe = pop_loadset('filename', Probe_epoch_fileNames{subj}, 'filepath', probe_epoch_filepath);

    % Display subject ID
    subject_id_maintenance = EEG_maintenance.subject;
    subject_id_probe = EEG_probe.subject;
    
    % Ensure subject IDs match
    if ~strcmp(subject_id_maintenance, subject_id_probe)
        error(['Subject ID mismatch between Maintenance and Probe data: ' subject_id_maintenance ' vs ' subject_id_probe]);
    end

    subject_ids{end+1} = subject_id_maintenance;
    
    %% Memory Maintenance Period: Theta Power Analysis

    % Compute ERSP for maintenance period within the theta range
    [ersp_maintenance, ~, ~, times] = newtimef(EEG_maintenance.icaact, EEG_maintenance.pnts, time_window_maintenance, ...
                                               EEG_maintenance.srate, [3 0.5], 'freqs', theta_freq_range, ...
                                               'nfreqs', length(theta_freq_range), 'timesout', 200, ...
                                               'baseline', NaN, 'plotersp', 'off');
    
    % Extract theta power (5-7 Hz) during maintenance period
    theta_power_maintenance{subj} = mean(ersp_maintenance, 1);

    %% Probe Response Period: Theta Power Analysis

    % Compute ERSP for probe period at 3 Hz
    [ersp_probe_3Hz, ~, ~, times_probe] = newtimef(EEG_probe.icaact, EEG_probe.pnts, time_window_probe, ...
                                                   EEG_probe.srate, [3 0.5], 'freqs', probe_3Hz, ...
                                                   'nfreqs', 1, 'timesout', 200, ...
                                                   'baseline', NaN, 'plotersp', 'off');

    % Extract theta power at 3 Hz during probe period
    theta_power_probe_3Hz{subj} = ersp_probe_3Hz;

    % Compute ERSP for probe period at 6 Hz
    [ersp_probe_6Hz, ~, ~, ~] = newtimef(EEG_probe.icaact, EEG_probe.pnts, time_window_probe, ...
                                         EEG_probe.srate, [3 0.5], 'freqs', probe_6Hz, ...
                                         'nfreqs', 1, 'timesout', 200, ...
                                         'baseline', NaN, 'plotersp', 'off');

    % Extract theta power at 6 Hz during probe period
    theta_power_probe_6Hz{subj} = ersp_probe_6Hz;

    %% Smoothing and Visualization of Time Courses

    % Apply smoothing to the theta power time courses
    smoothed_theta_6Hz = smoothdata(theta_power_probe_6Hz{subj}, 'movmean', smoothing_window);
    smoothed_theta_3Hz = smoothdata(theta_power_probe_3Hz{subj}, 'movmean', smoothing_window);

    % Visualize the results for each subject
    figure;
    subplot(2, 1, 1);
    plot(times_probe, smoothed_theta_6Hz);
    title(['Smoothed Theta Power (6 Hz) - Probe Response (Subject ' subject_ids{subj} ')']);
    xlabel('Time (ms)');
    ylabel('Log Power (6 Hz)');
    xlim(time_window_probe);
    ylim([-3 3]); % Adjust based on data range
    grid on;
    
    subplot(2, 1, 2);
    plot(times_probe, smoothed_theta_3Hz);
    title(['Smoothed Theta Power (3 Hz) - Probe Response (Subject ' subject_ids{subj} ')']);
    xlabel('Time (ms)');
    ylabel('Log Power (3 Hz)');
    xlim(time_window_probe);
    ylim([-3 3]); % Adjust based on data range
    grid on;
end

% Save the results for further analysis
save(fullfile(pathToEEGLAB, '7_spectral_analysis', 'theta_power_analysis.mat'), 'theta_power_maintenance', 'theta_power_probe_3Hz', 'theta_power_probe_6Hz', 'subject_ids');

disp('Spectral analysis of memory maintenance and probe responses complete.');
