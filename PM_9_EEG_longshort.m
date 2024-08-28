%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Long and Short Letter Presentations
Comparison Analysis
Action: Compare single-trial and mean activities of components between Long and Short letter presentations.
Purpose: To determine if presentation length affects spectral power changes.
Figures: Fig. 7 compares short versus long letter presentation.
Results Section: Short presentation control condition.
%}

%% Long vs Short Letter Presentations

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contains custom functions
addpath('utils');

% Paths to Long and Short letter presentation data
long_presentation_filepath = fullfile(pathToEEGLAB, '2_epoch_data/long_presentation');
short_presentation_filepath = fullfile(pathToEEGLAB, '2_epoch_data/short_presentation');

% Get a list of Long and Short letter presentation .set files
Long_epoch_files = dir(fullfile(long_presentation_filepath, '*.set'));
Short_epoch_files = dir(fullfile(short_presentation_filepath, '*.set'));

% Extract the filenames
Long_epoch_fileNames = {Long_epoch_files.name};
Short_epoch_fileNames = {Short_epoch_files.name};

% Display the list of files
disp('Long Presentation Files:');
disp(Long_epoch_fileNames);
disp('Short Presentation Files:');
disp(Short_epoch_fileNames);

% Check if there are 9 subjects for both conditions
if length(Long_epoch_fileNames) ~= 9 || length(Short_epoch_fileNames) ~= 9
    error('Mismatch in the number of Long and Short presentation files. Ensure there are 9 files for each condition.');
end

% Initialize variables for storing results
component_activities_long = [];
component_activities_short = [];
subject_ids = {};

% Loop through the subjects
for subj = 1:length(Long_epoch_fileNames)
    % Load the Long and Short presentation data
    EEG_long = pop_loadset('filename', Long_epoch_fileNames{subj}, 'filepath', long_presentation_filepath);
    EEG_short = pop_loadset('filename', Short_epoch_fileNames{subj}, 'filepath', short_presentation_filepath);

    % Display subject ID
    subject_id_long = EEG_long.subject;
    subject_id_short = EEG_short.subject;
    
    % Ensure subject IDs match
    if ~strcmp(subject_id_long, subject_id_short)
        error(['Subject ID mismatch between Long and Short presentations: ' subject_id_long ' vs ' subject_id_short]);
    end

    subject_ids{end+1} = subject_id_long;
    
    % Concatenate the Long and Short presentation data for ICA
    EEG_combined = pop_mergeset(EEG_long, EEG_short);
    
    % Perform ICA on the combined data
    EEG_combined = pop_runica(EEG_combined, 'extended', 1, 'stop', 1e-7);
    
    % Extract components whose locations match the clustered Long-presentation components
    % Assuming a custom function find_matching_components is available for this
    matching_components = find_matching_components(EEG_combined, EEG_long);

    % Store single-trial and mean activities of matching components for both conditions
    component_activities_long{subj} = EEG_long.icaact(matching_components, :, :);
    component_activities_short{subj} = EEG_short.icaact(matching_components, :, :);
end

% Compare single-trial and mean activities between Long and Short conditions
for subj = 1:length(subject_ids)
    long_activity = component_activities_long{subj};
    short_activity = component_activities_short{subj};
    
    % Calculate mean activities for both conditions
    mean_long_activity = mean(long_activity, 3);
    mean_short_activity = mean(short_activity, 3);
    
    % Calculate differences in spectral power
    power_diff = mean_long_activity - mean_short_activity;
    
    % Visualize the comparison for each subject
    figure;
    subplot(2, 1, 1);
    plot(mean_long_activity);
    title(['Mean Activity - Long Presentation (Subject ' subject_ids{subj} ')']);
    xlabel('Component');
    ylabel('Activity');
    
    subplot(2, 1, 2);
    plot(mean_short_activity);
    title(['Mean Activity - Short Presentation (Subject ' subject_ids{subj} ')']);
    xlabel('Component');
    ylabel('Activity');
    
    % Optionally, visualize the power difference
    figure;
    plot(power_diff);
    title(['Difference in Mean Activity (Long - Short) (Subject ' subject_ids{subj} ')']);
    xlabel('Component');
    ylabel('Power Difference');
end

% Save results for further analysis
save(fullfile(pathToEEGLAB, '6_comparisons', 'long_vs_short_activities.mat'), 'component_activities_long', 'component_activities_short', 'subject_ids');

disp('Comparison of Long vs Short letter presentations complete.');
