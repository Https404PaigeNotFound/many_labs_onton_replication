%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Long and Short Letter Presentations
Comparison Analysis
Action: Apply Independent Component Analysis (ICA) to combined Long and Short letter presentation data, then compare ERSPs of specific independent components (ICs) between conditions.
Purpose: To determine if presentation length affects spectral power changes, specifically in theta and low-beta bands, by focusing on the activity of independent components related to cognitive processing.
Figures: Fig. 7 compares ERSP differences between short versus long letter presentations for specific ICs.
Results Section: Short presentation control condition.
%}


%% Long vs Short Letter Presentations based on ERSPs of Independent Components

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAB as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contains custom functions
addpath('utils');

% Paths to Long and Short letter presentation EEG data
long_presentation_filepath = fullfile(pathToEEGLAB, '2_epoch_data', 'long_presentation');
short_presentation_filepath = fullfile(pathToEEGLAB, '2_epoch_data', 'short_presentation');

% Initialize variables for storing results
ersp_long_all = struct();
ersp_short_all = struct();
times_all = [];
freqs_all = [];
subject_ids = {};

% Get a list of Long and Short letter presentation .set files
Long_epoch_files = dir(fullfile(long_presentation_filepath, '*.set'));
Short_epoch_files = dir(fullfile(short_presentation_filepath, '*.set'));

% Check if there are matching numbers of subjects for both conditions
if length(Long_epoch_files) ~= length(Short_epoch_files)
    error('Mismatch in the number of Long and Short presentation files.');
end

% Loop through the subjects
for subj = 1:length(Long_epoch_files)
    % Load the Long and Short presentation EEG data
    EEG_long = pop_loadset('filename', Long_epoch_files{subj}, 'filepath', long_presentation_filepath);
    EEG_short = pop_loadset('filename', Short_epoch_files{subj}, 'filepath', short_presentation_filepath);

    % Ensure subject IDs match
    subject_id_long = EEG_long.subject;
    subject_id_short = EEG_short.subject;
    if ~strcmp(subject_id_long, subject_id_short)
        error(['Subject ID mismatch between Long and Short presentations: ' subject_id_long ' vs ' subject_id_short]);
    end
    
    subject_ids{end+1} = subject_id_long;
    
    % Concatenate the Long and Short presentation data for ICA
    EEG_combined = pop_mergeset(EEG_long, EEG_short);
    
    % Perform ICA on the combined data
    EEG_combined = pop_runica(EEG_combined, 'extended', 1, 'stop', 1e-7);
    
    % Identify the independent components (ICs) of interest (e.g., fmθ cluster)
    % Check the function
    ICs_of_interest = identify_components_of_interest(EEG_combined);

    % Loop through each independent component of interest
    for ic = 1:length(ICs_of_interest)
        % Calculate ERSPs for the Long presentation trials
        [ersp_long, times, freqs] = pop_newtimef(EEG_combined, 0, ICs_of_interest(ic), ...
            [-1000 1998], [3 0.5], 'baseline', NaN, 'plotersp', 'off', ...
            'plotitc', 'off', 'plotphase', 'off', ...
            'title', ['Long Presentation ERSP (IC ' num2str(ICs_of_interest(ic)) ')']);
        
        % Calculate ERSPs for the Short presentation trials
        [ersp_short, ~, ~] = pop_newtimef(EEG_combined, 0, ICs_of_interest(ic), ...
            [-1000 1998], [3 0.5], 'baseline', NaN, 'plotersp', 'off', ...
            'plotitc', 'off', 'plotphase', 'off', ...
            'title', ['Short Presentation ERSP (IC ' num2str(ICs_of_interest(ic)) ')']);
        
        % Store the ERSP data
        ersp_long_all.(subject_id_long).(['IC_' num2str(ICs_of_interest(ic))]) = ersp_long;
        ersp_short_all.(subject_id_short).(['IC_' num2str(ICs_of_interest(ic))]) = ersp_short;
    end
    
    % Store times and frequencies assuming consistency across subjects
    times_all = times;  
    freqs_all = freqs;  
end

% Compare ERSPs between Long and Short conditions for each IC
for subj = 1:length(subject_ids)
    subject_id = subject_ids{subj};
    
    IC_names = fieldnames(ersp_long_all.(subject_id));
    for ic_idx = 1:length(IC_names)
        IC_name = IC_names{ic_idx};
        
        long_ersp = ersp_long_all.(subject_id).(IC_name);
        short_ersp = ersp_short_all.(subject_id).(IC_name);
        
        % Calculate differences in ERSPs
        ersp_diff = long_ersp - short_ersp;
        
        % Visualize the ERSP difference for each IC
        figure;
        imagesc(times_all, freqs_all, ersp_diff);
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
        title(['ERSP Difference (Long - Short) - ' IC_name ' (Subject ' subject_id ')']);
        colorbar;
        
        % Optional: Save the ERSP difference images
        saveas(gcf, fullfile(pathToEEGLAB, '6_comparisons', ['ersp_diff_' IC_name '_' subject_id '.png']));
    end
end

% Save ERSP results for further analysis
save(fullfile(pathToEEGLAB, '6_comparisons', 'long_vs_short_ersp_differences.mat'), ...
    'ersp_long_all', 'ersp_short_all', 'subject_ids', 'times_all', 'freqs_all');

disp('Comparison of Long vs Short letter presentations based on ERSPs of Independent Components complete.');
