%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Component Power Spectra
Power Spectrum Calculation
Action: Calculate power spectra for each task epoch by averaging FFT spectra computed using data window lengths of 512 points without zero-padding.
Purpose: To analyse the frequency components of brain activity during different task conditions.
Figure:  Fig. 1C shows power spectra of individual components in the fmθ cluster. Fig. 1D shows power spectra at the Fz electrode.
Results Section: Frontal midline theta cluster.
%}

%%

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');

% paths to each epoch
probe_epoch_filepath = fullfile(pathToEEGLAB, '/2_epoch_data/probe');
memorise_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/memorise');
ignore_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/ignore');
fixation_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/fixation');
maintenance_epoch_filepath = fullfile(pathToEEGLAB, '2_epoch_data/maintenance');

%%

% EXAMPLE with probe epochs
% Get a list of all .set files in the directory
Probe_epoch_files = dir(fullfile(probe_epoch_filepath, '*.set'));
Probe_epoch_fileNames = {Probe_epoch_files.name}; % Extract the names of the files
disp(Probe_epoch_fileNames); % Display the list of .set files

% Initialise arrays to hold the spectra for all subjects
all_probe_spectra = [];
all_subject_ids = [];

% Loop over all participants 
for subj = 1:length(Probe_epoch_fileNames)
    % Load subject-specific epoch data
    EEG_Probe_epoch = pop_loadset('filename', Probe_epoch_fileNames{subj}, 'filepath', probe_epoch_filepath);
    disp('Data loaded')
   
    % Compute independent component activations with custom function
    EEG_Probe_epoch = compute_ica_activations(EEG_Probe_epoch);
    disp('independent component activations computed')
    
    % Compute the power spectra using FFT with zero-padding to 512 points
    winLength = 512; % Window length of 512 points (without zero-padding)
    [probe_spectra, probe_freqs] = spectopo(EEG_Probe_epoch.icaact, 0, EEG_Probe_epoch.srate, 'component', 'all', 'winsize', winLength, 'overlap', 0, 'plot', 'off');
    disp('power spectra computed')
    
    % Debugging: Check the output of spectopo
    if subj == 1
        disp('Spectopo output (first subject):');
        disp(probe_spectra);
        disp(probe_freqs);
    end

    % Convert power to log power
    log_probe_spectra = 10 * log10(probe_spectra);
    disp('power to log power converted ')
    
    % Normalise the spectra by subtracting mean log power from single-trial log power
    mean_log_probe_spectra = mean(log_probe_spectra, 2); % Calculate mean log power across trials for each component and frequency
    probe_spectra_norm = log_probe_spectra - mean_log_probe_spectra; % Subtract mean log power
    disp('spectra normalised')
    
    % Restrict the frequency range to 2-30 Hz
    freq_idx = probe_freqs >= 2 & probe_freqs <= 30;
    probe_freqs_norm = probe_freqs(freq_idx);
    probe_spectra_norm = probe_spectra_norm(:, freq_idx, :);
    disp('frequncy restricted')
    
    % Aggregate data across subjects
    all_probe_spectra = cat(4, all_probe_spectra, probe_spectra_norm); % Concatenate along the 4th dimension (subjects)
    all_subject_ids = [all_subject_ids, subj];
    disp('data aggregated')
end

%%
% Check if only one subject is present
if length(Probe_epoch_fileNames) == 1
    disp('Only one subject data file loaded. Skipping ANOVA.');
    % Visualise the component power spectra for the single subject
    numComponents = size(probe_spectra_norm, 1);
    figure;
    for comp = 1:numComponents
        subplot(ceil(sqrt(numComponents)), ceil(sqrt(numComponents)), comp);
        plot(probe_freqs_norm, squeeze(mean(probe_spectra_norm(comp, :, :), 3))); % Plot average across trials
        title(['IC ' num2str(comp)]);
        xlabel('Frequency (Hz)');
        ylabel('Normalised Power (dB)');
        xlim([2 30]); % Limit the x-axis to 2-30 Hz
    end
    sgtitle('Normalised Independent Component Power Spectra (Probe Epoch)');
else
    disp('Let us do stats!')
    % Perform statistical analysis to compare the spectra between subjects
    % Reshape data for ANOVA (flattening)
    all_probe_spectra_flat = reshape(all_probe_spectra, size(all_probe_spectra, 1), size(all_probe_spectra, 2), size(all_probe_spectra, 3) * size(all_probe_spectra, 4));

    % Define the frequency range and components for analysis
    freq_range = 2:30; % 2-30 Hz
    numComponents = size(all_probe_spectra_flat, 1);

    % Initialise arrays to hold ANOVA results
    anova_results = cell(numComponents, 1);

    % Perform ANOVA for each component across the frequency range
    for comp = 1:numComponents
        data = squeeze(all_probe_spectra_flat(comp, :, :)); % Data for the current component
        data = data'; % Transpose to have frequencies in columns and subjects in rows
        p_values = zeros(size(data, 2), 1); % Initialise p-values array

        for freq = 1:size(data, 2)
            % Perform one-way ANOVA for the current frequency across subjects
            p_values(freq) = anova1(data(:, freq), all_subject_ids, 'off');
        end

        % Store ANOVA results
        anova_results{comp} = p_values;
    end

    % Visualise the ANOVA results
    figure;
    for comp = 1:numComponents
        subplot(ceil(sqrt(numComponents)), ceil(sqrt(numComponents)), comp);
        plot(freq_range, anova_results{comp});
        title(['IC ' num2str(comp) ' ANOVA p-values']);
        xlabel('Frequency (Hz)');
        ylabel('p-value');
        xlim([2 30]); % Limit the x-axis to 2-30 Hz
        ylim([0 0.05]); % p-value threshold for significance
    end
    sgtitle('ANOVA Results for Independent Component Power Spectra');

    % Identify significant frequencies (e.g., p < 0.05)
    significant_freqs = cell(numComponents, 1);
    for comp = 1:numComponents
        significant_freqs{comp} = find(anova_results{comp} < 0.05);
    end

    disp('ANOVA analysis complete. Significant frequencies identified.');
end


disp('End of Component power spectra script')
