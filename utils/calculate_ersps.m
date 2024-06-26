%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Requires Signal Processing Toolbox
% Could use EEGLAB's newtimef instead of the spectrogram function - will
% run slower.
%%
function EEG = calculate_ersps(EEG, EEG_fixation, baseline_window, ersp_window, ersp_step, num_bootstraps, visualize)
    % Function to calculate ERSPs for each component and assess statistical significance
    % Inputs:
    %   EEG - EEGLAB EEG structure with ICA activations
    %   EEG_fixation - EEGLAB EEG structure of the fixation period for baseline
    %   baseline_window - vector specifying the start and end of the baseline period [start, end] in seconds
    %   ersp_window - window length for the FFT (in samples)
    %   ersp_step - step size for the moving window (in samples)
    %   num_bootstraps - number of bootstrap resamples for significance testing
    %   visualize - boolean flag to indicate whether to visualize the results
    % Outputs:
    %   EEG - EEGLAB EEG structure with added fields for ERSPs and significance

    disp('The calculate_ersps function requires Signal Processing Toolbox add-on.')

    % Number of components and trials
    num_components = size(EEG.icaact, 1);
    num_trials = size(EEG.data, 3);
    srate = EEG.srate;

    % Debugging: Print the sizes of the arrays
    disp(['Number of components: ' num2str(num_components)]);
    disp(['Number of trials: ' num2str(num_trials)]);
    disp(['Sample rate: ' num2str(srate)]);

    % Step 1: Preprocessing
    % Pre-calculate the baseline power from the fixation period
    [s, f, t] = spectrogram(EEG_fixation.icaact(1, :), ersp_window, ersp_step, [], srate); % Example spectrogram to get dimensions
    baseline_power = zeros(num_components, length(f));
    
    for comp = 1:num_components
        % Transform to Spectrographic Image using spectrogram
        [s, f, t] = spectrogram(EEG_fixation.icaact(comp, :), ersp_window, ersp_step, [], srate);
        % Calculate log power
        log_power = log(abs(s) .^ 2);
        % Step 2: Baseline Normalization
        % Compute the mean baseline spectrum during the central 4 seconds of the fixation period
        baseline_indices = find(t >= baseline_window(1) & t <= baseline_window(2));
        if ~isempty(baseline_indices)
            baseline_power(comp, :) = mean(log_power(:, baseline_indices), 2);
        else
            error('Baseline indices are empty. Check baseline window.');
        end
    end

    % Pre-allocate storage for ERSPs and significance masks
    EEG.icaact_ersp = cell(num_components, 1);
    EEG.icaact_ersp_sig = cell(num_components, 1);

    % Initialize time and frequency variables
    ERSP_times = [];
    ERSP_freqs = [];
    
    % Loop through each component
    for comp = 1:num_components
        ersp_all_trials = [];
        %disp(['Component: ' num2str(comp)])
        % Loop through each trial
        for trial = 1:num_trials
            %disp(['Trial: ' num2str(trial)])
            data = EEG.icaact(comp, :, trial);
            % Step 1: Preprocessing
            % Transform to Spectrographic Image using spectrogram
            [s, f, t] = spectrogram(data, ersp_window, ersp_step, [], srate);
            if isempty(ERSP_times)
                ERSP_times = t;
                ERSP_freqs = f;
            end
            % Calculate log power
            log_power = log(abs(s) .^ 2);
            % Step 3: ERSP Calculation
            % Baseline Subtraction: Subtract the trial-mean log baseline power from the log power
            ersp = bsxfun(@minus, log_power, baseline_power(comp, :)');
            ersp_all_trials = cat(3, ersp_all_trials, ersp);
        end
        
        % Step 5: Averaging
        % Compute the mean ERSP across trials
        mean_ersp = mean(ersp_all_trials, 3);

        % Debugging: Check statistics of the mean_ersp
        disp(['Component ' num2str(comp) ' mean ERSP stats:']);
        disp(['Min: ' num2str(min(mean_ersp(:)))]);
        disp(['Max: ' num2str(max(mean_ersp(:)))]);
        disp(['Mean: ' num2str(mean(mean_ersp(:)))]);
        disp(['Std Dev: ' num2str(std(mean_ersp(:)))]);

        EEG.icaact_ersp{comp} = mean_ersp;

        if num_trials > 1
            % Step 6: Statistical Significance
            % Bootstrap resampling to estimate baseline variability
            bootstrap_distributions = zeros(num_bootstraps, length(f), length(t));
            for b = 1:num_bootstraps
                resample_indices = randi([1 size(baseline_power, 2)], [1 size(baseline_power, 2)]);
                resampled_baseline_power = baseline_power(comp, resample_indices);
                bootstrap_distributions(b, :, :) = mean(log(resampled_baseline_power), 2);
            end

            % Calculate threshold for significance
            bootstrap_means = mean(bootstrap_distributions, 1);
            bootstrap_stds = std(bootstrap_distributions, 0, 1);
            significance_threshold = squeeze(bootstrap_means + 2.58 * bootstrap_stds); % 99% confidence interval

            % Identify significant changes in power
            sig_mask = mean_ersp > significance_threshold;
            EEG.icaact_ersp_sig{comp} = sig_mask;
        end
    end

    % Store the time and frequency information
    EEG.times = ERSP_times;
    EEG.freqs = ERSP_freqs;

    % Visualization
    if visualize
        figure;
        for comp = 1:num_components
            subplot(num_components, 1, comp);
            if ~isempty(EEG.icaact_ersp{comp})
                imagesc(EEG.times, EEG.freqs, EEG.icaact_ersp{comp});
                set(gca, 'YDir', 'normal');
                title(['Component ' num2str(comp) ' ERSP']);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                colorbar;
                if num_trials > 1 && isfield(EEG, 'icaact_ersp_sig') && ~isempty(EEG.icaact_ersp_sig{comp})
                    hold on;
                    contour(EEG.times, EEG.freqs, EEG.icaact_ersp_sig{comp}, [1 1], 'LineColor', 'k');
                end
            else
                % Debugging message
                disp(['No ERSP data for Component ' num2str(comp)]);
                text(0.5, 0.5, ['No data for Component ' num2str(comp)], 'HorizontalAlignment', 'center');
            end
        end
    end

    disp('ERSPs and significance calculated for each component');
end
