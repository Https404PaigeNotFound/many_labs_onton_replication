function ersp = calculate_single_trial_ersps(EEG)
    % Define parameters for time-frequency analysis
    freq_range = [1 50]; % 1 to 50 Hz
    time_range = [-1000 2000]; % -1000 to 2000 milliseconds
    num_freqs = 50; % 1-Hz intervals
    num_times = 200; % Adjusted number of times based on data length and sampling rate

    % Initialize the ERSP matrix with dynamic dimensions
    ersp_initialized = false;

    % Initialize the ERSP matrix to avoid output argument errors
    ersp = [];

    % Compute ERSP for each component and trial
    for comp = 1:size(EEG.icaact, 1)
        for trial = 1:size(EEG.data, 3)
            fprintf('Computing ERSP for component %d, trial %d\n', comp, trial);
            
            % Compute ERSP with adjusted parameters
            try
                [ersp_temp, ~, ~, ~] = newtimef(EEG.icaact(comp, :, trial), EEG.pnts, time_range, EEG.srate, ...
                    [3 0.5], 'freqs', freq_range, 'nfreqs', num_freqs, 'timesout', num_times, 'winsize', 256); % Reduced winsize to 256

                if ~ersp_initialized
                    [num_freqs, num_times] = size(ersp_temp);
                    ersp = zeros(size(EEG.icaact, 1), num_freqs, num_times, size(EEG.data, 3));
                    ersp_initialized = true;
                end

                ersp(comp, :, :, trial) = ersp_temp;
            catch ME
                fprintf('Error computing ERSP for component %d, trial %d: %s\n', comp, trial, ME.message);
                continue;
            end
        end
    end

    % Ensure ERSP is assigned even if no valid data was processed
    if isempty(ersp)
        ersp = zeros(size(EEG.icaact, 1), num_freqs, num_times, size(EEG.data, 3));
    end
end
