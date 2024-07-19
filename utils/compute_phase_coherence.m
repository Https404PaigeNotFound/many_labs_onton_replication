%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to compute phase coherence using wavelet transform
function [coherence, phase_diff] = compute_phase_coherence(data1, data2, wavelet_params)
    frequencies = wavelet_params.frequencies;
    num_cycles = wavelet_params.num_cycles;
    num_frequencies = length(frequencies);
    num_timepoints = size(data1, 2);
    num_trials = size(data1, 1);
    
    coherence = zeros(num_frequencies, num_timepoints);
    phase_diff = zeros(num_frequencies, num_timepoints, num_trials);
    
    for f = 1:num_frequencies
        freq = frequencies(f);
        cycle = num_cycles(f);
        
        % Compute wavelet transform for each trial
        wavelet_data1 = wavelet_transform(data1, freq, cycle);
        wavelet_data2 = wavelet_transform(data2, freq, cycle);
        
        % Calculate phase difference
        phase_diff(f, :, :) = angle(wavelet_data1 ./ wavelet_data2);
        
        % Calculate phase coherence
        coherence(f, :) = abs(mean(exp(1i * phase_diff(f, :, :)), 3));
    end
end
