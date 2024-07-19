%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to perform wavelet transform
function wavelet_data = wavelet_transform(data, frequency, num_cycles)
    % Create wavelet
    time = -1:1/1000:1; % Define time vector (e.g., -1 to 1 second)
    wavelet = (exp(2 * 1i * pi * frequency * time) .* exp(-time .^ 2 / (2 * (num_cycles / (2 * pi * frequency)) ^ 2)));
    wavelet = wavelet / sum(abs(wavelet)); % Normalize wavelet
    
    % Initialize wavelet data
    num_trials = size(data, 1);
    num_timepoints = size(data, 2);
    wavelet_data = zeros(num_trials, num_timepoints);
    
    % Convolve wavelet with data
    for trial = 1:num_trials
        wavelet_data(trial, :) = conv(squeeze(data(trial, :)), wavelet, 'same');
    end
end