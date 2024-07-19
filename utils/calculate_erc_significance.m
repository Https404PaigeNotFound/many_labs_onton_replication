%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to calculate ERC significance using bootstrap
function significant_erc = calculate_erc_significance(erc_matrix, baseline_data, wavelet_params, num_bootstraps, p_threshold)
    num_components = size(erc_matrix, 1);
    num_frequencies = size(erc_matrix, 3);
    num_timepoints = size(erc_matrix, 4);
    
    % Initialize significance matrix
    significant_erc = zeros(size(erc_matrix));
    
    % Bootstrap method
    for b = 1:num_bootstraps
        % Resample baseline data
        resampled_data = baseline_data(:, randi(size(baseline_data, 2), 1, size(baseline_data, 2)), :);
        
        % Calculate phase coherence for resampled data
        erc_baseline = zeros(num_components, num_components, num_frequencies, num_timepoints);
        for comp1 = 1:num_components
            for comp2 = comp1+1:num_components
                [coherence, ~] = compute_phase_coherence(squeeze(resampled_data(comp1, :, :)), squeeze(resampled_data(comp2, :, :)), wavelet_params);
                erc_baseline(comp1, comp2, :, :) = coherence;
                erc_baseline(comp2, comp1, :, :) = coherence; % Symmetric matrix
            end
        end
        
        % Compare with actual ERC values
        significant_erc = significant_erc + (erc_matrix > erc_baseline);
    end
    
    % Calculate p-values and determine significance
    significant_erc = significant_erc / num_bootstraps < p_threshold;
end