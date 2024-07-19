%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to visualize significant ERC results
function visualise_erc_results(significant_erc, wavelet_params, condition_names, participant_id)
    frequencies = wavelet_params.frequencies;
    num_conditions = length(condition_names);
    
    % Define time vector (assuming a specific sampling rate and length)
    time_vector = linspace(-1, 2, size(significant_erc{1}, 4));  % Adjust the range as necessary
    
    for cond_idx = 1:num_conditions
        condition = condition_names{cond_idx};
        sig_erc = significant_erc{cond_idx};
        
        % Create figure for each condition
        figure;
        num_components = size(sig_erc, 1);
        
        for comp1 = 1:num_components
            for comp2 = comp1+1:num_components
                subplot(num_components, num_components, (comp1-1)*num_components + comp2);
                imagesc(time_vector, frequencies, squeeze(sig_erc(comp1, comp2, :, :)));
                set(gca, 'YDir', 'normal');
                title(['Significant ERC: ' condition ' (' num2str(comp1) '-' num2str(comp2) ')']);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                colorbar;
            end
        end
        
        % Save the figure
        saveas(gcf, ['ERC_Significant_' condition '_' participant_id '.png']);
    end
end
