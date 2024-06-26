%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===Function to compute ICA activation
function EEG = compute_ica_activations(EEG)
% Compute independent component activations from an EEGLAB dataset
% Input:
%   EEG: EEGLAB dataset structure containing the data and ICA matrices
% Output:
%   EEG: Updated EEGLAB dataset structure with the computed icaact field

    % Check if the required fields exist
    if ~isfield(EEG, 'data') || ~isfield(EEG, 'icaweights') || ~isfield(EEG, 'icasphere')
        error('EEG structure must contain data, icaweights, and icasphere fields.');
    end

    % Get the channel data
    Y = EEG.data;
    [~, n_timepoints, n_epochs] = size(Y);
    n_components = size(EEG.icaweights, 1);
    
    % Initialize the icaact field
    EEG.icaact = zeros(n_components, n_timepoints, n_epochs);
    
    % Compute the activations for each epoch
    for epoch_i = 1:n_epochs
        Y_epoch = squeeze(Y(:,:,epoch_i)); % Extract the data for the current epoch
        
        % Compute the ICA activations
        EEG.icaact(:,:,epoch_i) = EEG.icaweights * EEG.icasphere * Y_epoch;
        
        % Debugging information
        if epoch_i == 1
            disp('First epoch ICA activations computed:');
            %disp(EEG.icaact(:,:,epoch_i));
        end
    end
    
    % Final check
    if ~isempty(EEG.icaact)
        disp('ICA activations computed and populated in EEG.icaact');
    else
        disp('ICA activations computation failed.');
    end
end
