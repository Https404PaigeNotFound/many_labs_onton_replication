%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function anterior_components = get_anterior_components(EEG)
    % This function returns the indices of frontal components.
    % Customise this based on your criteria or predefine the indices.
    anterior_components = find(strcmp({EEG.chanlocs.labels}, 'Fz') | ...
                              strcmp({EEG.chanlocs.labels}, 'Fpz') | ...
                              strcmp({EEG.chanlocs.labels}, 'F1') | ...
                              strcmp({EEG.chanlocs.labels}, 'F2'));
end
