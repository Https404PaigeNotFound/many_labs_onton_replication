%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EEG_filtered = filter_epochs_by_load(EEG, load)
    % Function to filter EEG epochs based on the specified memory load
    % Inputs:
    %   EEG - EEGLAB dataset structure containing the epochs
    %   load - The specified memory load (0-5)
    % Output:
    %   EEG_filtered - EEGLAB dataset structure containing the filtered epochs

    % Define the memorise and ignore event types with memory loads
    memoriseLetterEvents = {'s30', 's31', 's32', 's33', 's34', 's35', ...
                            's50', 's51', 's52', 's53', 's54', 's55', ...
                            's70', 's71', 's72', 's73', 's74', 's75'};
                        
    ignoreLetterEvents = {'s40', 's41', 's42', 's43', 's44', 's45', ...
                          's60', 's61', 's62', 's63', 's64', 's65', ...
                          's80', 's81', 's82', 's83', 's84', 's85'};
    
    % Convert the specified load to string for matching
    load_str = num2str(load);

    % Find the indices of the epochs that match the specified memory load
    event_indices = [];
    for i = 1:length(EEG.event)
        event_type = EEG.event(i).type;
        if any(strcmp(event_type, memoriseLetterEvents)) || any(strcmp(event_type, ignoreLetterEvents))
            % Check if the last character matches the specified load
            if strcmp(event_type(end), load_str)
                event_indices = [event_indices, EEG.event(i).epoch];
            end
        end
    end

    % Keep only the unique epoch indices
    event_indices = unique(event_indices);

    % Select the epochs that match the specified memory load
    EEG_filtered = pop_select(EEG, 'trial', event_indices);
end
