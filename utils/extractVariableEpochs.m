%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===Function to extract variable-length epochs
function EEG_out = extractVariableEpochs(EEG, startEvents, endEvents)
    % Find event latencies
    numEvents = length(EEG.event);
    startLatencies = NaN(numEvents, 1);  % Preallocate with NaN
    endLatencies = NaN(numEvents, 1);    % Preallocate with NaN
    startCount = 0;
    endCount = 0;
    
    for i = 1:numEvents
        if ismember(EEG.event(i).type, startEvents)
            startCount = startCount + 1;
            startLatencies(startCount) = EEG.event(i).latency;
        elseif ismember(EEG.event(i).type, endEvents)
            endCount = endCount + 1;
            endLatencies(endCount) = EEG.event(i).latency;
        end
    end
    
    startLatencies = startLatencies(1:startCount);  % Trim to actual size
    endLatencies = endLatencies(1:endCount);        % Trim to actual size
    
    % Preallocate epochs array with an estimated size
    maxEpochs = min(startCount, endCount); % Maximum possible number of epochs
    epochs = NaN(maxEpochs, 2);
    epochCount = 0;
    
    for i = 1:startCount
        % Find the closest end event after the start event
        endLatency = endLatencies(find(endLatencies > startLatencies(i), 1));
        if ~isempty(endLatency)
            epochCount = epochCount + 1;
            epochs(epochCount, :) = [startLatencies(i) endLatency];
        end
    end
    
    epochs = epochs(1:epochCount, :);  % Trim to actual size
    
    % Extract epochs
    EEG_out = pop_select(EEG, 'point', epochs);
    %EEG_out = pop_select(EEG, 'point', epochs(:, 1), epochs(:, 2));
end