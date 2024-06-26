%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===Function to remove incorrect trials
function EEG_out = removeIncorrectTrials(EEG, startTrigger, correctTrigger, endTrigger)
    numEvents = length(EEG.event);
    startLatencies = NaN(numEvents, 1);  % Preallocate with NaN
    correctLatencies = NaN(numEvents, 1);% Preallocate with NaN
    endLatencies = NaN(numEvents, 1);    % Preallocate with NaN
    startCount = 0;
    correctCount = 0;
    endCount = 0;
    
    % Find the latencies of the start, correct, and end events
    for i = 1:numEvents
        if strcmp(EEG.event(i).type, startTrigger)
            startCount = startCount + 1;
            startLatencies(startCount) = EEG.event(i).latency;
        elseif strcmp(EEG.event(i).type, correctTrigger)
            correctCount = correctCount + 1;
            correctLatencies(correctCount) = EEG.event(i).latency;
        elseif strcmp(EEG.event(i).type, endTrigger)
            endCount = endCount + 1;
            endLatencies(endCount) = EEG.event(i).latency;
        end
    end
    
    startLatencies = startLatencies(1:startCount);  % Trim to actual size
    correctLatencies = correctLatencies(1:correctCount); % Trim to actual size
    endLatencies = endLatencies(1:endCount);        % Trim to actual size
    
    % Preallocate segmentsToRemove with an estimated size
    maxSegments = endCount; % Maximum possible number of segments to remove
    segmentsToRemove = NaN(maxSegments, 2);
    segmentCount = 0;
    
    % Find segments to remove
    for i = 1:endCount
        % Find the nearest preceding start event
        startLatency = startLatencies(find(startLatencies < endLatencies(i) & ...
            (startLatencies > max([correctLatencies(correctLatencies < endLatencies(i)); 0])), 1, 'last'));
        if ~isempty(startLatency)
            segmentCount = segmentCount + 1;
            segmentsToRemove(segmentCount, :) = [startLatency endLatencies(i)];
        end
    end
    
    segmentsToRemove = segmentsToRemove(1:segmentCount, :);  % Trim to actual size
    
    % Remove segments
    EEG_out = EEG;
    for i = size(segmentsToRemove, 1):-1:1
        EEG_out = pop_select(EEG_out, 'nopoint', segmentsToRemove(i, 1):segmentsToRemove(i, 2));
    end
end