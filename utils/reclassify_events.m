%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===Function to remove incorrect trials
function EEG = reclassify_events(EEG)

    % Define trial start and end events
    trialStartEvent = 's15';  % rest - at start if the trails
    trialEndEvents = {'s88', 's89'};  % correct or incorrect response

    %{
    %%%%%%%%% JUST TO CHECK IF CHANGE TRIGGERS WORKS
    % save event data to a csv file for inspection
    % Check if EEG.event is not empty
    if ~isempty(EEG.event)
        % Initialize variables to store event details
        eventType = {EEG.event.type}';
        eventLatency = num2cell([EEG.event.latency]');
        eventDuration = num2cell([EEG.event.duration]');
        
        % Check if all events have a 'duration' field; if not, adjust accordingly
        if isfield(EEG.event, 'duration')
            % Combine into a table (if 'duration' is available for all events)
            eventsTable = table(eventType, eventLatency, eventDuration, 'VariableNames', {'Type', 'Latency', 'Duration'});
        else
            % If 'duration' is not a universal field, create a table without it
            eventsTable = table(eventType, eventLatency, 'VariableNames', {'Type', 'Latency'});
        end
        
        % Define your output CSV file name
        outputCSVFileName = 'Before_reclassify_EEGEvents.csv';
        
        % Write the table to a CSV file
        writetable(eventsTable, outputCSVFileName);
        
        disp(['Events have been exported to ' outputCSVFileName]);
    else
        disp('No events found in the EEG structure.');
    end
    %%%%%%%%%%%%%%%%
    %}
        
    % Initialize new event structure
    new_events = EEG.event;

    % Define memorize and ignore letter event types
    memoriseLetterEvents = {'s30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 's50', 's51', 's52', 's53', 's54', 's55', 's56', 's57', 's70', 's71', 's72', 's73', 's74', 's75', 's76', 's77'};
    ignoreLetterEvents = {'s40', 's41', 's42', 's43', 's44', 's45', 's46', 's47', 's60', 's61', 's62', 's63', 's64', 's65', 's66', 's67', 's80', 's81', 's82', 's83', 's84', 's85', 's86', 's87'};

    % Loop through each event
    for i = 1:length(EEG.event)
        if strcmp(EEG.event(i).type, trialStartEvent)
            trial_start_idx = i;
            
            % Find the corresponding trial end
            for j = trial_start_idx:length(EEG.event)
                if ismember(EEG.event(j).type, trialEndEvents)
                    trial_end_idx = j;
                    break;
                end
            end
            
            % Identify events within this trial
            trial_events = EEG.event(trial_start_idx:trial_end_idx);
            
            % Initialize counters for memorize and ignore letters
            memorize_count = 0;
            ignore_count = 0;
            
            % Reclassify events
            for k = 1:length(trial_events)
                event_type = trial_events(k).type;
                
                if ismember(event_type, memoriseLetterEvents)  % Memorize letters
                    memorize_count = memorize_count + 1;
                    new_event_type = sprintf('s3%d', memorize_count - 1);  % Zero-indexed position
                    trial_events(k).type = new_event_type;
                elseif ismember(event_type, ignoreLetterEvents)  % Ignore letters
                    ignore_count = ignore_count + 1;
                    new_event_type = sprintf('s4%d', ignore_count - 1);  % Zero-indexed position
                    trial_events(k).type = new_event_type;
                end
            end
            
            % Update the main event structure with the reclassified events
            new_events(trial_start_idx:trial_end_idx) = trial_events;
        end
    end

    % Update EEG.event with the new events
    EEG.event = new_events;
    %[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    %{
    %%%%%%%%% JUST TO CHECK IF CHANGE TRIGGERS WORKS
    % save event data to a csv file for inspection
    % Check if EEG.event is not empty
    if ~isempty(EEG.event)
        % Initialize variables to store event details
        eventType = {EEG.event.type}';
        eventLatency = num2cell([EEG.event.latency]');
        eventDuration = num2cell([EEG.event.duration]');
        
        % Check if all events have a 'duration' field; if not, adjust accordingly
        if isfield(EEG.event, 'duration')
            % Combine into a table (if 'duration' is available for all events)
            eventsTable = table(eventType, eventLatency, eventDuration, 'VariableNames', {'Type', 'Latency', 'Duration'});
        else
            % If 'duration' is not a universal field, create a table without it
            eventsTable = table(eventType, eventLatency, 'VariableNames', {'Type', 'Latency'});
        end
        
        % Define your output CSV file name
        outputCSVFileName = 'After_reclassify_EEGEvents.csv';
        
        % Write the table to a CSV file
        writetable(eventsTable, outputCSVFileName);
        
        disp(['Events have been exported to ' outputCSVFileName]);
    else
        disp('No events found in the EEG structure.');
    end
    %%%%%%%%%%%%%%%%
    %}

end