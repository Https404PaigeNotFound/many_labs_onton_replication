%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===Function to remove incorrect trials
function EEG = reclassify_events(EEG)


    %{
    %%%%%%%%% JUST TO CHECK IF CHANGE TRIGGERS WORKS
    % save event data to a csv file for inspection
    % Check if EEG.event is not empty
    if ~isempty(EEG.event)
        % Initialize variables to store event details
        eventType = {EEG.event.type}';
        eventLatency = num2cell([EEG.event.latency]');
        eventDuration = num2cell([EEG.event.duration]');
        eventCount = num2cell([EEG.event.urevent]');
        
        % Check if all events have a 'duration' field; if not, adjust accordingly
        if isfield(EEG.event, 'duration')
            % Combine into a table (if 'duration' is available for all events)
            eventsTable = table(eventType, eventLatency, eventDuration, eventCount, 'VariableNames', {'Type', 'Latency', 'Duration', 'Count'});
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

    
    % Define Trigger Events 
    blockRestEvent = {'s19', 's29'}; % block_start % s19 = short, s29 = long 
    interTrialRestEvent = {'s15'}; % rest
    fixationCrossEvent = {'s10'}; % fixation_cross 
    memoriseLetterEvents = {'s30', 's31', 's32', 's33', 's34', 's35', 's36', 's37', 's50', 's51', 's52', 's53', 's54', 's55', 's56', 's57', 's70', 's71', 's72', 's73', 's74', 's75', 's76', 's77'}; % trail_3_letters  
    % The first number 3,5 or 7 indicates memory load and the following number is letter position (zero-indexed). 
    ignoreLetterEvents = {'s40', 's41', 's42', 's43', 's44', 's45', 's46', 's47', 's60', 's61', 's62', 's63', 's64', 's65', 's66', 's67', 's80', 's81', 's82', 's83', 's84', 's85', 's86', 's87'}; % trail_3_letters  
    % The first number 4,6 or 8 indicates (memory load+1) just so can differentiate between the memorise and ignore. The following number is letter position (zero-indexed). 
    variableMaintenancePeriodEvent = {'s11'}; % memorise 
    probeEvents = {'s38', 's39', 's58', 's59', 's78', 's79'}; % probe % s38 means 3 memory load and probe letter was in set and s39 is same memory load but probe letter NOT in set  
    correctAnswerEvent = {'s88'}; % response_validation  
    wrongAnswerEvent = {'s89'}; % response_validation 
    blockEndEvent = {'s91', 's92'}; % blockEnd % s91 = short, s92 = long 
    endEvent = {'s99'}; % End 

    % Initialise new event fields
    for i = 1:length(EEG.event)
        EEG.event(i).AccLoad = 'n/a';
        EEG.event(i).MainLoad = 'n/a';
    end

    % Init counters
    accumulated_count = 0;  % Tracks number of memorise letters shown
    current_main_load = 'n/a';  % Tracks MainLoad for the trial
    last_memorise_count = 0;  % Tracks latest memorise count
    last_event = 'n/a';

    for i = 1:length(EEG.event)
        eventType = EEG.event(i).type;

        % Reset at Fixation or Rest (New Trial)
        if ismember(eventType, fixationCrossEvent) || ismember(eventType, interTrialRestEvent)
            EEG.event(i).AccLoad = 'n/a';
            EEG.event(i).MainLoad = 'n/a';
            accumulated_count = 0;
            last_memorise_count = 0;
            current_main_load = 'n/a';
            last_event = EEG.event(i).AccLoad;

        % Memorise Letter Events - Increments Load
        elseif ismember(eventType, memoriseLetterEvents)
            EEG.event(i).AccLoad = sprintf('L%d_Memorise', accumulated_count);
            last_event = EEG.event(i).AccLoad;
            mainLoadValue = eventType(2); % Extract first digit (memory load)
            current_main_load = sprintf('L%s_Maintenance', mainLoadValue);
            EEG.event(i).MainLoad = current_main_load; % Add Maintenance Load
            last_memorise_count = accumulated_count; % Update memorise count
            accumulated_count = accumulated_count + 1; % Increment memorise count

        % Ignore Letter Events 
        elseif ismember(eventType, ignoreLetterEvents)
            mainLoadValue = regexp(eventType, '\d', 'match'); % Extract all digits
            mainLoadValue = mainLoadValue{1}; % Get the first extracted digit
            mainLoadValue = str2double(mainLoadValue) - 1; % Convert to integer and subtract 1
            current_main_load = sprintf('L%d_Maintenance', mainLoadValue); 
            EEG.event(i).MainLoad = current_main_load; % Add Maintenance Load
            % If no memorise letters have been shown yet, use L0_Ignore and what to do if the last_event was L0_Memorise
            if strcmp(last_event, 'L0_Memorise')
                EEG.event(i).AccLoad = 'L1_Ignore';
            elseif last_memorise_count == 0
                EEG.event(i).AccLoad = 'L0_Ignore';
                last_event = EEG.event(i).AccLoad;
            else
                EEG.event(i).AccLoad = sprintf('L%d_Ignore', accumulated_count);
                last_event = EEG.event(i).AccLoad;
            end

        % Maintenance Event - Uses Current Main Load
        elseif ismember(eventType, variableMaintenancePeriodEvent)
            EEG.event(i).AccLoad = 'n/a';
            EEG.event(i).MainLoad = current_main_load;
            last_event = EEG.event(i).AccLoad;

        % Probe, Response, and Rest - Set to 'n/a'
        elseif ismember(eventType, probeEvents) || ismember(eventType, correctAnswerEvent) || ismember(eventType, wrongAnswerEvent) || ismember(eventType, interTrialRestEvent)
            EEG.event(i).AccLoad = 'n/a';
            EEG.event(i).MainLoad = 'n/a';
            last_event = EEG.event(i).AccLoad;
        end
    end

    %{
    %%%%%%%%% JUST TO CHECK IF CHANGE TRIGGERS WORKS
    % save event data to a csv file for inspection
    % Check if EEG.event is not empty
    if ~isempty(EEG.event)
        % Initialize variables to store event details
        eventType = {EEG.event.type}';
        eventLatency = num2cell([EEG.event.latency]');
        eventDuration = num2cell([EEG.event.duration]');
        eventAccLoad = {EEG.event.AccLoad}';
        eventMainLoad = {EEG.event.MainLoad}';

        
        
        % Check if all events have a 'duration' field; if not, adjust accordingly
        if isfield(EEG.event, 'duration')
            % Combine into a table (if 'duration' is available for all events)
            eventsTable = table(eventType, eventLatency, eventDuration, eventCount, eventAccLoad, eventMainLoad, 'VariableNames', {'Type', 'Latency', 'Duration', 'Count', 'AccLoad', 'MainLoad'});
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