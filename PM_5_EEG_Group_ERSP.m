%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Group ERSP Statistics
Statistical Testing
Action: Test ERSP values at each time/frequency point for significance across subjects using binomial probability.
Purpose: To identify significant spectral changes across subjects.
Figures: Fig. 6 shows significant mean power changes across subjects.
Results Section: Mean spectral power change.
%}
%%

% Group ERSP statistics

% Set variables
clear;
pathToEEGLAB = pwd; % Sets path to EEGLAD as the current working directory

% Change to EEGLAB directory and start EEGLAB
cd(pathToEEGLAB);
eeglab;

% Add the subdirectory to the path which contain custom functions
addpath('utils');


%%

% Load ERSP results
load(fullfile(pathToEEGLAB, '4_memory_load', 'memory_load_results.mat'), 'ersp_results');

%%

% Define significance threshold
p_threshold = 0.01;

% Extract unique participant identifiers
unique_participants = fieldnames(ersp_results);

% Define memory loads
memory_loads = 0:5;

% Initialise variables for group ERSP statistics
group_ersp_memorise = cell(1, length(memory_loads));
group_ersp_ignore = cell(1, length(memory_loads));
group_ersp_diff = cell(1, length(memory_loads));
sig_ersp_memorise = cell(1, length(memory_loads));
sig_ersp_ignore = cell(1, length(memory_loads));
sig_ersp_diff = cell(1, length(memory_loads));

%%

% Perform binomial probability test for each memory load
for load = memory_loads
    % Initialize matrices to store ERSP data
    ersp_memorise_data = [];
    ersp_ignore_data = [];
    
    % Loop through participants to collect ERSP data
    for p = 1:length(unique_participants)
        participant_id = unique_participants{p};
        
        % Get the ERSP data for the current memory load
        ersp_memorise = ersp_results.(participant_id).memorise{load + 1};
        ersp_ignore = ersp_results.(participant_id).ignore{load + 1};
        
        % Convert the ERSP data to numeric if necessary
        if ~isnumeric(ersp_memorise.icaact_ersp)
            ersp_memorise.icaact_ersp = cell2mat(ersp_memorise.icaact_ersp);
        end
        if ~isnumeric(ersp_ignore.icaact_ersp)
            ersp_ignore.icaact_ersp = cell2mat(ersp_ignore.icaact_ersp);
        end
        
        % Store the ERSP data
        ersp_memorise_data = cat(4, ersp_memorise_data, ersp_memorise.icaact_ersp);
        ersp_ignore_data = cat(4, ersp_ignore_data, ersp_ignore.icaact_ersp);
    end
    
    % Calculate mean ERSP for memorise and ignore conditions
    group_ersp_memorise{load + 1} = mean(ersp_memorise_data, 4);
    group_ersp_ignore{load + 1} = mean(ersp_ignore_data, 4);
    
    % Calculate ERSP difference
    group_ersp_diff{load + 1} = group_ersp_memorise{load + 1} - group_ersp_ignore{load + 1};
    
    % Binomial probability test for significance across subjects
    N = size(ersp_memorise_data, 4); % Number of participants
    sig_ersp_memorise{load + 1} = zeros(size(group_ersp_memorise{load + 1}));
    sig_ersp_ignore{load + 1} = zeros(size(group_ersp_ignore{load + 1}));
    sig_ersp_diff{load + 1} = zeros(size(group_ersp_diff{load + 1}));
    
    for t = 1:size(group_ersp_memorise{load + 1}, 2)
        for f = 1:size(group_ersp_memorise{load + 1}, 1)
            k_memorise = sum(ersp_memorise_data(f, t, :) > 0, 4); % Number of significant memorise results
            k_ignore = sum(ersp_ignore_data(f, t, :) > 0, 4); % Number of significant ignore results
            k_diff = sum((ersp_memorise_data(f, t, :) - ersp_ignore_data(f, t, :)) > 0, 4); % Number of significant difference results
            
            % Calculate binomial probabilities
            P_memorise = binocdf(k_memorise, N, p_threshold, 'upper');
            P_ignore = binocdf(k_ignore, N, p_threshold, 'upper');
            P_diff = binocdf(k_diff, N, p_threshold, 'upper');
            
            % Identify significant time/frequency points
            if P_memorise < p_threshold
                sig_ersp_memorise{load + 1}(f, t) = 1;
            end
            if P_ignore < p_threshold
                sig_ersp_ignore{load + 1}(f, t) = 1;
            end
            if P_diff < p_threshold
                sig_ersp_diff{load + 1}(f, t) = 1;
            end
        end
    end
end

%%

% Define time and frequency ranges for visualisation
time_range = linspace(-1, 2, size(group_ersp_memorise{1}, 3));
freq_range = linspace(1, 30, size(group_ersp_memorise{1}, 2));

% Define the number of components to visualize per figure
components_per_figure = 5;

%{
% Visualise group ERSP results
for load = memory_loads
    num_components = size(group_ersp_memorise{load + 1}, 1);
    num_figures = ceil(num_components / components_per_figure);
    
    for fig_num = 1:num_figures
        figure;
        for comp = 1:components_per_figure
            comp_idx = (fig_num - 1) * components_per_figure + comp;
            if comp_idx > num_components
                break;
            end
            
            subplot(components_per_figure, 3, (comp - 1) * 3 + 1);
            imagesc(time_range, freq_range, squeeze(group_ersp_memorise{load + 1}(comp_idx, :, :)));
            set(gca, 'YDir', 'normal');
            title(['Group Memorise Load ' num2str(load) ' Component ' num2str(comp_idx)]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
            hold on;
            [T, F] = meshgrid(time_range, freq_range);
            Z = squeeze(sig_ersp_memorise{load + 1}(comp_idx, :, :));
            if all(size(T) == size(Z)) && all(size(F) == size(Z))
                contour(T, F, Z, [1 1], 'LineColor', 'k');
            end
            
            subplot(components_per_figure, 3, (comp - 1) * 3 + 2);
            imagesc(time_range, freq_range, squeeze(group_ersp_ignore{load + 1}(comp_idx, :, :)));
            set(gca, 'YDir', 'normal');
            title(['Group Ignore Load ' num2str(load) ' Component ' num2str(comp_idx)]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
            hold on;
            Z = squeeze(sig_ersp_ignore{load + 1}(comp_idx, :, :));
            if all(size(T) == size(Z)) && all(size(F) == size(Z))
                contour(T, F, Z, [1 1], 'LineColor', 'k');
            end
            
            subplot(components_per_figure, 3, (comp - 1) * 3 + 3);
            imagesc(time_range, freq_range, squeeze(group_ersp_diff{load + 1}(comp_idx, :, :)));
            set(gca, 'YDir', 'normal');
            title(['Group Difference Load ' num2str(load) ' Component ' num2str(comp_idx)]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            colorbar;
            hold on;
            Z = squeeze(sig_ersp_diff{load + 1}(comp_idx, :, :));
            if all(size(T) == size(Z)) && all(size(F) == size(Z))
                contour(T, F, Z, [1 1], 'LineColor', 'k');
            end
        end
        
        % Optionally, save the figure
        saveas(gcf, ['Group_ERSP_Load_' num2str(load) '_Figure_' num2str(fig_num) '.png']);
    end
end

disp('Group ERSP statistics and visualisation complete.');

%}
