%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGLAB 2024.0 | MATLAB R2024a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ICs_of_interest = identify_components_of_interest(EEG)

    % Initialize an empty array to hold indices of ICs of interest
    ICs_of_interest = [];

    % Set criteria for theta band (e.g., 4-8 Hz) and a dipole location criterion
    theta_range = [4 8];
    acc_location = [0 40 40]; % Approximate MNI coordinates for ACC
    max_dist = 15; % Maximum distance from ACC in mm

    % Loop through each independent component
    for ic = 1:size(EEG.icaweights, 1)
        % Get the power spectrum of the IC
        [spectra, freqs] = spectopo(EEG.icaact(ic, :, :), 0, EEG.srate);

        % Check if there's a peak in the theta band
        theta_power = mean(spectra(freqs >= theta_range(1) & freqs <= theta_range(2)));

        % Perform dipole fitting and check the location
        EEG = pop_dipfit_settings(EEG, 'hdmfile', 'standard_BEM/standard_vol.mat', ...
            'coordformat', 'MNI', 'mrifile', 'standard_mri.mat', ...
            'chanfile', 'standard_1005.elc', 'coord_transform', [0 0 0 0 0 0 1 1 1], ...
            'chansel', 1:EEG.nbchan);
        
        EEG = pop_multifit(EEG, ic, 'dipoles', 1, 'rmse', 10, 'threshold', 100);
        dipole_loc = EEG.dipfit.model(ic).posxyz;

        % Calculate distance from ACC
        dist_from_acc = sqrt(sum((dipole_loc - acc_location).^2));
        
        % Check if the IC meets the criteria
        if theta_power > threshold && dist_from_acc < max_dist
            ICs_of_interest = [ICs_of_interest, ic];
        end
    end
end
