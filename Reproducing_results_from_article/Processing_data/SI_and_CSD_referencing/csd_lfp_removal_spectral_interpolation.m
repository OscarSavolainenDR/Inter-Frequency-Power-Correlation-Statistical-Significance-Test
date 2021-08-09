function [CSD_ref_signal] = csd_lfp_removal_spectral_interpolation(session,raw_Broadband,show_figures,figures_cutoff)

    % Input:
    % session: session index, a positive scalar value
    % raw_Broadband: raw neural data
    % show_figures: boolean variable, if true: show figures, if false: don't.
    % figures_cutoff: how many time samples to use in a figure

    % Output: Broadband from central electrode (5th column in
    % raw_Broadband) minus the mean LFP from the 4 surrounding electrodes.
    % Intended to remove the "non-local" contributions to the LFP.
    
    % Initialise
    Inter_Scale_Stat_Sig_Sabes_Parameters;

    % Note: Line noise is at 60 Hz
    L = floor(floor(length(raw_Broadband)/Fs_neural*60)/60*Fs_neural); % to make sure the number of samples is equal to a integer-multiple of the signal length
    nb_channels = min(size(raw_Broadband));
%     line_FFT_xlim(raw_Broadband(1:L,sample_channel),Fs_neural,[50 9000]); % to plot original Fourier series

    % Set up Spectral Interpolation
    frq = Fs_neural*linspace(0,1,L);
    [~,freq_res_index] = min(abs(frq-freq_res_peak_search)); % how the freq_res range translates into indices
    Y = fft(raw_Broadband(1:L,sample_channel)); % FFT
    P2 = abs(Y); % FFT abs

    % Find index of peak power around expected harmonic
    [~,gg] = min(abs(frq-harmonic*60)); % find closest measured frequency to theoretical line-noise harmonic
    [~,ff] = max(P2(gg-ceil(freq_res_index/2):gg+ceil(freq_res_index/2))); % get max power in area around that frequency, presumably corresponds to the harmonic. Using sample_channel: a channel with strong harmonics, based on observation of the Fourier.
    ff = ff + gg -ceil(freq_res_index/2); % index of the harmonic
    kk = frq(ff); % frequency of the harmonic

    harmonics_vector = kk/harmonic : kk/harmonic : kk/harmonic * floor(Fs_neural/60);
    harmonics_vector_width = ones(length(harmonics_vector),1) * harmonic_width; % how many Hz to throw away from either side of the harmonic peak
    harmonics_vector_neigh = ones(length(harmonics_vector),1) * harmonic_neighbours; % from how many Hz to take the average for the interpolation. Doesn't really matter, all the data gets pre-whitened later anyway (flat Fourier spectrum).

    % Frequencies to interpolate, iterate through them (could parallelise)
    f2int = zeros(length(harmonics_vector),2); % how thick harmonic is
    f4int = zeros(length(harmonics_vector),4); % how manys fs to use in averaging the power (doesn't matter since we prewhiten later anyway)
    for i = 1:length(harmonics_vector)
        f2int(i,:) = [harmonics_vector(i)-harmonics_vector_width(i) harmonics_vector(i)+harmonics_vector_width(i)];
    end
    % Frequencies used for interpolation
    for i = 1:length(harmonics_vector)
        f4int(i,:) = [f2int(i,1)-harmonics_vector_neigh(i) f2int(i,:) f2int(i,2)+harmonics_vector_neigh(i)];
    end
    
    % Interpolate 60Hz (and harmonics)
    Y = zeros(L,nb_channels);
    for i = 1:length(harmonics_vector)/2 % (/2 since the Fourier spectra is mirrored around the middle)
        % Find indices of frequencies to interpolate
        [~,gg] = min(abs(frq-f2int(i,1))); % find closest measured frequency to theoretical line-noise harmonic (left of harmonic)
        [~,gg2] = min(abs(frq-f2int(i,2)));
        smpl2int = [gg gg2]; % range of frequencies that will be interpolated

        % Find indices of frequencies used for interpolation
        [~,gg] = min(abs(frq-f4int(i,1))); % find closest measured frequency to theoretical line-noise harmonic
        [~,gg4] = min(abs(frq-f4int(i,4)));
        smpl4int = [gg smpl2int gg4]; % range of frequencies that will be used for interpolation
        
        for channel = 1:nb_channels
            % To add random phase shifts (negative for conjugates), preserve DC
            % offset. Here we calculate random phases.
            rnd_theta = -pi + (2*pi).*rand(1,smpl2int(2)-smpl2int(1)+1);

            if i == 1 % To avoid over-writing Y each time
                Y(:,channel) = fft(raw_Broadband(1:L,channel)); % FFT
            end
            P2 = abs(Y(:,channel)); % FFT abs
            
            % Interpolation magnitude
            x = [mean(P2(smpl4int(1):smpl4int(2))) mean(P2(smpl4int(3):smpl4int(4)))]; % 2 meaned values of either side of interpolation
            coefficients = polyfit([smpl2int(1), smpl2int(2)], [x(1), x(2)], 1); % 1st order polynomial to span the interpolation "gap"
            temp = coefficients(1)*[smpl2int(1):smpl2int(2)] + coefficients(2); % a*x + b form, where x is the harmonic gap to interpolate over

            % Patch spectra with interpolation
            Y(smpl2int(1):smpl2int(2),channel) = temp.*exp(1i*rnd_theta); % interpolated magnitude with randomised phase
            Y(L-smpl2int(2):L-smpl2int(1),channel) = flip(temp.*exp(-1i*rnd_theta)); % same for negative frequencies
        end
    end
    
    % There may be individual, non-line noise harmonics that need to be
    % removed. component is a cell matrix contained in
    % "Inter_Scale_Stat_Sig_Sabes_Parameters;", that contains the details
    % of individual components interpolated over in this dataset.
    if session < numel(component)+1
        Y = filter_out_individual_component(Y, Fs_neural, component{session}, component_neigbours);
    end
    % Removve a common component at 55-65 Hz.
     Y = filter_out_individual_component(Y, Fs_neural, [55 65], component_neigbours);
    
    % Get IFFT of interpolated spectrum
    raw_Broadband = single(ifft(Y,L,'symmetric'));


    %% Lowpass the Spectrally Interplated data to get the LFPs
    data_temp_LFP = single(zeros(size(raw_Broadband)));
    counter = 1; %CAR_ref_LFP = zeros(max(size(raw_Broadband)),1);
    
    % Doing the filtering in groups because of memory limitations.
    group_division = 5;
    for group = 1:ceil(min(size(raw_Broadband))/group_division)
        if group < ceil(min(size(raw_Broadband))/group_division)
            group_vector = counter:counter+group_division-1;
            data_temp_LFP(:,group_vector) = single(lowpass(raw_Broadband(:,group_vector),LFP_cutoff,Fs_neural));
        else
            group_vector = counter:min(size(raw_Broadband));
            data_temp_LFP(:,group_vector) = single(lowpass(raw_Broadband(:,group_vector),LFP_cutoff,Fs_neural));
        end
        counter = counter + group_division;
    end
    
    %% Do the CSD operation, taking the mean of the four electrodes around the referenced one as the references
    CSD_ref_signal = zeros(size(raw_Broadband));
    for y = 1:10
        for x = 1:10

            % Get electrode IDs given xy-coordinates
            [electrode_ID_ref] = return_sabes_utah_electrode_ID([y,x]);
            [electrode_ID_ref_1] = return_sabes_utah_electrode_ID([y+1,x]);
            [electrode_ID_ref_2] = return_sabes_utah_electrode_ID([y-1,x]);
            [electrode_ID_ref_3] = return_sabes_utah_electrode_ID([y,x+1]);
            [electrode_ID_ref_4] = return_sabes_utah_electrode_ID([y,x-1]);

            if electrode_ID_ref == 0 % if we indexed ground, the channel doesn't exist
                continue
            end
            
            % Electrodes we'll use to reference the central electrode
            ref_vec = [electrode_ID_ref_1; electrode_ID_ref_2; electrode_ID_ref_3; electrode_ID_ref_4];
            not_borders = logical((ref_vec ~= -1) .* (ref_vec ~= 0)); % electrodes that exist, i.e. not over the edge of the array (happens if referenced electrode is on the edge of the array), and not ground
            
            % Average the LFPs to get the reference signal
            ref_signal = mean(data_temp_LFP(:,ref_vec(not_borders)),2);
            
            % Obtain CSD-referenced data
            CSD_ref_signal(:,electrode_ID_ref) = raw_Broadband(:,electrode_ID_ref) - ref_signal;
            
            % Plot Fourier transform if desired
%             figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
%             line_FFT_xlim(normalize(CSD_ref_signal(:,electrode_ID_ref)),Fs_neural,[30 500]);
%             title(['Session ',num2str(rec),' - Channel ',num2str(electrode_ID_ref)])

  
        end
    end

    % Some plots
    if show_figures
        
        electrode_ID_ref = 10; % random channel

        % Fourier transform
        figure, line_FFT_xlim(CSD_ref_signal(:,electrode_ID_ref),Fs_neural,[50 8500]);
        title(['Electrode ',num2str(electrode_ID_ref)])
        ylim([0 1])

        figure, plot(raw_Broadband(1:figures_cutoff,electrode_ID_ref)); title('Original Broadband')
        legend('1','2','3','4','Central electrode')
        xlabel('Time (s)')
        ylabel('Amplitude')
        ylim([-1500 1500])

        figure, plot(ref_signal(1:figures_cutoff)); title('Ref LFP signal')
        legend('1','2','3','4','Central electrode')
        xlabel('Time (s)')
        ylabel('Amplitude')
        ylim([-1500 1500])

        figure, plot(CSD_ref_signal(1:figures_cutoff,electrode_ID_ref)); title('Central Broadband minus the LFPs from surrounding electrodes')
        xlabel('Time (s)')
        ylabel('Amplitude')
        ylim([-1500 1500])

    end
end
