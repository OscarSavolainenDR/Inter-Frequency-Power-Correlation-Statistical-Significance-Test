function [all_r_phase, all_f] = arbitrary_get_IAAFT_phase_ran(x_all,Fs,max_MCs_per_channel,thresh_for_corr_inclusion)
    all_r_phase = cell(length(x_all(1,:)),1);
    all_f = cell(length(x_all(1,:)),1);

%     clipping = 0.001;
    for syn_signal = 1:length(x_all(1,:))
        clear compressed_r_phase
        % Which syn signal
        x = x_all(:,syn_signal);
        random_seed = 0;
        size(x)

        %% Perform Continuous Wavelet Transform (CWT) on x
        hh = tic;
        [S_neural,f,coi] = cwt(x,'morse',Fs); % S: Wavelet Coefficients, f: scale vector, coi : Cone of Influence.
        CWT_freq_limits = [min(f) max(f)];
        perform_1_CWT_time = toc(hh);
        fprintf(['Time to do 1 CWT was ',num2str(perform_1_CWT_time),' s \n'])
        hh = tic;

        %% Making eveything outside COI zero
        S_neural = S_neural';
        g = length(f);
        for i = 1:length(coi)
            [ ~, ix ] = min( abs( f-coi(i) ) );
            S_neural(i,ix:g)= 0;
        end

        %% Calculate Wavelet Power Spectrum
        S_neural = abs(S_neural).^2;
        S_neural(S_neural==0)= NaN; % set zero values to NaN (including all values outside the COI).

        %% Get inter-scale correlations of S -> C (in paper it is 'r')
        % Note: Does not calculate the correlations that had fewer than
        % (thresh_for_corr_inclusion x 100)% of the data points in x in their
        % calculation. This is only to insure that enough data points are
        % included in each inter-scale correlation calculation.
        temp = isnan(S_neural); % whether a value is NaN  or not
        ind = sum(temp); % the number of NaN samples per scale
        f_nan_ratio = ind/length(x); % ratio of NaN to not-NaN points for each scale
        a = f_nan_ratio < thresh_for_corr_inclusion; % 1 = that scale has enough data points, 0 = doesn't.
        f(~a) = []; % remove frequencies from the vector that don't have enough samples in their calculation.
        S_neural(:,~a) = [];

        % Cut off all of the NaN values, makes the corrcoef calculation way, way, way easier (faster) than using the 'pairwise' option in corrcoef
        [~, min_snip_cutoff] = min(temp);
        S_neural(1:min_snip_cutoff(max(cumsum(a))),:) = [];
        [~, max_snip_cutoff] = max(isnan(S_neural));
        S_neural(max_snip_cutoff(max(cumsum(a))):end,:) = [];
        clear temp

        S_neural = single(S_neural);
        processing_time = toc(hh);
        fprintf(['Processing time was ',num2str(processing_time),' s \n'])

        %% Get neural ISCM
        hh = tic;
        r_orig = corrcoef(S_neural); % Gets inter-scale correlations, ignores NaN rows
        r_orig(1:5,1:5)
        orig_correlation_time = toc(hh);
        fprintf(['Correlation time was ',num2str(orig_correlation_time),' s \n'])

        % Save computation times
        %times.loading_time = loading_time;
        times.perform_1_CWT_time = perform_1_CWT_time;
        times.processing_time = processing_time;
        times.orig_correlation_time = orig_correlation_time;

        %% IAAFT Phase-randomisation MC simulation
        for MC_iteration = 1:max_MCs_per_channel

            if MC_iteration > max_MCs_per_channel % we've produced enough phase-ran ISCMs
              continue
            end

            %% IAAFT phase ran
            % This works by downsampling the data to nyquist_boost (here: 4) x Nyquist_frequency, and
            % performing the phase-randomisation at that sampling resolution. This
            % dramatically accelerates the process for lower-frequency scales.
            nyquist_boost = 4; % for downsampling the data during phase-ran, the higher the better resolution one gets. 1 is the minimum.
            hh = tic;
            rng(random_seed+MC_iteration);
            [S_neural,ran_times_temp,AAFT_iterations_temp] = AAFT_custom(S_neural, f, Fs_neural, nyquist_boost);
            phase_ran_time(MC_iteration,1) = toc(hh);
            fprintf(['Phase-ran time for Iteration ',num2str(MC_iteration),' was ',num2str(phase_ran_time(MC_iteration)),' s \n'])

            %% Calculate inter-scale corr matrix
            hh = tic;
            r_phase = corrcoef(S_neural); % Gets inter-scale correlations, ignores NaN rows
            %[x2] = reshape_2D_symmetric_matrix_to_1D(r_phase);
            %if isequal(length(x2),(p^2-p)/2)
            compressed_r_phase(:,:,MC_iteration) = r_phase;
            %else
             %   compressed_r_phase(:,:,MC_iteration) = zeros((p^2-p)/2,1);
            %end
            phase_ran_correlation_time(MC_iteration,1) = toc(hh);
            fprintf(['Correlation time for Iteration ',num2str(MC_iteration),' was ',num2str(phase_ran_correlation_time(MC_iteration)),' s \n'])

            times.phase_ran_correlation_time = phase_ran_correlation_time;
            times.phase_ran_time = phase_ran_time;

        end
        all_r_phase{syn_signal,1} = compressed_r_phase;
        all_f{syn_signal,1} = f;
    end
    
end