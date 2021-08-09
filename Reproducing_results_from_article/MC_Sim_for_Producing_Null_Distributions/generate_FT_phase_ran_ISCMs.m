
function generate_FT_phase_ran_ISCMs(PBS_array_index, cores, max_MCs_per_channel, MC_index_offset, nb_channels)

    % Calculates the ISCMs for the channels in a specified session-channel
    % pair, given by the indexing below.
    % Phase randomizes each WPS frequency band so we can test against
    % autocorr in the final test. Uses FT phase-randomisation. The
    % phase-randomisation-derived ISCMs are compressed, taking advantage of
    % the symmetrical nature of the correlation matrices.
    
    Inter_Scale_Stat_Sig_Sabes_Parameters;
    
    %% Indexing
    channel = rem(PBS_array_index,nb_channels);
    if channel == 0
        channel = nb_channels;
    end
    session = ceil(PBS_array_index/nb_channels);
    random_seed = (PBS_array_index-1) * 1000 + MC_index_offset;
    
    session_channel_seed = [session channel random_seed]


    %% Set path to results
    % HPC save path
    save_path = [save_FT_phase_ran_ISCMs,'/Session_',num2str(session)];
    mkdir(save_path)
    if exist([save_path,'/Channel_',num2str(channel),'.mat']) == 2
        fprintf(['Session ',num2str(session),'; Channel ',num2str(channel),' done \n'])
        try
            load([save_path,'/Channel_',num2str(channel),'.mat'])
            MC_local_offset = length(compressed_r_phase(1,:)); % how many phase-ran ISCMs have been produced
        catch
            MC_local_offset = 0; % file does not exist/ is corrupt, i.e. no phase-ran ISCMs have been produced
        end
    else
        MC_local_offset = 0;
    end
    MC_local_offset
    
    % If we have produced enough phase-ran ISCMs
    if MC_local_offset > max_MCs_per_channel-1
      return
    end
    
    %% Set paths and load SI-CSD processed neural data
    hh = tic;
    load([processed_data_folder,'/',processed_sessions{session},'.mat'])
    loading_time = toc(hh); 
    fprintf(['Loading time was ',num2str(loading_time),' s \n'])
    
    x = single(SI_CSD_Broadband(:,channel));
    size(x)
    x_sample_orig = x(1:10000);  

    %% Standardize and clip the data length
    neural_time = zeros(length(x),1);
    [x, n, computation_times_temp, neural_time] = shorten_and_clip(x, neural_time);
    size(x)
    clipping_counter = 1;
    
    % If the entire signal was clipped, re-do after having removed some
    % data from the beginning of the signal.
    while length(x) == 1
      fprintf('Re-doing the clipping \n')
      clipping_counter = clipping_counter + 1e4;
      x = single(SI_CSD_Broadband(clipping_counter:end,channel));
      neural_time = zeros(length(x),1);
      [x, n, computation_times_temp, neural_time] = shorten_and_clip(x, neural_time);
    end
    size(x)
    clear SI_CSD_Broadband neural_time
%     x_sample_pre_processed = x(1:10000);

    %% Perform Continuous Wavelet Transform (CWT) on x
    hh = tic;
    [S_neural,f,coi] = cwt(x,'morse',Fs_neural,'FrequencyLimits',CWT_freq_limits); % S: Wavelet Coefficients, f: scale vector, coi : Cone of Influence.
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
    times.loading_time = loading_time;
    times.perform_1_CWT_time = perform_1_CWT_time;
    times.processing_time = processing_time;
    times.orig_correlation_time = orig_correlation_time;
    
    %% IAAFT Phase-randomisation MC simulation
    for MC_iteration = 1+MC_local_offset:max_MCs_per_channel+MC_local_offset
    
        if MC_iteration > max_MCs_per_channel % we've produced enough phase-ran ISCMs
          continue
        end
        
        %% Phase-randomize neural scales, using FT phase-ran
        hh = tic;
        rng(random_seed+MC_iteration);
        S_neural = phaseRandomize_2(S_neural);
        phase_ran_time(MC_iteration,1) = toc(hh);
        fprintf(['Phase-ran time for Iteration ',num2str(MC_iteration),' was ',num2str(phase_ran_time(MC_iteration)),' s \n'])

        %% Calculate inter-scale corr matrix
        hh = tic;
        r_phase = corrcoef(S_neural); % Gets inter-scale correlations, ignores NaN rows
        [x2] = reshape_2D_symmetric_matrix_to_1D(r_phase);
        if isequal(length(x2),(p^2-p)/2)
            compressed_r_phase(:,MC_iteration) = x2;
        else
            compressed_r_phase(:,MC_iteration) = zeros((p^2-p)/2,1);
        end
        phase_ran_correlation_time(MC_iteration,1) = toc(hh);
        fprintf(['Correlation time for Iteration ',num2str(MC_iteration),' was ',num2str(phase_ran_correlation_time(MC_iteration)),' s \n'])
       
        times.phase_ran_correlation_time = phase_ran_correlation_time;
        times.phase_ran_time = phase_ran_time;
        if rem(MC_iteration,5) == 0
            % save([save_path,'/Channel_',num2str(channel),'.mat'],'n','clipping_counter','x_sample_pre_processed','x_sample_orig','times','compressed_r_phase','r_orig','f','channel','session','cores','PBS_array_index','max_MCs_per_channel','random_seed');
            save([save_path,'/Channel_',num2str(channel),'.mat'],'n','x_sample_orig','times','compressed_r_phase','r_orig','f','channel','session','cores','PBS_array_index','max_MCs_per_channel','random_seed');
        end
    end

    % 'x_sample_pre_processed','clipping_counter',
    hh = tic;
    save([save_path,'/Channel_',num2str(channel),'.mat'],'n','x_sample_orig','times','compressed_r_phase','r_orig','f','channel','session','cores','PBS_array_index','max_MCs_per_channel','random_seed');
    saving_time = toc(hh);
    fprintf(['Saving time was ',num2str(saving_time),' s \n'])

end