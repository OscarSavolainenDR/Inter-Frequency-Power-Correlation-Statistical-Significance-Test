% Script that iterates through the Sabes lab sessions and channels to 
% calculate the inter-scale correlation matrices. The function 
% calculate_inter_scale_correlation_matrix generalises to all time series.
% The CWT and inter-scale correlation calculations can consume significant 
% memory (10s of GBs of RAM).

% main_calc_stat_sig_inter_scale_corr_mat;

    sessions_vector = 1;
    channels_vector = 1:96;

    Inter_Scale_Stat_Sig_Sabes_Parameters; % initiliases the parameters used throughout the work

    save_results = true;

    % Iterate through sessions
    for session = sessions_vector

        save_file = [save_inter_scale_corrs,'\',processed_sessions{session},'_neural_ISCM.mat'];

        % If the inter-cale correlation matrix for the current session and 
        % channel pair has already been stored, then skip this iteration
        if exist(save_file) == 2
            fprintf(['Session ',processed_sessions{session},' has already been done \n'])
            continue
        end

        % Iterate through channels
        for channel = channels_vector 
            %% Load data (SI_CSD_Broadband)
            % Note: this takes more time than loading the data outside the
            % channel-for loop since it's the same data. However, it allows us to
            % save some memory by clearing SI_CSD_Broadband after the
            % specified channel has been accessed.
            hh = tic;
            processed_data_filename = [processed_sessions{session},'.mat'];
            load([processed_data_folder,'\',processed_data_filename])
            loading_time = toc(hh); 
            fprintf(['Loading time was ',num2str(loading_time),' s \n'])
            x = single(SI_CSD_Broadband(:,channel));
            clear SI_CSD_Broadband

            %% Finish pre-processing, i.e. standardise length (shorten), clip and pre-whiten
            %[x, n, computation_times_temp, neural_time] = shorten_clip_and_prewhiten(x, neural_time);

            %% Calculate the inter-scale correlation matrix ('C', in paper is 'r')
            [C,f,computation_times] = calculate_inter_scale_correlation_matrix(x, Fs_neural, CWT_freq_limits, thresh_for_corr_inclusion);
            computation_times{channel,1}.loading_time = loading_time;
            computation_times{channel,1}.whitening_time = computation_times_temp.whitening_time;

            neural_results{channel,1}.inter_scale_corr_matrix = C; % inter-scale correlation matrices for each channel
            neural_results{channel,1}.frequency_vector = f; % frequency vector for each channel
            neural_results{channel,1}.neural_signal_length = n; % number of samples in the neural signal
            neural_results{channel,1}.computation_times = computation_times; % computation times, for those who are interested
            neural_results{channel,1}.session = processed_sessions{session};
            neural_results{channel,1}.channel = channel;

        end
        
        %% Compress neural ISCMs into 2D matrix, cuts the memory by half
        [stored_neural_ISCMs] = compress_neural_ISCMs(channels_vector,neural_results,p);
        for channel = channels_vector
            neural_results{channel,1} = rmfield(neural_results{channel,1},'inter_scale_correlation_matrix');
        end
        neural_parameters = neural_results;

        %% Save results
        if save_results
            save(save_file,'neural_parameters','stored_neural_ISCMs')
            saving_time = toc(hh);
            fprintf(['Saving time was ',num2str(saving_time),' s \n'])
        end
    end

end