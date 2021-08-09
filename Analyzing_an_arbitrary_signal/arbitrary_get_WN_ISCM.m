function [all_r_WN, all_f_WN] = arbitrary_get_WN_ISCM(x_all,Fs,max_MCs_per_channel, thresh_for_corr_inclusion, clip_epsilon)

    all_r_WN = cell(length(x_all(1,:)),1);
    all_f_WN = cell(length(x_all(1,:)),1);
    % For x_all, columns represent individual signals, rows represent time
    % samples
   
    % We iterate through each signal, clip it, then get the WN_ISCMs
    for syn_signal = 1:length(x_all(1,:))
        
        % Clip
        x = rescale(x_all(:,syn_signal),0,1);
        signal_length =   find_matching_end(x,clip_epsilon);
        x = x(1:signal_length);

        % All to get thresholded f vector
        [~,f,~] = cwt(x,'morse',Fs);
        CWT_freq_limits = [min(f) max(f)];
        [~,f,~] = calculate_inter_scale_correlation_matrix(x, Fs, CWT_freq_limits, thresh_for_corr_inclusion);

        r_WN = zeros(length(f),length(f),max_MCs_per_channel);

        %all_f = cell(length(x_all(1,:)),1);

        for MC_iteration = 1:max_MCs_per_channel

            rnd_seed = MC_iteration;
            [WN_ISCM_individual,f2] = arbitrary_generate_WN_ISCMs(rnd_seed,signal_length, Fs, CWT_freq_limits,thresh_for_corr_inclusion);
            try
                r_WN(:,:,MC_iteration) = WN_ISCM_individual;
            catch
                f
                f2
                r_WN(:,:,MC_iteration) = WN_ISCM_individual;
            end
            %all_f(,1) = f;
            %else
             %   compressed_r_phase(:,:,MC_iteration) = zeros((p^2-p)/2,1);
            %end
            %phase_ran_correlation_time(MC_iteration,1) = toc(hh);
            %fprintf(['Correlation time for Iteration ',num2str(MC_iteration),' was ',num2str(phase_ran_correlation_time(MC_iteration)),' s \n'])
        end
        all_r_WN{syn_signal,1} = r_WN;
        all_f_WN{syn_signal,1} = f;
    end
    
end