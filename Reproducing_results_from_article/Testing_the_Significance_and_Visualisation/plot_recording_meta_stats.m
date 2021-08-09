
% Function for plotting meta-statistical information about the statistically significant inter-
% frequency correlation results across all analysed recordings.

function plot_recording_meta_stats(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters)

    Inter_Scale_Stat_Sig_Sabes_Parameters; % initialise parameters used in this work

    % Parameters for this script
    show_figures = meta_parameters.show_figures; 
    try_load_C_flag = meta_parameters.try_load_C_flag; % if it fails to load each specified neural ISCM, it continues if this flag is true. If false, it stops when it fails. False is better for testing, true is better for skipping over missing results.
    try_plot_flag = meta_parameters.try_plot_flag; % if it fails to plot the ISCM, then if 'true', it continues. If 'false', it doesn't.
    local_colorbar_bool = meta_parameters.local_colorbar_bool; % whether to plot the colorbar individually for each ISCM
    plot_abs_corr = meta_parameters.plot_abs_corr; % whether to plot the correlations (false) or the absolute values of the correlations (true). The absolute values are less informative, but make the color scheme more intuitive.
    colormap_scheme = meta_parameters.colormap_scheme; % what MATLAB colormap scheme to use
    second_colormap_scheme = meta_parameters.second_colormap_scheme; % what MATLAB colormap scheme to use for meta-percentage plot
    
    % Initialise parameters
    H_stored = zeros(p);
    all_stored_c_temp = zeros(p,p,length(sessions_vector)*length(channels_vector));
    pair_counter = 0;

    % Choose which phase-ran ISCMs are used
    if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'FT') % use FT phase-ran nulls
        save_phase_ran_ISCMs = save_FT_phase_ran_ISCMs;
    elseif strcmp(use_IAAFT_or_FT_phase_ran_nulls,'IAAFT') % use IAAFT phase-ran nulls
        save_phase_ran_ISCMs = save_IAAFT_phase_ran_ISCMs;
    else
        error('use_IAAFT_or_FT_phase_ran_nulls should be a string equal to FT or IAAFT')
    end
    
    
    %% Load white noise ISCMs
    length_MC_counter = 0;
    for length_MC = length_MC_vector % 350, 400 or 500
        length_MC_counter = length_MC_counter + 1;
        load([save_white_noise_ISCMs,'\All_White_Noise_MCs_Length_',num2str(length_MC),'.mat'])
        for mc_iteration = 1:length(WN_ISCM(1,:))
            C_noises(:,:,mc_iteration) = reshape_1D_vector_to_2D_symmetric_matrix(WN_ISCM(:,mc_iteration));
        end
        white_noise_ISCM_all{length_MC_counter} = C_noises; % all recording lengths  
    end
    clear C_noises WN_ISCM




    %% Iterate through listed channels
    % Plots the (significant values of the) inter-scale correlation matrices 
    % of each session
    chan_counter = 0;
    for channel = channels_vector
        chan_counter = chan_counter + 1;
        session_counter = 0;
        for session = sessions_vector % iterate through the listed sessions in each channel
            session_counter = session_counter + 1;
            pair_counter = pair_counter + 1;
            if rem(pair_counter,10) == 0
                fprintf(['Session-channel pair ',num2str(pair_counter),' out of ',num2str(length(channels_vector)*length(sessions_vector)),'\n'])
            end
            counter = session_counter;

            if try_plot_flag
                try

                    % Load phase-ran ISCMs (and neural ISCM) 
                    [r_neural,phase_ran_ISCMs,random_seed_vector,f_return,n] = load_phase_randomised_null(session,channel, save_phase_ran_ISCMs, try_load_C_flag,use_IAAFT_or_FT_phase_ran_nulls);
                    if isempty(f_return) % if the phase-ran null does not exist for the session-channel pair
                        phase_ran_ISCMs = zeros(p,p,2);
                        r_neural = [];
                        nb_MCs_phase_ran(session,channel) = 0;
                        if isempty(f) % if f has never been initialised before
                            f = zeros(length(p),1);
                        end               
                    else
                        f = f_return; % the vector of frequencies in the CWT
                        nb_MCs_phase_ran(session,channel) = length(phase_ran_ISCMs(1,1,:));
                    end
 
                    %% Load white noise ISCMs                    
                    % Check recording length, see which standardised lengh
                    % it corresponds to
                    session_lengths = get_standardised_length(length_MC_vector, minimum_time_gap, round(n/Fs_neural));
                    length_MC_counter = find(session_lengths == length_MC_vector); % gives the index for what part of the C_noises_all matrix to index for the given session length
             
                    % Check the standardisation of length went correcty, if
                    % not, skip this session-channel pair
                    if isempty(length_MC_counter)
                        continue
                    end

                    % Index appropriate white-noise ISCMs
                    white_noise_ISCMs = white_noise_ISCM_all{length_MC_counter}; % index correct MC length
                    nb_MCs = length(white_noise_ISCMs(1,1,:));
                    
                    %% Get null distribution
                    % Sum mean white-noise ISCM and phase-ra nISCM
                    % distributions to get null distribution
                    null_distributions = phase_ran_ISCMs + mean(white_noise_ISCMs,3) - eye(p); % deleting the identity matrix is optional, since diagonal elements are never considered in the stat test anyway.
                    
                    
                    %% Test the significance of the neural ISCM elements
                    % Returns the indicator function of the significant relationships, the test
                    % statistic matrix, and the threshold for significant results.
                    hh = tic;
                    [H,test_statistics,thresh] = multiple_testing_procedure_Cai_Liu(r_neural, null_distributions,f,FDR_alpha,Cai_Liu_thresh_vector,G,show_figures);
                    fprintf(['Stat test took ',num2str(toc(hh)),' s \n'])
                    r_neural(~logical(H)) = NaN; % set non-significant values to NaN
                    
                    % If we want the eplots to use the absolute value of
                    % the correlation
                    if plot_abs_corr
                        r_neural = abs(r_neural);
                    end

                    % Store resuts from testing
                    all_stored_c_temp(:,:,pair_counter) = r_neural; % for all recordings
                    H_stored = H_stored + H;

                catch
                end
            else

                % Load phase-ran ISCMs (and neural ISCM) 
                [r_neural,phase_ran_ISCMs,random_seed_vector,f_return,n] = load_phase_randomised_null(session,channel, save_phase_ran_ISCMs, try_load_C_flag,use_IAAFT_or_FT_phase_ran_nulls);
                if isempty(f_return) % if the phase-ran null does not exist for the session-channel pair
                    phase_ran_ISCMs = zeros(p,p,2);
                    r_neural = [];
                    nb_MCs_phase_ran(session,channel) = 0;
                    if isempty(f) % if f has never been initialised before
                        f = zeros(length(p),1);
                    end               
                else
                    f = f_return; % the vector of frequencies in the CWT
                    nb_MCs_phase_ran(session,channel) = length(phase_ran_ISCMs(1,1,:));
                end

                %% Load white noise ISCMs                    
                % Check recording length, see which standardised lengh
                % it corresponds to
                session_lengths = get_standardised_length(length_MC_vector, minimum_time_gap, round(n/Fs_neural));
                length_MC_counter = find(session_lengths == length_MC_vector); % gives the index for what part of the C_noises_all matrix to index for the given session length

                % Check the standardisation of length went correcty, if
                % not, skip this session-channel pair
                if isempty(length_MC_counter)
                    continue
                end

                % Index appropriate white-noise ISCMs
                white_noise_ISCMs = white_noise_ISCM_all{length_MC_counter}; % index correct MC length
                nb_MCs = length(white_noise_ISCMs(1,1,:));

                %% Get null distribution
                % Sum mean white-noise ISCM and phase-ra nISCM
                % distributions to get null distribution
                null_distributions = phase_ran_ISCMs + mean(white_noise_ISCMs,3) - eye(p); % deleting the identity matrix is optional, since diagonal elements are never considered in the stat test anyway.


                %% Test the significance of the neural ISCM elements
                % Returns the indicator function of the significant relationships, the test
                % statistic matrix, and the threshold for significant results.
                hh = tic;
                [H,test_statistics,thresh] = multiple_testing_procedure_Cai_Liu(r_neural, null_distributions,f,FDR_alpha,Cai_Liu_thresh_vector,G,show_figures);
                fprintf(['Stat test took ',num2str(toc(hh)),' s \n'])
                r_neural(~logical(H)) = NaN; % set non-significant values to NaN

                % If we want the eplots to use the absolute value of
                % the correlation
                if plot_abs_corr
                    r_neural = abs(r_neural);
                end

                % Store resuts from testing
                all_stored_c_temp(:,:,pair_counter) = r_neural; % for all recordings
                H_stored = H_stored + H;

            end
        end
    end

    %% Plot percentiles and quartiles
    fig = figure('Renderer', 'painters', 'Position', [150 150 800 500]);

    stored_c_temp_2 = all_stored_c_temp;
    for i = 1:p
        for j = 1:p
            temp = sort(abs(squeeze(stored_c_temp_2(i,j,:))));
            cc = sum(isnan(temp));
            temp(length(temp)-cc+1:end) = [];
            temp = [nan(cc,1); temp];
            stored_c_temp_2(i,j,:) = temp;
        end
    end

    % 1st percentile
    first_percentile_index = floor(1/100*length(stored_c_temp_2(i,j,:)));
    if first_percentile_index == 0
        first_percentile_index = 1;
    end
    subplot(2,3,1)
    xxx = abs(squeeze(stored_c_temp_2(:,:,first_percentile_index)));
    min_C_temp = min(min(xxx));
    max_C_temp = max(max(xxx));
    plot_2D_colormap_nan_ignore(xxx, f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title('(a) 1^{st} Percentile')

    % 1st quartile
    first_quartile_index = floor(25/100*length(stored_c_temp_2(i,j,:)));
    subplot(2,3,2)
    xxx = abs(squeeze(stored_c_temp_2(:,:,first_quartile_index)));
    min_C_temp = min(min(xxx));
    max_C_temp = max(max(xxx));
    plot_2D_colormap_nan_ignore(xxx, f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title('(b) 1^{st} Quartile')

    % Median
    median_index = floor(50/100*length(stored_c_temp_2(i,j,:)));
    subplot(2,3,3)
    xxx = abs(squeeze(stored_c_temp_2(:,:,median_index)));
    min_C_temp = min(min(xxx));
    max_C_temp = max(max(xxx));
    plot_2D_colormap_nan_ignore(xxx, f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title('(c) Median')

    % 3rd quartile
    third_quartile_index = floor(75/100*length(stored_c_temp_2(i,j,:)));
    subplot(2,3,4)
    xxx = abs(squeeze(stored_c_temp_2(:,:,third_quartile_index)));
    min_C_temp = min(min(xxx));
    max_C_temp = max(max(xxx));
    plot_2D_colormap_nan_ignore(xxx, f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title('(d) 3^{rd} Quartile')

    % 99th percentile
    nine_nine_index = floor(99/100*length(stored_c_temp_2(i,j,:)));
    subplot(2,3,5)
    xxx = abs(squeeze(stored_c_temp_2(:,:,nine_nine_index)));
    min_C_temp = min(min(xxx));
    max_C_temp = max(max(xxx));
    plot_2D_colormap_nan_ignore(xxx, f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title('(e) 99^{th} Percentile')

    if ~local_colorbar_bool
        % Plot global colorbar
        hp3 = get(subplot(2,3,3),'Position');
        hcb=colorbar('Position', [hp3(1)+hp3(3)+0.02 hp3(2)  0.02  0.32]);
        set(get(hcb,'label'),'string','| Correlation |');
        hcb.Ticks = linspace(0, 1, 11);
        hcb.TickLabels = num2cell(0:0.1:1); 
        caxis([0 1])
    end 

    subplot(2,3,6)
    plot_2D_colormap_nan_ignore(H_stored/pair_counter*100, f, f, 'log', second_colormap_scheme, local_colorbar_bool, plot_abs_corr, [0 100]);
    xticks([1 10 100 1000])
    yticks([1 10 100 1000])
    title(['(f) Percentage of Recordings',newline,'with significant correlations'])

    if ~local_colorbar_bool
        % Plot global colorbar
        colormap(gca,second_colormap_scheme)
        hp4 = get(subplot(2,3,6),'Position');
        hcb=colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.02  hp4(2)+hp4(3)]);
        set(get(hcb,'label'),'string','Percentage of Recordings');
        hcb.Ticks = linspace(0, 1, 11) ;
        hcb.TickLabels = num2cell(0:10:100) ;
    end

    % Plot figure global axis labels
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    han.Title.Visible='on';
    ylabel(han,'Frequency (Hz)');
    xlabel(han,'Frequency (Hz)');
    % title(han,['Inter-Scale Correlation Matrices for Channel ',num2str(channel),', with Statistically ',...
    % 'Insignificant Relationships shown in Black',newline])


    %% How many recordings have sub-10 Hz correlations with above-100 Hz?
    fig2 = figure;

    % Frequency vector indices
    [~,sub_10_index] = min(abs(f-10)); sub_10_index = sub_10_index - 1; % below 10 Hz
    [~,sup_100_index] = min(abs(f-100)); sup_100_index = sup_100_index + 1; % above 100 Hz

    for i = 1:(length(sessions_vector)*length(channels_vector))
        temp = squeeze(all_stored_c_temp(:,:,i)); % inter-scale corr. matrix from the index recording
        if sum(sum(temp==0)) == numel(temp) % if the session-channel wasn't loaded
            prop_sig_corrs(i,1) = NaN;
            continue
        end
        temp(temp==0) = NaN;

        % Count proportion of non-NaN low-high entries in the recording = proportion of stat. sig low-high relationships
        prop_sig_corrs(i,1) = sum(sum(~isnan(temp(sub_10_index:end,1:sup_100_index))))/numel(temp(sub_10_index:end,1:sup_100_index));  
    end

    x = sort(prop_sig_corrs);
    histogram(x,[0:0.01:1])
    ylabel(['Number of Recordings that achieved',newline,'the Specified Ratio of Significant',newline,'Correlations'])
    xlabel(['Ratio of Significant Low-High Correlations vs. All Tested Low-High Correlations'])
%     yticks([0:10:130])
    % grid


    %% Plot histogram of how many relationships of all kinds were significant in how many recordings
    fig3 = figure;

    for i = 1:(length(sessions_vector)*length(channels_vector))
        temp = squeeze(all_stored_c_temp(:,:,i));
        if sum(sum(temp==0)) == numel(temp) % if the session-channel wasn't loaded
             proportion_of_sig_relationships(i,1) = NaN;
            continue
        end
        temp(temp==0) = NaN;

        % Count proportion of non-NaN entries in the recording = proportion of stat. sig relationships
        proportion_of_sig_relationships(i,1) = sum(sum(~isnan(temp)))/numel(temp); 
    end

    x = sort(proportion_of_sig_relationships);
    histogram(x,[0:0.05:1])
    ylabel('Number of Recordings')
    xlabel(['Ratio of Significant Correlations vs. All Tested Correlations'])

    %% Plot positive v negative ratio
    fig3 = figure;

    % for i = 1:(length(session_vector)*length(channel_vector))
    %     temp = squeeze(all_stored_c_temp(:,:,i));
    %     if sum(sum(temp==0)) == numel(temp) % if the session-channel wasn't loaded
    %          proportion_of_sig_relationships(i,1) = NaN;
    %         continue
    %     end
    %     temp(temp==0) = NaN;
    % 
    %     % Count proportion of non-NaN entries in the recording = proportion of stat. sig relationships
    %     proportion_of_sig_relationships(i,1) = nansum(nansum(temp>0))/numel(temp); 
    % end

    figure
    x = all_stored_c_temp(~isnan(all_stored_c_temp));
    histogram(x,[-0.5:0.02:1])
    ylabel(['Number of Significant',newline,'Correlations'])
    xlabel(['Correlation Values'])

end
