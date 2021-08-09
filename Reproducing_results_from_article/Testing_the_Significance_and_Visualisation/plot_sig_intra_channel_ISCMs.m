% Function for plotting the intra-channel inter-frequency correlation
% matrices, after statistical testing. This plots a 3x5 subplot array, for
% the 15 sessions that were analysed in this work. This script iterates
% through the channel vector, i.e. if channel_vector has 20 elements, then
% there will be 20 figures. Black parts in the plots represent
% non-statistically significant correlations. Each channel and session pair
% is tested individually, and the FDR is not controlled across all plots
% simulataneously.
%
% Note: Edit the subplot configuration if not all 15 sessions are analysed.

function plot_sig_intra_channel_ISCMs(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters)


    Inter_Scale_Stat_Sig_Sabes_Parameters; % initialise parameters used in this work

    % Parameters for this script, give them the shorthand in a very inefficient way because I'm lazy
    show_figures = meta_parameters.show_figures; 
    save_stat_sig_inter_scale_corr_plots = meta_parameters.save_stat_sig_inter_scale_corr_plots; % whether to save the plots of the statistically significant ISCMs to folder "Figures_stat_sig_ISCMs", specified in "Inter_Scale_Stat_Sig_Sabes_Parameters".
    try_load_C_flag = meta_parameters.try_load_C_flag; % if it fails to load each specified neural ISCM, it continues if this flag is true. If false, it stops when it fails. False is better for testing, true is better for skipping over missing results.
    try_plot_flag = meta_parameters.try_plot_flag; % if it fails to plot the ISCM, then if 'true', it continues. If 'false', it doesn't.
    local_colorbar_bool = meta_parameters.local_colorbar_bool; % whether to plot the colorbar individually for each ISCM
    plot_abs_corr = meta_parameters.plot_abs_corr; % whether to plot the correlations (false) or the absolute values of the correlations (true). The absolute values are less informative, but make the color scheme more intuitive.
    colormap_scheme = meta_parameters.colormap_scheme; % what MATLAB colormap scheme to use
    
    % Choose which phase-ran ISCMs are used
    if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'FT') % use FT phase-ran nulls
        save_phase_ran_ISCMs = save_FT_phase_ran_ISCMs;
    elseif strcmp(use_IAAFT_or_FT_phase_ran_nulls,'IAAFT') % use IAAFT phase-ran nulls
        save_phase_ran_ISCMs = save_IAAFT_phase_ran_ISCMs;
    else
        error('use_IAAFT_or_FT_phase_ran_nulls should be a string equal to FT or IAAFT')
    end
    
    % Initialise some values
    all_test_stat_stored = zeros(p,p,length(sessions_vector)*length(channels_vector));
    all_stored_c_temp = zeros(p,p,length(sessions_vector)*length(channels_vector));

    % Check
    if max(sessions_vector) > numel(processed_sessions)
        error(['There are elements in the session_vector that do not correspond to the sessions listed in processed_sessions'])
    end

    % If the white noise ISCMs haven't been collated into the same file yet, do so
    if exist([save_white_noise_ISCMs,'\All_White_Noise_MCs_Length_',num2str(length_MC_vector(3)),'.mat']) ~= 2 % check if the MCs haven't been placed into the same file
        format_WN_ISCM_data;
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
    clear C_noises WI_ISCM

    
    %% Iterate through session-channel pairs
    % Plots the (significant values of the) inter-scale correlation matrices 
    % of each session
    chan_counter = 0;
    counter = 0; 
    for channel = channels_vector
        chan_counter = chan_counter + 1;
        clear stored_c_temp
        fig = figure('Renderer', 'painters', 'Position', [50 50 1150 600]);
        session_counter = 0;
        max_C_temp = 0; min_C_temp = 1; % two values used to calibrate the global colorar bar per channel
        for session = sessions_vector % iterate through the listed sessions in each channel
            session_counter = session_counter + 1;
            counter = counter + 1;

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
                    session_length = get_standardised_length(length_MC_vector, minimum_time_gap, round(n/Fs_neural));
                    length_MC_counter = find(session_length == length_MC_vector); % gives the index for what part of the C_noises_all matrix to index for the given session length
             
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

                    %% Store resuts from testing
                    test_statistics(~logical(H)) = 0;
                    all_test_stat_stored(:,:,counter) = test_statistics;

                    stored_c_temp(:,:,session_counter) = r_neural; % for each channel in that session, for plotting
                    all_stored_c_temp(:,:,counter) = r_neural; % for all recordings

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
                session_length = get_standardised_length(length_MC_vector, minimum_time_gap, round(n/Fs_neural));
                length_MC_counter = find(session_length == length_MC_vector); % gives the index for what part of the C_noises_all matrix to index for the given session length

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

                %% Store resuts from testing
                test_statistics(~logical(H)) = 0;
                all_test_stat_stored(:,:,counter) = test_statistics;

                stored_c_temp(:,:,session_counter) = r_neural; % for each channel in that session, for plotting
                all_stored_c_temp(:,:,counter) = r_neural; % for all recordings

            end
            min_C_temp = min([min(min(r_neural)) min_C_temp]);
            max_C_temp = max([max(max(r_neural)) max_C_temp]);
        end

        for session = sessions_vector % iterate through the listed sessions in each channel, for plotting

            if try_plot_flag
                try
                    % Setting colarmap limits
                    if local_colorbar_bool % if we're treating each recording independently
                        min_C_temp = min(min(squeeze(stored_c_temp(:,:,session))));
                        max_C_temp = max(max(squeeze(stored_c_temp(:,:,session))));
                    end

                    % Plotting stat. sig elements of the inter-frequency matrices
                    hh = tic;
                    subplot(3,5,session)
                    plot_2D_colormap_nan_ignore(squeeze(stored_c_temp(:,:,session)), f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
                    xticks([1 10 100 1000])
                    yticks([1 10 100 1000])
                    title(['Session ',num2str(session),' - Chan. ',num2str(channel)])
    %                 title(['Session ',num2str(session),' - Chan. ',num2str(channel),'; Nulls: ',num2str(nb_MCs_phase_ran(session,channel))])
%                     title(['Sess. ',num2str(session),' - Chan. ',num2str(channel),'; Nulls: ',num2str(nb_MCs_phase_ran(session,channel))])   
                    drawnow;
                    fprintf(['Plotting took ',num2str(toc(hh)),' s \n'])

                catch
                end
            else
                % Setting colormap limits
                if local_colorbar_bool % if we're treating each recording independently
                    min_C_temp = min(min(squeeze(stored_c_temp(:,:,session))));
                    max_C_temp = max(max(squeeze(stored_c_temp(:,:,session))));
                end

                % Plotting stat. sig elements of the inter-frequency matrices
                hh = tic;
                subplot(3,5,session)
                plot_2D_colormap_nan_ignore(squeeze(stored_c_temp(:,:,session)), f, f, 'log', colormap_scheme, local_colorbar_bool, plot_abs_corr, [min_C_temp max_C_temp]);
                xticks([1 10 100 1000])
                yticks([1 10 100 1000])
                title(['Session ',num2str(session),' - Chan. ',num2str(channel)])
%                 title(['Sess. ',num2str(session),' - Chan. ',num2str(channel),'; Nulls: ',num2str(nb_MCs_phase_ran(session,channel))])   
                drawnow;
                fprintf(['Plotting took ',num2str(toc(hh)),' s \n'])
            end
        end

        % Plot global colorbar
        if ~local_colorbar_bool 
            hp4 = get(subplot(3,5,15),'Position');
            hcb=colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(2)+hp4(3)*5.7]);
            % title(hcb,'Correlation')
            if ~plot_abs_corr
                set(get(hcb,'label'),'string','Correlation');
            else
                set(get(hcb,'label'),'string','| Correlation |');
            end
            caxis([min_C_temp max_C_temp])
        end

        % Plot figure axis labels
        han=axes(fig,'visible','off'); 
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        han.Title.Visible='on';
        ylabel(han,'Frequency (Hz)');
        xlabel(han,'Frequency (Hz)');
        title(han,['Inter-Scale Power Correlation Matrices for Channel ',num2str(channel),', with Statistically ',...
        'Insignificant Relationships shown in Black',newline])
        drawnow;

        % Can save figures to a folder "save_stat_sig_inter_scale_plots"
        % if desired
        if save_stat_sig_inter_scale_corr_plots
            [~,~] = export_fig([Figures_stat_sig_ISCMs,'\Channel_',num2str(channel)],'-jpg','-m12');
            close all
        end
        
%         if rem(channel,20) == 0
%             close all
%         end
    end

end