% Function for plotting and saving the time domain signals
% of the processed data. This can help in identifying line noise,
% artifacts and pathological activity.

function save_time_series_of_processed_data(sessions_vector,channels_vector, grid_bool)

    % Initialise parameters and directories
    Inter_Scale_Stat_Sig_Sabes_Parameters;

    % Plot time series, one figure for each channel-session pair
    if ~grid_bool
        % Iterate through selected sessions
        directory_marker = 0;
        for session = sessions_vector

            % Load data from session
            load([processed_data_folder,'\',processed_sessions{session},'.mat'])

            % Iterate through selected channels
            for chan = channels_vector

                % Create figure directories
                if directory_marker == 0
                    mkdir([time_domain_plots,'\Channel_',num2str(chan)])
                end

                % Check if the figure exists already
                if isfile([time_domain_plots,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg']) %Session_',num2str(session)
                    fprintf([num2str( [chan session]),' done already \n'])
                    continue
                end

                % Print channel-session pair, and do some formatting
                fprintf(['Channel: ',num2str(chan),'; Session: ',num2str(session),'\n'])
                temp_name = processed_sessions{session};
                temp_name = strrep(temp_name,'_','-');

                % Plot and save time series
                figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
                plot(neural_time,SI_CSD_Broadband(:,chan))
                xlabel('Time (s)'); ylabel('Amplitude')
                title([temp_name,'   (Session ',num2str(session),')   -   Channel ',num2str(chan)])
                [~,~] = export_fig([time_domain_plots,'\Channel_',num2str(chan),'\',processed_sessions{session}],'-jpg','-m6');

        %         return
                close all

            end
            directory_marker = 1;
        end
        
    % Plot time series as 3x5 grid of subplots, where all sessions from 
    % the same channel are plotted in the same figure
    elseif grid_bool
        session_counter = 0;
        for session = session_vector
            session_counter = session_counter + 1;

            % Load data from session
            load([processed_data_folder,'\',processed_sessions{session},'.mat'])

            % Iterate through selected channels
            for chan = channel_vector

                % Check if the figure exists already
                if isfile([time_domain_plots,'\Channel_',num2str(chan),'.jpg'])
                    fprintf([num2str([chan session]),' done already \n'])
                    continue
                end

                fprintf(['Doing channel ',num2str(chan),'; Session ',num2str(session),' \n'])
      
                if session_counter == 1
                    figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
                end

                figure(chan); % index correct figure
                subplot(3,5,session_counter)
                plot(neural_time,SI_CSD_Broadband(:,chan)); % plot time domain

                % Save figure
                if session_counter == numel(sessions_vector)
%                     temp_name = processed_sessions{session};
%                     temp_name = strrep(temp_name,'_','-');
                    title(['Session ',num2str(session)])
                    drawnow;
                    [~,~] = export_fig([time_domain_plots,'\Channel_',num2str(chan)],'-jpg','-m6');
                end
            end
        end
        close all
    end
       
end