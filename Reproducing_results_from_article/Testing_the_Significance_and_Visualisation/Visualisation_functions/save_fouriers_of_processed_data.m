% Script for plotting and saving the Fourier domains
% of the processed data. This can help in identifying line noise,
% artifacts and pathological activity.

function save_fouriers_of_processed_data(session_vector,channel_vector,broadband_range,mid_range, grid_bool)

    Inter_Scale_Stat_Sig_Sabes_Parameters;

    % Plot Fouriers individually, one figure for each channel-session pair
    % for high-res
    if ~grid_bool
        % Iterate through selected sessions
        directory_marker = 0;
        for session = session_vector

            % Load data from session
            load([processed_data_folder,'\',processed_sessions{session},'.mat'])

            % Iterate through selected channels
            for chan = channel_vector

                % Create figure directories
                if directory_marker == 0
                    mkdir([FT_plots_CSD_Broadband_folder_midrange,'\Channel_',num2str(chan)])
                    mkdir([FT_plots_CSD_Broadband_folder_broadband,'\Channel_',num2str(chan)])
                end

                % Check if the figure exists already
                if isfile([FT_plots_CSD_Broadband_folder_midrange,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg'])
                    fprintf([num2str( [chan session]),' done already \n'])
                    continue
                elseif isfile([FT_plots_CSD_Broadband_folder_broadband,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg'])
                    fprintf([num2str( [chan session]),' done already \n'])
                    continue
                end

                % Print channel-session pair, and do some formatting
                fprintf(['Channel: ',num2str(chan),'; Session: ',num2str(session),'\n'])
                temp_name = processed_sessions{session};
                temp_name = strrep(temp_name,'_','-');

                % Plot Fourier, broadband
                figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
                line_FFT_xlim(normalize(SI_CSD_Broadband(:,chan)),Fs_neural, broadband_range); % plot Fourier
                title([temp_name,'   (Session ',num2str(session),')   -   Channel ',num2str(chan)])
                [~,~] = export_fig([FT_plots_CSD_Broadband_folder_broadband,'\Channel_',num2str(chan),'\',processed_sessions{session}],'-jpg','-m6'); % save figure

                % Plot Fourier, medium frequencies
                figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
                line_FFT_xlim(normalize(SI_CSD_Broadband(:,chan)),Fs_neural, mid_range); % plot Fourier
                title([temp_name,'   (Session ',num2str(session),')   -   Channel ',num2str(chan)])
                [~,~] = export_fig([FT_plots_CSD_Broadband_folder_midrange,'\Channel_',num2str(chan),'\',processed_sessions{session}],'-jpg','-m6'); % save figure

                close all

            end
            directory_marker = 1;
        end
    
    % Plot as 5x3 grid, easier to compare to ISCMs.
    else 
        
        for freq_range_counter = 1:2 % broadband or mid-range frequencies
            
            if freq_range_counter == 1 % If midrange frequencies
                freq_range = mid_range;
                fourier_directory = FT_plots_CSD_Broadband_folder_midrange;
            elseif freq_range_counter == 2 % If broadband frequencies
                freq_range = broadband_range;
                fourier_directory = FT_plots_CSD_Broadband_folder_broadband;
            end       

            session_counter = 0;
            for session = session_vector
                session_counter = session_counter + 1;
               
                % Load data from session
                load([processed_data_folder,'\',processed_sessions{session},'.mat'])

                % Iterate through selected channels
                for chan = channel_vector
                    
                    % Check if the figure exists already
                    if isfile([fourier_directory ,'\Channel_',num2str(chan),'.jpg'])
                        fprintf([num2str( [chan session]),' done already \n'])
                        continue
                    end
                    
                     if freq_range_counter == 1 % If midrange frequencies
                        fprintf(['Doing channel ',num2str(chan),'; Session ',num2str(session),'; Midrange frequencies \n'])
                    elseif freq_range_counter == 2 % If broadband frequencies
                        fprintf(['Doing channel ',num2str(chan),'; Session ',num2str(session),'; Midrange frequencies \n'])
                    end       
                
                    if session_counter == 1
                        figure('Renderer', 'painters', 'Position', [50 50 1200 600]);
                    end

                    figure(chan); % index correct figure
                    subplot(3,5,session_counter)
                    line_FFT_xlim(normalize(SI_CSD_Broadband(:,chan)),Fs_neural, freq_range); % plot Fourier
                    
                    % Save figure
                    if session_counter == numel(session_vector)
%                         temp_name = processed_sessions{session};
%                         temp_name = strrep(temp_name,'_','-');
                        title(['Session ',num2str(session)])
                        drawnow;
                        [~,~] = export_fig([fourier_directory,'\Channel_',num2str(chan)],'-jpg','-m6');
                    end
                end
            end
            close all
        end
    end 
end