function plot_Fouriers_all(sessions_vector,channels_vector)

    Inter_Scale_Stat_Sig_Sabes_Parameters;
  
    for chan = channels_vector
        for session = sessions_vector
            for freq_range_counter = 1:2 % broadband or mid-range frequencies

                % Plot images
                if freq_range_counter == 1 % Mid frequencies
                    img = imread([FT_plots_CSD_Broadband_folder_midrange,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg']);
                    figure
                    image(img);
                elseif freq_range_counter == 2 % Broadband frequencies
                    img = imread([FT_plots_CSD_Broadband_folder_broadband,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg']);
                    figure
                    image(img);
                end
                drawnow
            end
        end
    end

end
