function plot_time_series_all(sessions_vector,channels_vector, grid_bool)

    Inter_Scale_Stat_Sig_Sabes_Parameters;
    
    % Plot time series as 3x5 grid of subplots, where all sessions from 
    % the same channel are plotted in the same figure
    if grid_bool
        for chan = channels_vector
            
            % Plot images
            img = imread([time_domain_plots,'\Channel_',num2str(chan),'.jpg']);
            figure
            image(img);
        end
        
    % Plot time series, one figure for each channel-session pair
    elseif ~grid_bool
        
        for chan = channels_vector
            for session = sessions_vector
                % Plot images
                img = imread([time_domain_plots,'\Channel_',num2str(chan),'\',processed_sessions{session},'.jpg']);
                figure
                image(img);
            end
        end
        
    end
end