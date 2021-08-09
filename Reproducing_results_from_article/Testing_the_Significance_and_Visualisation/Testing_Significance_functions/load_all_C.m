function [all_C,session_lengths,f] = load_all_C(path_C, session_vector, processed_sessions, channel_vector, try_flag)
 
    session_counter = 0;
    for session = session_vector
        session_counter = session_counter + 1;
        chan_counter = 0;
        
        % Extract inter-scale correlation matrices
        if try_flag
            try
                load([path_C,'\',processed_sessions{session},'_neural_ISCM.mat'])
                [neural_results] = decompress_neural_ISCMs(stored_neural_ISCMs,neural_parameters,channel_vector); % decompress the neural ISCMs
                f = neural_results{1,1}.frequency_vector;
                for channel = channel_vector
                    chan_counter = chan_counter + 1;
                    all_C{session_counter,chan_counter} = neural_results{channel,1}.inter_scale_correlation_matrix;
                end
            catch
                 fprintf(['Session ',num2str(session),' (',processed_sessions{session},') failed \n'])
                 for channel = channel_vector
                    chan_counter = chan_counter + 1;
                    all_C{session_counter,chan_counter} = [];
                 end
            end

        else
            load([path_C,'\',processed_sessions{session},'_neural_ISCM.mat']) 
            [neural_results] = decompress_neural_ISCMs(stored_neural_ISCMs,neural_parameters,channel_vector); % decompress the neural ISCMs
            f = neural_results{1,1}.frequency_vector;
            for channel = channel_vector
                chan_counter = chan_counter + 1;
                all_C{session_counter,chan_counter} = neural_results{channel,1}.inter_scale_correlation_matrix;
            end    
        end
        
        % Get standardised recording lengths
        if try_flag
            try
                session_lengths(session_counter,1) = neural_results{1,1}.neural_signal_length_n; % using channel 1, but any channel will do
            catch
                session_lengths(session_counter,1) = 0;
            end
        else
            session_lengths(session_counter,1) = neural_results{1,1}.neural_signal_length_n; % using channel 1, but any channel will do
        end
    end


end
