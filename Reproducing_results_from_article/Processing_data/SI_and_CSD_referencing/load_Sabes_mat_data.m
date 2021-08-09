function [raw_Broadband, neural_time] = load_Sabes_mat_data(raw_neural_data_path,intro_cutoff)
    
    % Loads the neural data with time stamps, and cuts off the beginning as specified by intro_cutoff.
    load(raw_neural_data_path)
    raw_Broadband = double(stored_data.neural_recordings(intro_cutoff:end,:));
    neural_time = stored_data.time(intro_cutoff:end);
end