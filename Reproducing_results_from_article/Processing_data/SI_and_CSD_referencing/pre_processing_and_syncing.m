function [SI_CSD_Broadband, neural_time, behav_time, finger_pos, cursor_pos, target_pos] = pre_processing_and_syncing(filename, save_data_folder, raw_neural_data_path, recording_index, behavioral_data_path, show_figures, sync_neural_and_behavioral)
    % This function goes through the entire pre- and post-processing to produce
    % the inter-scale correlation matrices from the raw neural data.
    % However, due to the size of the Sabes lab dataset, it
    % requires significant memory and so is impractical for most computers. This
    % function may serve more as a guide, and may require deconstruction
    % depending on the user's computational constraints.
    % If sync_neural_and_behavioral == true, then this function also syncs
    % the neura and behavioral data. If one is not interested in the
    % behavioral data, sync_neural_and_behavioral can be set to false to
    % save memory and time.
    

    % Initialise
    Inter_Scale_Stat_Sig_Sabes_Parameters;
    behav_time = [];
    finger_pos = [];
    cursor_pos = [];
    target_pos = [];
    
  
    %% Pre-processing
    %fprintf('Raw data has not yet been processed. \nThe neural data will be Spectrally Interpolated, Current Source \nDensity (CSD) Referenced, and then synced with the Behavioral data. \n')

    % Load raw data, have as function
    [raw_Broadband, neural_time] = load_Sabes_mat_data(raw_neural_data_path,intro_cutoff);

    if sync_neural_and_behavioral
        % Load Behavioral data
        load(behavioral_data_path);
    end

    % Spectral Interpolation (SI) and Current Source Density (CSD) Referencing
    [SI_CSD_Broadband] = csd_lfp_removal_spectral_interpolation(recording_index,raw_Broadband,show_figures,figures_cutoff);
    neural_time(length(SI_CSD_Broadband)+1:end) = []; % I know how the data is processed, only parts from the end of the neural data have been removed
    clear raw_Broadband

    if sync_neural_and_behavioral
        % Synchronises the Neural and Behavioral data
        % NOTE: if one wants to sync the Neural and Behavioral data, this function can be run.
        [neural_time, SI_CSD_Broadband, behav_time, finger_pos, cursor_pos, target_pos] = Syncing_data(neural_time, SI_CSD_Broadband, t, finger_pos, cursor_pos, target_pos, Fs_neural, Fs_behav);
    end

    % Save SI CSD-referenced data
    if sync_neural_and_behavioral
        save([save_data_folder,'\',filename],'SI_CSD_Broadband','neural_time','behav_time', 'finger_pos', 'cursor_pos', 'target_pos','-v7.3')
    else
        save([save_data_folder,'\',filename],'SI_CSD_Broadband','neural_time','-v7.3');
    end

end
