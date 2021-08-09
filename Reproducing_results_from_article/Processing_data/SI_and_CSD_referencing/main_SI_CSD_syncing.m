% Function for processing the raw neural data (.mat files) from Sabes lab. 
% Spectrally Interpolates the line noise and Current Source Density 
% references the data. It is also optional to sync the Neural and 
% Behavioral data.

function main_SI_CSD_syncing(sessions_vector,sync_neural_and_behavioral)

    Inter_Scale_Stat_Sig_Sabes_Parameters;

    show_figures = false;

    mkdir(processed_data_folder) % where to save the processed data
    d_neural = dir(raw_neural_data_folder); d_neural(1:2) = []; % neural sessions  

    for session_index = sessions_vector
        hh = tic;
        neural_path = [d_neural(session_index).folder,'\',processed_sessions{session_index},'.mat'];
        behav_path = [behavioral_data_folder,'\',processed_sessions{session_index},'.mat'];
        
        if ~isfile([processed_data_folder,'\',processed_sessions{session_index},'.mat'])
            [~,~,~,~,~,~] = pre_processing_and_syncing([processed_sessions{session_index},'.mat'], processed_data_folder, neural_path, session_index, behav_path, show_figures, sync_neural_and_behavioral);
            fprintf(['Took ',num2str(toc(hh)),' s \n'])
        else
            fprintf([processed_sessions{session_index},'.mat exists already, no processing required \n'])
        end
    end

end