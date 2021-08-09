% The FT phase-randomised ISCMs are stored in batches of 10 (or in one case
% 16) channels, due to memory constraints on Zendodo. This function
% collates all of the data into the same folder. The FT phase-ransomised
% ISCM .zip files should be loaded into the same folder, the path of which
% is specified in the 'save_FT_phase_ran_ISCMs' variable in the 
% "Inter_Scale_Stat_Sig_Sabes_Parameters" script.
% This function will collate all of the unzipped folders into the same
% directory, and delete the original seperate folders. The .zip files can
% then be deleted seperately by the user.

function merge_FT_phase_ran_folders(sessions_vector,channels_vector)

    Inter_Scale_Stat_Sig_Sabes_Parameters;

    flag_1 = false;
    for session = sessions_vector
        for channel = channels_vector
            if ~isfile([save_FT_phase_ran_ISCMs,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat'])
                flag_1 = true;
            end
        end
    end
    
    counter = 0; flag_4 = 0; i = 1;
    while flag_4 == 0 
        if counter + 5 > max(channels_vector)
            channels_sub{i} = counter+1:max(channels_vector);
            flag_4 = 1;
        else
            channels_sub{i} = counter+1:counter + 5;
        end
        counter = counter + 5;
        i = i + 1;
    end


    if flag_1
        fprintf('Merging FT phase-ransomised null folders \n')
        for iteration = 1:i-1
%             iteration

            for session = sessions_vector
%                 session
                try
                    for channel = channels_sub{iteration}
    %                     channel
                        data_folder = [save_FT_phase_ran_ISCMs,'\Channels_',num2str(min(channels_sub{iteration})),'_to_',num2str(max(channels_sub{iteration})),'_FT_distributions'];
                        data_sub_folder = [data_folder,'\Session_',num2str(session)];
                        data_file = [data_sub_folder,'\Channel_',num2str(channel),'.mat'];

                        save_path = [save_FT_phase_ran_ISCMs,'\Session_',num2str(session)];
                        if exist(save_path) ~= 7
                            mkdir(save_path)
                        end
%                         data_file
%                         save_path
                        copyfile(data_file,save_path)
                        delete(data_file)
                    end
                catch
                end
                
            end
            [~,~] = cmd_rmdir(data_folder);
        end

    end
    
    for session = sessions_vector
        for channel = channels_vector
            save_path = [save_FT_phase_ran_ISCMs,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat'];
            if ~isfile(save_path) 
                fprintf(['Session ',num2str(session),' - Channel ',num2str(channel),' missing \n'])
%                 save_path
%                 error(['Session ',num2str(session),' - Channel ',num2str(channel),' missing'])
            end
        end
    end
end