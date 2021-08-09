
function collate_IAAFT_nulls(session,channel,phase_ran_results_folder,try_load_flag)
    
    % Collates the IAAFT nulls across all of the different parallel runs, for the
    % specified channel-session pair.

    % Finds the maximum amount of parallel phase-ran runs per channel, assumes the
    % run number is at the end of the file name and preceded by an
    % underscore '_'.
    [phase_ran_results_folder,'\Session_',num2str(session)]
    d = dir([phase_ran_results_folder,'\Session_',num2str(session)]); d(1:2) = [];
    run = zeros(length(d),1);
    for i = 1:numel(d)
        s = d(i).name;
        counter = length(s)-5;
        run(i,1) = NaN;
        while isnan(run(i,1))
            counter = counter - 1;
            run(i,1) = str2double(s(end-(4+counter):end-4));
        end
    end

    % Iterate across parallel runs
    x = [];
    random_seed_vector = [];
    for run = 1:max(run)
%         if try_load_flag

            try 
                phase_file = [phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'_run_',num2str(run),'.mat'];
                load(phase_file)

                % The parallel IAAFT runs (different stored files) are
                % combined together, here we figure out the indices
                if isempty(x)
                    counter = 0;
                else
                    counter = length(x(1,:));
                end
                if exist('r_phase')
                    counter_2 = length(compressed_r_phase(1,:));
                else
                    counter_2 = 0;
                end
                x = [x compressed_r_phase];
%                 x(:,counter+1:counter+counter_2) = compressed_r_phase;
                delete(phase_file)
                clear counter_2 r_phase
    %                 random_seed_vector = [random_seed_vector random_seed];
            catch
                fprintf(['Failed to load Session ' ,num2str(session),' - Channel ',num2str(channel),' - Run ',num2str(run),'\n'])
            end

%         else % if we have the try load flag = false
%             phase_file = [phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'_run_',num2str(run),'.mat'];
%             load(phase_file)
% 
%             % Decompress data (it is stored as a 2D matrix, where each column
%             % represents 1 phase-ran/white-noise ISCM, and redundant 
%             % info due to the symmetrical  nature of the ISCM is 
%             % thrown away.
%             for MC_iteration = 1:length(compressed_r_phase(1,:))
%                 r_phase(:,:,MC_iteration) = reshape_1D_vector_to_2D_symmetric_matrix(compressed_r_phase(:,MC_iteration));
%             end
% 
%             % The parallel IAAFT runs (different stored files) are
%             % combined together, here we figure out the indices
%             if isempty(r_phase_collated)
%                 counter = 0;
%             else
%                 counter = length(r_phase_collated(1,1,:));
%             end
%             if exist('r_phase')
%                 counter_2 = length(r_phase(1,1,:));
%             else
%                 counter_2 = 0;
%             end
% 
%             % Collate IAAFTs across different parallel runs
%             r_phase_collated(:,:,counter+1:counter+counter_2) = r_phase;
%             clear counter_2 r_phase
%             delete(phase_file)
%             random_seed_vector = [random_seed_vector random_seed];
%         end
    end
%     phase_file = [phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'_run_',num2str(run),'.mat'];
    compressed_r_phase = x;
    save([phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat'],'n','S_neural_extract','clipping_counter','x_sample_pre_processed','x_sample_orig','times','compressed_r_phase','r_orig','f','channel','session','cores','PBS_array_index_orig','MCs_per_channel','random_seed_vector');
        
end
