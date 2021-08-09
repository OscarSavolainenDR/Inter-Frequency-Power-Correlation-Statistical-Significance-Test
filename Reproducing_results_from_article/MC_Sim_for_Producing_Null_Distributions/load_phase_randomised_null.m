function [r_orig,r_phase_collated,random_seed_vector,f,n] = load_phase_randomised_null(session,channel, phase_ran_results_folder, try_load_flag, use_IAAFT_or_FT_phase_ran_nulls)

    %% If IAAFT phase-ran, the parallel runs need to be collated
    if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'IAAFT')  
        % If session-channel pair files have not been collated, do so 
        if ~isfile([phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat']) 
            collate_IAAFT_nulls(session,channel,phase_ran_results_folder,try_load_flag);
        end  
    end

    %% Load data
    if try_load_flag
        try 
            phase_file = [phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat'];
            load(phase_file)

            % Decompress data (it is stored as a 2D matrix, where each column
            % represents 1 phase-ran/white-noise ISCM, and redundant 
            % info due to the symmetrical  nature of the ISCM is 
            % thrown away.
            for MC_iteration = 1:length(compressed_r_phase(1,:))
                r_phase_collated(:,:,MC_iteration) = reshape_1D_vector_to_2D_symmetric_matrix(compressed_r_phase(:,MC_iteration));
            end
            if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'FT') % just a naming convention
                random_seed_vector = random_seed;
            end
        catch
            fprintf(['Failed to load Session ' ,num2str(session),' - Channel ',num2str(channel),'\n'])
        end
    else
        phase_file = [phase_ran_results_folder,'\Session_',num2str(session),'\Channel_',num2str(channel),'.mat'];
        load(phase_file)

        % Decompress data (it is stored as a 2D matrix, where each column
        % represents 1 phase-ran/white-noise ISCM, and redundant 
        % info due to the symmetrical  nature of the ISCM is 
        % thrown away. Here we decmpress by reconstructing the 3D mtrix of
        % white-noise ISCMs, where each slice in the 3rd dimension is one
        % ISCM.
        for MC_iteration = 1:length(compressed_r_phase(1,:))
            r_phase_collated(:,:,MC_iteration) = reshape_1D_vector_to_2D_symmetric_matrix(compressed_r_phase(:,MC_iteration));
        end
        if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'FT') % just a naming convention
            random_seed_vector = random_seed;
        end
    end

    %% Throw away empty null ISCMs (where the procedure failed, and it's just a matrix of zeros)
    if ~isempty(r_phase_collated)
        size_r = size(r_phase_collated);
        r_rows = reshape(r_phase_collated,[],size_r(3))';

        r_phase_collated = unique(r_rows,'stable','rows');
        r_phase_collated = reshape(r_phase_collated',size_r(1),size_r(2),[]);
    else % if nothing loaded, add dumby variables
        r_orig = [];
        f = [];
        n = 503*24414;
    end
end
