function [WN_ISCM_individual,f] = generate_WN_ISCMs(seed_1,seed_2)

    % Generate inter-scale correlation matrices from white noise input, 
    % used in the statistical significance test. In this work, this function is 
    % called by main_WN_ISCM_MC_simulation, a BASH script that calls this
    % MATLAB function and another C function (that produces the PCG
    % numbers).
    %
    % NOTE: DESIGNED FOR LINUX.
    %
    % Inputs:
    % - seed_1: RNG state constant, equal to 1 in this work for all MC iterations
    % - seed_2: RNG sequence seed   
    %    
    % Implicit input:
    % - random_numbers_file: a .txt file containing the PCG random numbers. It's name is [seed_2 value].txt
    %
    % Outputs: 
    % WN_ISCM_individual: inter-scale correlation matrix, produced from random white noise.
    % f: frequency vector.
    
    % Intitialise parameters used in this work
    Inter_Scale_Stat_Sig_Sabes_Parameters;
    save_results = 1;
    
    
    seed_2
    
    
    
    % File random numbers are saved to
    random_numbers_file = [random_numbers_folder,num2str(seed_2),'.txt'];
    random_numbers_file

    %% Get HPC save path
    % Observe seed_2 = n+v, from there derive if the MC is
    % of length 350, 400 or 500 s. Set the save path. This needs to match up with
    % length_MC_vector in "Inter_Scale_Stat_Sig_Sabes_Parameters.m".
    save_path = save_white_noise_ISCMs;
    if seed_2 > 8544500 && seed_2 < 9765519 % if index between 350 and 400, store as 350 s
        save_path = [save_path,'/Length_350'];
    elseif seed_2 > 9765500 && seed_2 < 9968625 % if index between 400 and ~410, store as 400 s
        save_path = [save_path,'/Length_400'];
    elseif seed_2 > 12206000 && seed_2 < 12300031 % if index between 500 and ~510, store as 500 s
        save_path = [save_path,'/Length_500'];
    else
      error('signal_length does not correspond to the approved ranges for 350, 400 or 500 s processes, sampled at Fs')
    end
    mkdir(save_path)
    
    % 8544922 samples -> 350 s @ at Fs of 24414... Hz
    % 9765625 samples -> 400 s.
    % 12207031 samples -> 500 s.


    %% If this WN_ISCM has already been calculated, delete the .txt file and return
    if isfile([save_path,'/',num2str(seed_2),'.mat'])
        delete(random_numbers_file)
        fprintf('Done already \n')
        return
    end
    
    %% Extract the PCG random numbers
    fileID = fopen(random_numbers_file,'r');
    formatSpec = '%d';
    x_noise = fscanf(fileID,formatSpec); % store PCG numbers as x_noises
    fclose(fileID);
    random_numbers_file % print PCG numbers file name
    delete(random_numbers_file) % delete the .txt file with PCG numbers
    x_noise(1:5)' % print some of the random numbers for fun
    
    %% Clip x_noise
    x_noise = rescale(x_noise,0,1);
    n_noise = find_matching_end(x_noise,clipping_epsilon);
    x_noise = x_noise(1:n_noise);
    
%     %% Pre-whiten x_noise
%     x_noise = whitening(x_noise, Fs_neural);
%     
    %% Calculate Inter-Scale Correlation Matrix (ISCM)
    [WN_ISCM_individual,f,computation_times] = calculate_inter_scale_correlation_matrix(x_noise, Fs_neural, CWT_freq_limits,thresh_for_corr_inclusion);
    
    %% Compress ISCM, take advantage of symmetrical nature
    [WN_ISCM_individual] = reshape_2D_symmetric_matrix_to_1D(WN_ISCM_individual);
    if ~isequal(length(WN_ISCM_individual),(p^2-p)/2)
        WN_ISCM_individual = zeros((p^2-p)/2,1);
    end
        
    %% Save results
    if save_results == 1
      hh= tic;
        save([save_path,'/',num2str(seed_2),'.mat'],'computation_times','WN_ISCM_individual','f','seed_1','seed_2','n_noise');
        saving_time = toc(hh);
        fprintf(['Saving time was ',num2str(saving_time),' s \n'])
        [save_path,'/',num2str(seed_2),'.mat'] % print save file name, for fun
    end

end
