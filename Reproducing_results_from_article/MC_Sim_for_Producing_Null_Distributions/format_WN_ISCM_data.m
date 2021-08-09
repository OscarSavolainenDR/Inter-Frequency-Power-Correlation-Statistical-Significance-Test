% format_WN_ISCM_data;

% Script for storing the white noise ISCMs into 
% the same files. The MC inter-scale correlation matrices for each length
% (e.g. 350, 400, 500 s) go into the same '.mat' file. 
% This script does not require paramaters, it uses the 
% 'save_white_noise_ISCMs' folder specified in 
% 'Inter_Scale_Stat_Sig_Sabes_Parameters.m'.

% Initialise
% Inter_Scale_Stat_Sig_Sabes_Parameters;


%% Iterate through MC lengths, i.e. for lengths of 350, 400 and 500 s
for length_MC = length_MC_vector

    clear PCG_parameters WN_ISCM
    
    %% Load WN_ISCMs, collate
    d = dir([save_white_noise_ISCMs,'\Length_',num2str(length_MC)]); d(1:2) = [];
    WN_ISCM = single(zeros(((p*p-p)/2),numel(d)));
    for i = 1:numel(d)
        load([save_white_noise_ISCMs,'\Length_',num2str(length_MC),'\',d(i).name])
        WN_ISCM(:,i) = WN_ISCM_individual;%reshape_1D_vector_to_2D_symmetric_matrix(WN_ISCM_individual);
        
        % Store phase-ran seed, in case someone wants to re-create the
        % results.
%         phase_rand_seed(i,1) = random_seed; % 
    end
    
    white_noise_ISCM_stats.mean = mean(WN_ISCM,3); % calculate prior
    white_noise_ISCM_stats.std = std(WN_ISCM,0,3); % calculate prior
   
    save([save_white_noise_ISCMs,'\All_White_Noise_MCs_Length_',num2str(length_MC),'.mat'],'WN_ISCM','white_noise_ISCM_stats')
    
end

clear channel session length_MC d phase_random_seed i WN_ISCM_individual WN_ISCM