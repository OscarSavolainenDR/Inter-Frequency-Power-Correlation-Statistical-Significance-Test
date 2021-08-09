function [WN_ISCM_individual,f] = arbitrary_generate_WN_ISCMs(rnd_seed,signal_length, Fs_neural, CWT_freq_limits,thresh_for_corr_inclusion)

    % Generate inter-scale correlation matrices from white noise input, 
    % used in the statistical significance test. This function uses
    % Mersenne Twisters to create the random numbers, using MATLAB's rand
    % function.
    
    % Inputs:
    % - rnd_seed: RNG seed   
    % - signal:length: how long the generated WN sequence is
    %    
    % Outputs: 
    % WN_ISCM_individual: inter-scale correlation matrix, produced from random white noise.
    % f: frequency vector.
    
    rng(rnd_seed);
    x_noise = rand(signal_length,1);
    
    %% Clip x_noise
    x_noise = rescale(x_noise,0,1);
    

    %% Calculate Inter-Scale Correlation Matrix (ISCM)
    [WN_ISCM_individual,f,computation_times] = calculate_inter_scale_correlation_matrix(x_noise, Fs_neural, CWT_freq_limits,thresh_for_corr_inclusion);
    
%     %% Compress ISCM, take advantage of symmetrical nature
%     [WN_ISCM_individual] = reshape_2D_symmetric_matrix_to_1D(WN_ISCM_individual);
%     if ~isequal(length(WN_ISCM_individual),(p^2-p)/2)
%         WN_ISCM_individual = zeros((p^2-p)/2,1);
%     end


end
