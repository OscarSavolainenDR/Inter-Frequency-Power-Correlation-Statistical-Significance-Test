function [x, n, computation_times, neural_time] = shorten_clip_and_prewhiten(x, neural_time)

    % Standardizes the length of, clips and pre-whitens the SI_CSD_Broadband
    % data.
    % 
    % Inputs: 
    % - SI_CSD_Broadband: the SI CSD referenced Broadband neural data.
    % - channel: the ID of the channel we are analysing, between 1 and 96
    % for the Sabes lab data.
    % - neural_time: the neural time stamps.
    %
    % Outputs:
    % - x: the neural time series from the chosen channel, after all of the
    % rest of the pre-processing (standardizing of length, clipping,
    % pre-whitening).
    % - n: the number of samples in x after the rest of the pre-processing.
    % - computation_times: how long certain operations took, for those who
    % are interested.
    % - neural_time: the neural time stamps.
    
    Inter_Scale_Stat_Sig_Sabes_Parameters; % initiliases the parameters used throughout the work

    %% Standardize the data length
    length_x_samples = length(x)/Fs_neural;
    shortened_index = get_standardised_length(length_MC_vector, minimum_time_gap, length_x_samples); % returns the closest smallest value in the vector of standardised lengths
    x(round(Fs_neural*(shortened_index+minimum_time_gap)):end) = []; % shorten

    % Have the recording have an even number of samples, which marginally simplifies the
    % Fourier calculation
    cutoff = length(x)-1;
    if rem(cutoff,2) == 0
        cutoff = cutoff + 1;
    end
    x = x(1:cutoff);
    
    %% Clip the end of the time series
    hh = tic;
    x = rescale(x,0,1);
    n = find_matching_end(x,clipping_epsilon); 
    x = x(1:n);
    neural_time = neural_time(1:n);
    
    %% Pre-whiten the time series
    x = whitening(x, Fs_neural); % pre-whiten data

    whitening_time = toc(hh);
    fprintf(['Whitening time was ',num2str(whitening_time),' s \n'])
    computation_times.whitening_time = whitening_time;
   
end
