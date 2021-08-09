function [r,f,computation_times] = calculate_inter_scale_correlation_matrix(x, Fs, freq_limits, thresh_for_corr_inclusion)
    
    % Calculates the inter-scale correlation matrices from a given time
    % series.
    %
    % Inputs: 
    % - x: the time series.
    % - Fs: the sampling frequency of x.
    % - freq_limits: a 2-element vector containing the lower and upper
    % limits, in that order, of the CWT operation, e.g. return the
    % coefficients between 0.1 and 1000 Hz.
    % - thresh_for_corr_inclusion: a value between 0 and 1, that specifies
    % the proportion of the number of samples in x that at the least need
    % to be involved in every inter-scale correlation calculation. Because
    % some values are outside of the Cone of Influence, some scales have
    % fewer reliable values that others. If a scale contains fewer than 
    % thresh_for_corr_inclusion * length(x) values inside the COI, then it
    % is not included in the inter-scale correlation matrix. The choice of
    % this value will depend largely on each user's desired number of
    % samples per inter-scale correlation calculation.
    %
    % Outputs:
    % - r: inter-scale correlation matrix of x.
    % - f: the vector of frequencies that are involved in the inter-scle
    % correlation matrix.
    % - computation_times: contains the computation times of each process,
    % for those who are interested.
    %
    % This function can be used for any time series, given the sampling
    % rate and desired frequency limits of the Wavelet Power Spectrum.
    % 
    % All code is available at XXX.
    

    %% Perform Continuous Wavelet Transform (CWT) on x
    hh = tic;
    [S,f,coi] = cwt(x,'morse',Fs,'FrequencyLimits',freq_limits); % S: Wavelet Coefficients, f: scale vector, coi : Cone of Influence.
    perform_1_CWT_time = toc(hh);
    fprintf(['Time to do 1 CWT was ',num2str(perform_1_CWT_time),' s \n'])
    hh = tic;
    S = single(S);
    

    %% Making eveything outside COI zero
    S = S';
    rotation_time = toc(hh);
    fprintf(['Rotation time was ',num2str(rotation_time),' s \n'])
    hh = tic;

    g = length(f);
    for i = 1:length(coi)
        [ ~, ix ] = min( abs( f-coi(i) ) );
        S(i,ix:g)= 0;
    end
    processing_time = toc(hh);
    fprintf(['Processing time was ',num2str(processing_time),' s \n'])
    hh = tic;

    %% Calculate Wavelet Power Spectrum
    S = abs(S).^2;
    S(S==0)= NaN; % set zero values to NaN (including all values outside the COI).
    abs_nan_time = toc(hh);
    fprintf(['Abs Nan time was ',num2str(abs_nan_time),' s \n'])
    hh = tic;

    %% Get inter-scale correlations of S -> C (in paper it is 'r')
    % Note: Does not calculate the correlations that had fewer than
    % (thresh_for_corr_inclusion x 100)% of the data points in x in their
    % calculation. This is only to insure that enough data points are
    % included in each inter-scale correlation calculation.
    temp = isnan(S); % whether a value is NaN  or not
    ind = sum(temp); % the number of NaN samples per scale
    f_nan_ratio = ind/length(x); % ratio of NaN to not-NaN points for each scale
    a = f_nan_ratio < thresh_for_corr_inclusion; % 1 = that scale has enough data points, 0 = doesn't.
    f(~a) = []; % remove frequencies from the vector that don't have enough samples in their calculation.
    S(:,~a) = [];
    
    % Cut off all of the NaN values, makes the corrcoef calculation way, way, way easier (faster) than using the 'pairwise' option in corrcoef
    [~, min_snip_cutoff] = min(temp);
    S(1:min_snip_cutoff(max(cumsum(a))),:) = [];
    [~, max_snip_cutoff] = max(isnan(S));
    S(max_snip_cutoff(max(cumsum(a))):end,:) = [];
    clear temp
    
    fprintf(['Signal length reduced from ',num2str(length(x)),' to ',num2str(length(S)),' samples \n'])

    processing_time = toc(hh);
    fprintf(['Processing time was ',num2str(processing_time),' s \n'])
    
    hh = tic;
    r = corrcoef(S); % Gets inter-scale correlations, ignores NaN rows (e.g. ignores samples outside the COI).
    correlation_time = toc(hh);
    fprintf(['Correlation time was ',num2str(correlation_time),' s \n'])

    % Store computation times
    computation_times.perform_1_CWT_time = perform_1_CWT_time;
    computation_times.rotation_time = rotation_time;
    computation_times.processing_time= processing_time;
    computation_times.abs_nan_time = abs_nan_time;
    computation_times.correlation_time = correlation_time;

end
