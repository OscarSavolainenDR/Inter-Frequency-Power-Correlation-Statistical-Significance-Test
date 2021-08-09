%% Calculate Multiple testing threshold
function [H0,T,thresh] = multiple_testing_procedure_Cai_Liu(C,C_noises,f,fdr_alpha,thresh_vector,G,show_figures)

%     Function that performs the mutiple testing, given the genuine and Monte
%     Carlo inter-scale correlation matrices.
% 
%     Inputs:
%     - C: the 2D neural inter-scale correlation matrix.
%     - C_noises: the 3D null distribution inter-scale correlation matrix.
%     - f: the vector of frequencies.
%     - fdr_alpha: the proportion at which the FDR is controlled, between 0 and 1.
%     - thresh_vector: a vector containing the values at which the number of
%       rejected hypotheses is calculated. This is initialised in 
%       'Inter_Scale_Stat_Sig_Sabes_Parameters.m' so it does not have to be
%       computed each time this function is called.
%     - G: tail of normal distribution used in thresholding calculation,
%       see Ci and Liu (2017) for details.
%     - show_figures: a boolean, whether the map of rejected hypotheses is
%       plotted or not
%
%     Outputs:
%     - H0: the binary map of rejected hypotheses.
%     - T: the matrix of test statistics.
%     - thresh: the threshold above which the hypotheses were rejected.

    if ~isequal(size(C_noises),[length(f) length(f) length(C_noises(1,1,:))])
        error('C_noises should be of dimensions [length(f) length(f) nb_of_MCs]')
    end

    %% Initialise
    p = length(f);
    estimated_nulls = (p^2-p)/2;
    
    %% Calculate test statistic by standardizing C by C_noises
    T = (C -mean(C_noises,3))./std(C_noises,0,3);
    T(isnan(T)) = 0;
    
    %% Calculating threshold range, see Cai (2017)
    eye_mat = zeros(p);
    for i = 1:p
        eye_mat(i:p,i) = 1;
    end
    T(~logical(eye_mat)) = 0;
    R = zeros(length(thresh_vector),1);
    for i = 1:length(thresh_vector)
        R(i,1) = sum(sum(abs(T)>thresh_vector(i))); % rejected hypotheses for each threshold level
        if R(i,1) == 0
            R(i,1) = 1;
        end
    end
   
    %% Calculate threshold that controls the FDR at the given level
    temp = G.*estimated_nulls./R' < fdr_alpha;
    thresh = thresh_vector(min(find(temp==1)));
    if isempty(thresh)
        thresh = 2*sqrt(log(p));
    end

    %% Reject hypotheses above the threshold
    H0 = abs(T)>thresh; % get rejected hypotheses map (one-sided)
    H0 = H0 + H0'; % mirror rejected hypotheses map
    H0 = H0>0; % returns logical array

    %% Plot map of statistically significant correlations
%     if show_figures
%         figure
%         plot_2D_colormap_nan_ignore(single(H0),f,f,'log')
%         title('Rejected H0')
%         xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)')
%     end
end
