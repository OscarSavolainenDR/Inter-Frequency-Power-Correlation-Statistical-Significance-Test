%% Arbitrary data ISPCM analysis
% This script can analyse multiple signals in matrix x_all, where each column
% represents its own signal. Any arbitrary signal can be tested as a column of x_all, or 
% the script could be adapted to allow signals of different lengths to be
% tested at the same time. It tends to run slowly for signals larger than
% 1e5 samples, in which case using a High Performance Cluster is
% recommended. 

% Two example signals are generated and analysed: a double sinusoidal
% signal enveloped using a Hanning window, and a white noise signal.

clearvars
% close all

%% Inititialise
L = 5e4;
Fs = 1e3;
t = 0/Fs:1/Fs:(L-1)/Fs;

% Temporal limits
t1_lims = round([0.2*L 0.21*L]);
t1 = t(t1_lims(1):t1_lims(2));
t2_lims = round([0.5*L 0.53*L]);
t2 = t(t2_lims(1):t2_lims(2));


all_r = cell(2,1);
all_f = cell(2,1);
signal_index = 1;

thresh_for_corr_inclusion = 0.1;
noise_level = 0.1;
clip_epsilon = 0.01;

FDR_alpha = 0.001;
max_MCs_per_channel_WN = 100; % larger values give more accurate estimates of the mean white noise ISPCM
max_MCs_per_channel_PR = 100; % larger values give more accurate estimates of the variance of the null distributions

%% Sinusoidal
f1_sin_harm = 50;
f2_sin_harm = 10;
x_sin = sin(2*pi*f1_sin_harm*(t)) + sin(2*pi*f2_sin_harm*(t));

env_big = hanning(round(L/50));
env_big = [zeros(round(24*L/50),1); env_big; zeros(10+round(25*L/50),1)];
if length(env_big) > L
    env_big(L+1:end) = [];
end

x_sin = x_sin' .* env_big + rand(L,1) * noise_level;
x_sin = x_sin - mean(x_sin);

title_str{signal_index} = 'Sinusoid harmonics';
[r,f] = arbitrary_plots(t,x_sin,Fs, title_str{signal_index}, thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;
signal_index = signal_index + 1;


%% White noise
x_white = rand(L,1) * noise_level;
x_white = x_white - mean(x_white);

% [x_white, n, computation_times, white_time] = arb_shorten_clip_and_prewhiten(x_white, t, Fs);

title_str{signal_index} = 'White noise';
[r,f] = arbitrary_plots(t,x_white,Fs, title_str{signal_index}, thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;
signal_index = signal_index + 1;


%% Collate all signals
x_all = [x_sin, x_white];

fprintf('\n')
disp('Press a key to determine statistical significance of ISPCMs')
pause;

%% Generate null distributions
[white_noise_ISCMs, all_f_WN] = arbitrary_get_WN_ISCM(x_all,Fs,max_MCs_per_channel_WN,thresh_for_corr_inclusion, clip_epsilon);
[all_r_phase, all_f_PR] = arbitrary_get_FT_phase_ran(x_all,Fs,max_MCs_per_channel_PR,thresh_for_corr_inclusion, clip_epsilon);

for i = 1:length(all_f)
    if length(all_f_WN{i}) ~= length(all_f_PR{i})
        fprintf('Something wrong, WN f length does not match PR f length \n')
        break
    end
    if length(all_f_WN{i}) ~= length(all_f{i})
        fprintf('Something wrong, WN f length does not match local f length \n')
        break
    end
end

%% Stat test
close all
for syn_signal_index = 1:length(x_all(1,:))
    phase_ran_ISCMs = all_r_phase{syn_signal_index};
    white_noise_ISCMs_used = white_noise_ISCMs{syn_signal_index};
    f = all_f{syn_signal_index};
    r = all_r{syn_signal_index};

    % Params
    p = length(f);
    ap = 2*log(log(p));
    bp = sqrt(4*log(p)-ap);
    Cai_Liu_thresh_vector = 0:0.001:bp;
    G = 2-2*normcdf(Cai_Liu_thresh_vector);
    show_figures = true;
    plot_abs_corr = false;

    % Get null distribution
    % Sum mean white-noise ISCM and phase-ran ISCM
    % distributions to get null distribution
    null_distributions = phase_ran_ISCMs + mean(white_noise_ISCMs_used,3) - eye(p); % deleting the identity matrix is optional, since diagonal elements are never considered in the stat test anyway.

    % Test the significance of the neural ISCM elements
    % Returns the indicator function of the significant relationships, the test
    % statistic matrix, and the threshold for significant results.
    hh = tic;
    [H,test_statistics,thresh] = multiple_testing_procedure_Cai_Liu(r, null_distributions,f,FDR_alpha,Cai_Liu_thresh_vector,G,show_figures);
    fprintf(['Stat test took ',num2str(toc(hh)),' s \n'])
    r(~logical(H)) = NaN; % set non-significant values to NaN
    
    if sum(sum(isnan(r))) ~= numel(r)
        [~,~] = arbitrary_plots(t,x_all(:,syn_signal_index),Fs, title_str{syn_signal_index},thresh_for_corr_inclusion, clip_epsilon);
        colorbar_bool = true;
        plot_abs_corr = false;
        caxes = [min(min(r)) max(max(r))];
        
        figure,
        plot_2D_colormap_nan_ignore(r, f, f, 'log',jet, colorbar_bool, plot_abs_corr, caxes)
    %     colorbar
        xlabel('Frequency (Hz)')
        ylabel('Frequency (Hz)')
        title(['ISPCM ',title_str{syn_signal_index}])
        xticks([0.1 1 10 100])
        yticks([0.1 1 10 100])
    else
        
        [~,~] = arbitrary_plots(t,x_all(:,syn_signal_index),Fs, title_str{syn_signal_index},thresh_for_corr_inclusion, clip_epsilon);
        colorbar_bool = true;
        plot_abs_corr = false;
        
        figure,
        plot_2D_colormap_nan_ignore(r, f, f, 'log',jet, colorbar_bool, plot_abs_corr, [-1 1])
    %     colorbar
        xlabel('Frequency (Hz)')
        ylabel('Frequency (Hz)')
        title(['ISPCM ',title_str{syn_signal_index}])
        xticks([0.1 1 10 100])
        yticks([0.1 1 10 100])
        title([title_str{syn_signal_index},' failed'])
    end
end
