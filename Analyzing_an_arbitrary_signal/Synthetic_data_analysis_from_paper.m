%% Synthetic data ISPCM analysis
% This script can analyse multiple signals in matrix x_all, where each column
% represents its own signal. Here we generate and then analyse an envelope signal,
% double sinusoid, a double pulse, a double chirp, a square wave, and 
% white noise. The ISPCM null distributions are produced for each signal. 
% Any arbitrary signal can be tested as a column of x_all, or 
% the script could be adapted to allow signals of different lengths to be
% tested at the same time. It tends to run slowly for signals larger than
% 1e5 samples, in which case using a High Performance Cluster is
% recommended. 

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


all_r = cell(4,1);
all_f = cell(4,1);
signal_index = 1;

thresh_for_corr_inclusion = 0.1;
noise_level = 0.1;
clip_epsilon = 0.01;

FDR_alpha = 0.001;
max_MCs_per_channel_WN = 250; % larger values give more accurate estimates of the mean white noise ISPCM
max_MCs_per_channel_PR = 150; % larger values give more accurate estimates of the variance of the null distributions


%% Envelope
x_env = zeros(L,1);

% 1st envelope
env_1 = hanning(t1_lims(2)-t1_lims(1)+1);

% 2nd envelope
env_2 = hanning(t2_lims(2)-t2_lims(1)+1);

% Set
x_env(t1_lims(1):t1_lims(2)) = env_1;
x_env(t2_lims(1):t2_lims(2)) = env_2;
x_env = x_env + rand(L,1) *noise_level;
x_env = x_env - mean(x_env);
title_str{signal_index} = 'Envelope';

[r,f] = arbitrary_plots(t,x_env,Fs, title_str{signal_index}, thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;

% Fourier plot
figure, line_FFT_xlim(x_env,Fs,[0 Fs/2]);
title(['Frequency-domain ',title_str{signal_index}])

signal_index = signal_index + 1;

%% Pulses
x_pulse = zeros(L,1);

% Carrier freq
% Pulses
f1_1 = 30;
f1_2 = 5;

x_pulse(t1_lims(1):t1_lims(2)) = x_env(t1_lims(1):t1_lims(2)) .* sin(2*pi*f1_1*(t1-t(t1_lims(1))))';
x_pulse(t2_lims(1):t2_lims(2)) = x_env(t2_lims(1):t2_lims(2)) .* sin(2*pi*f1_2*(t2-t(t2_lims(1))))';
x_pulse = x_pulse + rand(L,1) *noise_level;
x_pulse = x_pulse -mean(x_pulse);

title_str{signal_index} = 'Double pulse';
[r,f] = arbitrary_plots(t,x_pulse,Fs, title_str{signal_index},thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;
signal_index = signal_index + 1;


%% Chirp
x_chirp = zeros(L,1);
% Carrier freq
% 1st Chirp
f1_1 = 50;
f0_1 = 20; % (beginning freq)
% 2nd Chirp 
f1_2 = 5;
f0_2 = 1; % (beginning freq)

% 1st chirp
t1_1 = t(t1_lims(2))-t(t1_lims(1));
y1 = chirp(t1-t1(1),f0_1,t1_1,f1_1,'quadratic',90);
% 2nd chirp
t1_2 = t(t2_lims(2))-t(t2_lims(1));
y2 = chirp(t2-t2(1),f0_2,t1_2,f1_2,'quadratic',90);

% Assign
x_chirp(t1_lims(1):t1_lims(2)) = x_env(t1_lims(1):t1_lims(2)) .* y1';
x_chirp(t2_lims(1):t2_lims(2)) = x_env(t2_lims(1):t2_lims(2)) .* y2';
x_chirp = x_chirp + rand(L,1) * noise_level;
x_chirp = x_chirp -mean(x_chirp);

title_str{signal_index} = 'Chirp';
[r,f] = arbitrary_plots(t,x_chirp,Fs, title_str{signal_index},thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;
signal_index = signal_index + 1;


%% Square waves
x_square = zeros(L,1);

% Combined square wave
f1 = 2*pi*0.005; % 0.01
t_combined = t(L/4+100:L/4+100 + Fs/f1 * 1); % 2

y_comb = square(2*pi*f1*(t_combined - t_combined(1)));

% Assign
x_square(L/4+100:L/4+100 + Fs/f1 * 1) = y_comb;
x_square = x_square + rand(L,1) * noise_level;
x_square = x_square - mean(x_square);

title_str{signal_index} = 'Square wave';
[r,f] = arbitrary_plots(t,x_square,Fs, title_str{signal_index}, thresh_for_corr_inclusion, clip_epsilon);
all_r{signal_index} = r;
all_f{signal_index} = f;
signal_index = signal_index + 1; 

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
x_sin = x_sin -mean(x_sin);

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
x_all = [x_env, x_pulse, x_chirp, x_square, x_sin, x_white]; %

% fprintf('\n')
% disp('Press a key to determine statistical significance of ISPCMs')
% pause;

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
    p = length(f); % THIS SHOULD BE FILLED IN MANUALLY AND WILL VARY FOR DIFFERENT APPLICATIONS. Number of scales for given CWT_freq_limits, can be discovered empirically for each time series (depends on length) and CWT decomposition, i.e. length of frequency vector returned by CWT, length of ISCM matrix (1 side).
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
