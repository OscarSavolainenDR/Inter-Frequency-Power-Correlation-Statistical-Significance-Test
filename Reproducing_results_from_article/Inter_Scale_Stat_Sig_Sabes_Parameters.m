%% Short script for saving hyper-parameters for the whole project
% Contains:
% - File directories;
% - Pre-processing parameters;
% - Post-processing parameters;
% - Multiple testing parameters;
% - Utah array map.

% To access the parameters, simply run
% "Inter_Scale_Stat_Sig_Sabes_Parameters;" in the script or command line.

%% Select either FT or IAAFT phase-randomisation
use_IAAFT_or_FT_phase_ran_nulls = ''; % 'FT' or 'IAAFT'

%% Directories
% Raw data
% Note: if downloaded results, this can be any random non-empty string.
raw_neural_data_folder = '';
behavioral_data_folder = ''; % If we're not interested in syncing the behvioral and neural data, any random string can be entered here

% Processed data
% Note: if downloaded results, this can be any random non-empty string.
processed_data_folder = '';

% Fourier plots. One could just have 1 folder, but it's useful to have 2 or more.
FT_plots_CSD_Broadband_folder_midrange = '';
FT_plots_CSD_Broadband_folder_broadband = '';

% Time domain plot
time_domain_plots = '';

% Inter-scale corrrelation matrices, and plots
save_white_noise_ISCMs = '';
save_FT_phase_ran_ISCMs = '';
save_IAAFT_phase_ran_ISCMs = '';
Figures_stat_sig_ISCMs = '';

% Random white noise 
% Note: if downloaded results, this can be any random non-empty string.
random_numbers_folder = ''; % folder to save PCG random numbers to, during white-noise IFCM generation. It should be the same folder with the PCG C functions.

% Test to make sure they are all filled in
if isempty(processed_data_folder) || ...
        isempty(raw_neural_data_folder) || ...
        isempty(behavioral_data_folder) || ...
        isempty(FT_plots_CSD_Broadband_folder_midrange) || ...
        isempty(FT_plots_CSD_Broadband_folder_broadband) || ...
        isempty(time_domain_plots) || ...
        isempty(save_white_noise_ISCMs) || ...
        isempty(save_FT_phase_ran_ISCMs) || isempty(save_IAAFT_phase_ran_ISCMs) || ...
        isempty(Figures_stat_sig_ISCMs) || isempty(random_numbers_folder)
    error('Ensure that all of the directories in "Inter_Scale_Stat_Sig_Sabes_Parameters" have been filled in.')
end

%% Processed filenames
% The processed files in this work
processed_sessions = {'indy_20160624_03',...
                      'indy_20160915_01',...
                      'indy_20160916_01',...
                      'indy_20160921_01',...
                      'indy_20160927_04',... % 5
                      'indy_20160927_06',...
                      'indy_20160930_02',...
                      'indy_20160930_05',...
                      'indy_20161005_06',...
                      'indy_20161006_02',... % 10
                      'indy_20161007_02',...
                      'indy_20161011_03',...
                      'indy_20161013_03',...
                      'indy_20161014_04',...
                      'indy_20161017_02'};   % 15
       
% Analysed sessions and channels
sessions = 1:15; % which sessions to analyse
channels = 1:96; % which channels to analyse


%% Pre-processing parameters used in this work
Fs_behav = 250; % Behavioral data sampling frequency
Fs_neural = 24414.0625; % Neural data sampling rate
intro_cutoff = 1e5; % Cut off the beginning of neural signal
figures_cutoff = 10000; % for plotting figures during CSD referencing
% raw_Broadband = double(stored_data.neural_recordings(intro_cutoff:end,:));

% Spectral Interpolation
sample_channel = 10; % in the Sabes dataset this channel has strong harmonics, can be used to zero-in on line frequency across the recording
harmonic = 1; % 150; % the bigger the harmonic the better the estimation of the line noise, however there is increased risk of 1) missing the harmonic and 2) it may not be there in the first place
harmonic_width = 2; % Hz, how wide a section to interpolate across for each line noise harmonic
harmonic_neighbours = 3; % Hz, how many Hz to use in calculating the edge values of the interpolation. Practically, this value is unimportant, as the data is later pre-whitened so the amplitude of the interpolation is insignificant.
component_neigbours = 3; % scalar value, how many frequencies to use in calculating the average value on either side of the interpolation in the single component interpolation. Practically, this value is unimportant, as the data is later pre-whitened so the amplitude of the interpolation is insignificant.
freq_res_peak_search = 5; % how far to look on either side of the predicted harmonic location for the peak (e.g. expected at (harmonic_index x 60), but it may be at ((harmonic_index x 60) +- freq_res_peak_search).

% Lowpass frequency for extracting LFP during Current Source Denity (CSD) referencing
LFP_cutoff = 500; % Used in LFP-CSD offset. The mean LFP from all other electrodes is removed from the Broadband of each electrode.
% LFP = lowpass(raw_Broadband(Broadband,LFP_cutoff,Fs);

% Individual components to spectrally interpolate over. The cell
% represents the session. The first entry is the smallest 
% frequency of the first component. The second entry is the largest 
% frequency of the first component. The third entry, if non-zero, is the 
% smallest frequency of the second component, etc. All entries are in Hz.
component{1,1} = [2845 2880];
component{2,1} = [2610 2660 2320 2360];
component{3,1} = [2630 2690];
component{4,1} = [2630 2680];
component{5,1} = [2710 2750 2360 2390]; % 5
component{6,1} = [2730 2760 2400 2420];
component{7,1} = [2250 2350];
component{8,1} = [2400 2440 2730 2770];
component{9,1} = [2800 2830 2440 2470];
component{10,1} = [2700 2760 2350 2410]; % 10
component{11,1} = [2220 2290];
component{12,1} = [2640 2690 2320 2360];
component{13,1} = [2675 2720];
component{14,1} = [2760 2790 2410 2450];
component{15,1} = [2710 2750 2360 2400]; % 15

clipping_epsilon = 0.01; % when clipping the data (re-scaled to between 0 and 1), clip the end of x so that the end is within clipping_epsilon of x[1].

% According to notes: The (raw Broadband) data are sampled at 24414.0625 Hz and are unfiltered,
% except for an anti-aliasing filter built-in to the recording amplifier: a 4th order low-pass
% with a roll-off of 24 dB per octave at 7.5 kHz, operating at the sampling rate.

%% Post-processing parameters
CWT_freq_limits = [0.1 8.5e3];
thresh_for_corr_inclusion = 0.1; % (1-this value) is the minimum ratio of non-Nan points, relative to NaN points, it needs to have to be included in the R matrix
FDR_alpha = 0.01; % FDR alpha

length_MC_vector = [350 400 500]; % lengths in s of the MC white noise processes (standardised recording lengths)
minimum_time_gap = 5; % the elements in MC_length_vector should be seperated by at least this much (in s). Depends on the MC run, but this allows us to do some indexing elsewhere.

length_MC_vector = sort(length_MC_vector);
if sum(diff(length_MC_vector)<minimum_time_gap)>0
    error(['The elements in MC_length_vector need to be seperated by at least ',num2str(minimum_time_gap ),' s'])
end

%% Multiple testing parameters
p = 163; % THIS SHOULD BE FILLED IN MANUALLY AND WILL VARY FOR DIFFERENT APPLICATIONS. Number of scales for given CWT_freq_limits, can be discovered empirically for each time series (depends on length) and CWT decomposition, i.e. length of frequency vector returned by CWT, length of ISCM matrix (1 side).
ap = 2*log(log(p));
bp = sqrt(4*log(p)-ap);
Cai_Liu_thresh_vector = 0:0.001:bp;
G = 2-2*normcdf(Cai_Liu_thresh_vector);


%% Utah electrode map
% NOTE: Uncomment out the following block to visualise the Utah array with
% labelled electrodes

% figure, hold on
% for i = 1:max(size(stored_data.electrode_Map))
%     text(double(stored_data.electrode_Map(i,1)),double(stored_data.electrode_Map(i,2)),num2str(i))
% end
% xlim([min(stored_data.electrode_Map(:,1)) max(stored_data.electrode_Map(:,1))])
% ylim([min(stored_data.electrode_Map(:,2)) max(stored_data.electrode_Map(:,2))])
% title('Electrode map')
% electrode_map_matrix =  [0 45 88 84 80 76 72 68 60 0;
%                         47 6  2  8  86 82 78 64 56 52;
%                         41 43 21 23 4  66 74 70 71 67;
%                         39 37 17 19 7  11 62 50 69 63;
%                         35 33 13 15 91 3  54 58 65 59;
%                         31 29 9  5  87 83 95 57 61 55;
%                         25 27 1  93 85 79 75 96 92 53;
%                         46 48 44 89 81 77 73 94 90 49;
%                         42 40 36 32 28 24 20 16 12 51;
%                         0  38 34 30 26 22 18 14 10 0 ];

