% Analyses the publicly available results from doi: 10.21203/rs.3.rs-329644/v1
%
% Tests the significance of the inter-frequency correlation 
% matrices and plots them, plots the Fourier plots, and plots the time domain
% plots. This script also plots some meta-statistical information about 
% how many recordings had how many significant inter-frequency correlations.
% 
% NOTE 1: The pubicly available results need to be downloaded from
% [url: https://zenodo.org/deposit/4399993#]  or [doi: 10.5281/zenodo.4399993].
% Alternatively they can be computed from the original '.nwb' Sabes lab files,
% using the Python code from directory "Formatting_the_NWB_data_into_MATLAB", and then running 
% "main_recreate_results_from_mat_Sabes_data.m". This may take a few days,
% even a week, and requires a large amount of RAM (10s of GB).
% Computing the MC null distributions also likely requires a High
% Performance Computing (HPC) cluster, as hundreds/thousands of
% white noise processes of a large amount of samples need to be generated
% and their inter-scale correlation matrices (ISCM) calculated, along with
% the ISCMs derived from the phase-randomised WPS.
% 
% NOTE 2: If downloaded, the .zip files should be placed in distinct folders,
% and then unzipped. After unzipping, the .zip files can be deleted.
% The relevant directories then need to be filled in in 
% "Inter_Scale_Stat_Sig_Sabes_Parameters.m". 
% The relevant directories include:
% - "raw_neural_data_folder" - where the raw Sabes lab .mat files are (if
% downloaded results, can enter random non-empty string).
% - "FT_plots_CSD_Broadband_folder_midrange" - where to save Fourier plots,
% midrange, to specify.
% - "FT_plots_CSD_Broadband_folder_broadband" - where to save Fourier
% plots, broadband, to specify.
% - "time_domain_plots" - where to save time domain plots, to specify
% - "processed_data_folder" - where to save the SI-CSD referenced neural data
% (if downloaded results, can enter random non-empty string).
% - "save_white_noise_ISCMs" - directory with downloaded white noise
% generated ISCMs, where one placed the unzipped data.
% - "save_FT_phase_ran_ISCMs" - directory with downloaded ISCMs generated from
%  FT phase-ransomised neural WPS, where one placed the unzipped data.
% - "save_IAAFT_phase_ran_ISCMs"  - directory with downloaded ISCMs generated from
%  IAAFT phase-ransomised neural WPS, where one placed the unzipped data.
% - "Figures_stat_sig_inter_scale_corr_mats" -  where to save the ISCM
% plots, to specify.
% - "random_numbers_folder" - folder from random numbers used for white noise
% (if downloaded results, can enter random non-empty string).

clearvars
close all

Inter_Scale_Stat_Sig_Sabes_Parameters;

%% Parameters
% Analysed sessions and channels
sessions_vector = 1:15; % can manually edit, which sessions to analyse. min: 1; max: 15; for the publicly available results.
channels_vector = 1:20; % can manually edit, which channels to analyse. min: 1; max: 96; for the publicly available results.


% NOTE: the sessions_vector and channel_vecor can be different for each
% called function, e.g. for the ISCM plotting and the meta-stats plotting

% Whether to plot the Fouriers and time series plots in a 5x3 grid of 
% subplots (true), or individually, each with their own figure (false).
grid_bool = true; 

% Paramaters for plotting and analyzing the results
meta_parameters.show_figures = false;
meta_parameters.try_load_C_flag = true; % if it fails to load each specified neural ISCM, it continues if this flag is true. If false, it stops when it fails. False is better for testing, true is better for skipping over missing results.
meta_parameters.try_plot_flag = true; % if it fails to plot the ISCM, then if 'true', it continues. If 'false', it doesn't.
meta_parameters.save_stat_sig_inter_scale_corr_plots = false; % whether to save the plots of the statistically significant ISCMs to folder "save_stat_sig_inter_scale_plots", specified in "Inter_Scale_Stat_Sig_Sabes_Parameters".
meta_parameters.local_colorbar_bool = false; % whether to plot the colorbar individually for each ISCM
meta_parameters.plot_abs_corr = false; % whether to plot the correlations (false) or the absolute values of the correlations (true). The absolute values are less informative, but make the color scheme more intuitive.
meta_parameters.colormap_scheme = jet; % what MATLAB colormap scheme to use
meta_parameters.second_colormap_scheme = copper; % what MATLAB colormap scheme to use for meta-percentage plot

%% Merge FT phase-ransomised ISCM folders
% Note: one should first unzip the FT phase-ran files, and the path in
% which they were stored as unzipped folders should be specified in the
% 'save_FT_phase_ran_ISCMs' variable in the
% "Inter_Scale_Stat_Sig_Sabes_Parameters script". This function will delete
% the separated folders after merging.
if strcmp(use_IAAFT_or_FT_phase_ran_nulls,'FT')
    merge_FT_phase_ran_folders(sessions_vector,channels_vector);
end


%% Analyse results
% Determine and plot the statisically significant elements of the neural
% inter-scale corrrelation matrices (remove close all at the beginning of
% this script if necessary)
plot_sig_intra_channel_ISCMs(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters);

% Determine and plot the meta-statistical information, across all analysed recordings,
% of the statisically significant elements of the neural inter-scale correlation 
% matrices (remove close all at the beginning of this script if necessary)
plot_recording_meta_stats(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters);


%% Plot neural data (Time and Fourier plots)
% Plot the SI-CSD processed neural signals (Fourier and time domains)
% Plot the Fouriers of the SI-CSD Broadband. Plotting a subset is 
% generally easier, or just looking at the image files directly in whatever
% file explorer app is used in your computer.
example_sessions = 1;
example_channels = 1:2;
plot_Fouriers_all(example_sessions,example_channels);

% Plot the time series of the SI-CSD Broadband. Plotting a subset is 
% generally easier, or just looking at the image files directly in whatever
% file explorer app is used in your computer.
plot_time_series_all(example_sessions,example_channels, grid_bool);