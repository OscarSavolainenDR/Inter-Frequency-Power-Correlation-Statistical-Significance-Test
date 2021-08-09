% Recreate the results from The Significance of Neural Inter-Frequency Correlations
% doi: 10.21203/rs.3.rs-329644/v1
%
% Ths script processes the raw broadband neural data, produces the inter-frequency correlation 
% matrices, tests their significance, and presents the results. It also 
% saves plots of the processed neural data (Fourier and time domain
% plots), and plots them.
%
% NOTE 1: Testing the significance requires the use of null inter-scale correlation
% null distributions derived from Monte Carlo simulations. The BASH code for
% the white-noise generated ISCMs is contained in folder
% "MC_Sim_for_Producing_Null_Distributions" in file
% "main_WN_ISCM_MC_simulation.pbs".
% If using FT phase-ransomisation, the BASH code for the FT 
% phase-randomisation derived ISCMs is in file
% "main_FT_ISCM_MC_simulation.pbs".
% If using IAAFT phase-ransomisation, the BASH code for the IAAFT 
% phase-randomisation derived ISCMs is in file
% "main_IAAFT_ISCM_MC_simulation.pbs".
%
% NOTE 2: One should also first run the Python code in 
% "Formatting_the_NWB_data_into_MATLAB" to format the .nwb files as .mat
% files. See main directory README.md
%
% NOTE 3: Some directories need to be filled in in the 
% 'Inter_Scale_Stat_Sig_Sabes_Parameters' script. The relevant directories include:
% - "raw_neural_data_folder" - where the raw Sabes lab .mat files were
% saved during the Python script .nwb-to-.mat formatting.
% - "FT_plots_CSD_Broadband_folder_midrange" - where to save Fourier plots,
% midrange frequencies.
% - "FT_plots_CSD_Broadband_folder_broadband" - where to save Fourier
% plots, broadband.
% - "time_domain_plots" - where to save SI-CSD processed neural time domain
% plots.
% - "processed_data_folder" - where to save the SI-CSD referenced neural data.
% - "save_white_noise_ISCMs" - where to save the white noise generated ISCMs.
% - "save_FT_phase_ran_ISCMs" - where to save the ISCMs generated from
%  FT phase-ransomised neural WPS, where one placed the unzipped data.
% - "save_IAAFT_phase_ran_ISCMs"  - where to save the ISCMs generated from
%  IAAFT phase-ransomised neural WPS, where one placed the unzipped data.
% - "Figures_stat_sig_inter_scale_corr_mats" - where to save the ISCM
% plots after stat testing.
% - "random_numbers_folder" - folder for random numbers used for WN ISCM
% generation, probably a HPC directory (the C code for RNG runs on Linux).
close all
clearvars

Inter_Scale_Stat_Sig_Sabes_Parameters;

%% Parameters for processing data, analysing results and plotting neural signals
% Analysed sessions and channels
sessions_vector = 1:15; % can manually edit, which sessions to analyse. min: 1; max: 15; for the publicly available results.
channels_vector = 1:20; % can manually edit, which channels to analyse. min: 1; max: 96; for the publicly available results.

% Parameters for plotting and analyzing the results
meta_parameters.show_figures = false;
meta_parameters.try_load_C_flag = true; % if it fails to load each specified neural ISCM, it continues if this flag is true. If false, it stops when it fails. False is better for testing, true is better for skipping over missing results.
meta_parameters.try_plot_flag = true; % if it fails to plot the ISCM, then if 'true', it continues. If 'false', it doesn't.
meta_parameters.save_stat_sig_inter_scale_corr_plots = true; % whether to save the plots of the statistically significant ISCMs to folder "save_stat_sig_inter_scale_plots", specified in "Inter_Scale_Stat_Sig_Sabes_Parameters".
meta_parameters.local_colorbar_bool = false; % whether to plot the colorbar individually for each ISCM
meta_parameters.plot_abs_corr = false; % whether to plot the correlations (false) or the absolute values of the correlations (true). The absolute values are less informative, but make the color scheme more intuitive.
meta_parameters.colormap_scheme = jet; % what MATLAB colormap scheme to use
meta_parameters.second_colormap_scheme = copper; % what MATLAB colormap scheme to use for meta-percentage plot

% Whether to plot the Fouriers and time series plots in a 5x3 grid of 
% subplots (true), or individually, each with their own figure (false).
grid_bool = false;

% The frequency ranges for the Fourier plots. It is generally useful to
% have 2 plots at different ranges, to see broadband and narrowband
% Fouriers.
% NOTE: The user can customize the below ranges as desired.
broadband_range = [30 8600]; % Hz, plots the Fourier in this range, from 1st element of the vector to the 2nd
mid_range = [20 500]; % Hz, also plots the Fourier in this range, from 1st element of the vector to the 2nd

%% Process the data
% Pre-process the neural data
% NOTE: One has the option to download the associated behavioral data
% (from https://zenodo.org/record/3854034#.YFzdL6_7Q2x [main files, not 
% Supplemental files]). Selecting "sync_neural_and_behavioral = true" syncs
% the neural and behavioral data, if one wants to do some behavioral
% decoding. If one is only looking at the neural data,
% "sync_neural_and_behavioral = false" is fine and saves memory and time.
sync_neural_and_behavioral = false;
main_SI_CSD_syncing(sessions_vector, sync_neural_and_behavioral);

% Calculate the null distributions. This involves transforming
% "main_WN_ISCM_MC_simulation.txt" into a PBS file, and running it in a HPC
% cluster. Or, code and run a script that calls function "generate_WN_ISCMs"
% appropriately to produce a sufficient number of WN ISCMs.
% Then, one needs to choose to use either FT or IAAFT phase-randomisation derived
% ISCMs. If FT, this involves transforming
% "main_WN_ISCM_MC_simulation.txt" into a PBS file, and running it in a HPC
% cluster. Or, code and run a script that calls function "generate_FT_ISCMs"
% appropriately to produce a sufficient number of FT ISCMs. If IAAFT, do
% the same for main_WN_ISCM_MC_simulation.txt / generate_IAAFT_ISCMs.
% Then, delete or uncomment the error statement below.
error('Produce the null distributions using the code form directory:\n"MC_Sim_for_Producing_Null_Distributions". If one\nis using the FT phase-randomisation method, set use_IAAFT_or_FT_phase_ran_nulls = "FT",\notherwise set it as "IAAFT". use_IAAFT_or_FT_phase_ran_nulls is located in Inter_Scale_Stat_Sig_Sabes_Parameters')


%% Plot and save significant inter-frequency correlation results
% Determine and plot the statisically significant elements of the neural
% inter-scale corrrelation matrices (remove close all at the beginning of
% this script if necessary)
plot_sig_intra_channel_ISCMs(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters);

% Determine and plot the meta-statistical information, across all analysed recordings,
% of the statisically significant elements of the neural inter-scale correlation 
% matrices (remove close all at the beginning of this script if necessary)
plot_recording_meta_stats(sessions_vector,channels_vector,use_IAAFT_or_FT_phase_ran_nulls,meta_parameters);


%% Save the SI-CSD processed neural signals (Fourier and time domains)
% Save the Fouriers of the SI-CSD processed Broadband data, for 2
% ranges, as it's useful to see both
save_fouriers_of_processed_data(sessions_vector,channels_vector,broadband_range,mid_range,grid_bool);

% Save the time series of the SI-CSD processed Broadband data
save_time_series_of_processed_data(sessions_vector,channels_vector,grid_bool);

%% Plot the SI-CSD processed neural signals (Fourier and time domains)
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

