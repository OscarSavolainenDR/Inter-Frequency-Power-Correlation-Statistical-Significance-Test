
# Script for extracting data from .nwb files and formatting it to .mat files, for use 
# with the Sabes lab raw broadband (Supplemental) neural dataset from
# https://zenodo.org/record/3854034#.YFyjfa_7Q2w.

##################################################################################
# MIT License
#
# Copyright (c) 2020 OscarSavolainen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###################################################################################

import h5py
from os import listdir
from extract_data import extract_data
import scipy.io


# Directories to data (need to be filled in)
Sabes_folder = r"path to folder with .nwb files"
results_folder_mat = r"path to folder where you want the .mat files"
x = listdir(Sabes_folder)


# Using for loop to extract.nwb data and save as .mat file
# raw_filename = 'indy_20160930_02.nwb' # example
for raw_filename in x:
    print("Formatting raw neural data from file: " + raw_filename)
    raw_filename_total = Sabes_folder + '\\' + raw_filename  # path to file
    (nwb_time, electrode_idx, electrode_map, raw_data) = extract_data(raw_filename_total)

    stored_data = {
        "time": nwb_time,
        "electrode_ID": electrode_idx,
        "electrode_Map": electrode_map,
        "neural_recordings": raw_data
    }

    scipy.io.savemat(results_folder_mat + "\\" + raw_filename.replace('.nwb', '') + ".mat",
                     mdict={'stored_data': stored_data})

# # Example plot
# import matplotlib.pyplot as plt
# num_steps = 100000
# y = stored_data.neural_recording[1:num_steps,1]
# t = stored_data.time[1:num_steps]
#
# plt.figure(0)
# plt.plot(t,y)
# plt.show()


