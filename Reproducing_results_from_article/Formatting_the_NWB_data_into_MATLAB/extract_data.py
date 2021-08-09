# Function for extracting data (the specified features) from .nwb files

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

def extract_data(raw_filename_inside):
    # Function for extracting data (the specified features) from .nwb files
    import h5py
    with h5py.File(raw_filename_inside, 'r') as f:
        nwb_time_inside = f['/acquisition/timeseries/broadband/timestamps'][:]
        electrode_idx_inside = f['/acquisition/timeseries/broadband/electrode_idx'][:]
        electrode_map_inside = f['/general/extracellular_ephys/electrode_map'][:]
        raw_data_inside = f['/acquisition/timeseries/broadband/data'][:]

    return nwb_time_inside, electrode_idx_inside, electrode_map_inside, raw_data_inside
