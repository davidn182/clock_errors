#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:45:04 2021
Clock drift is a program to correct the clock errors of obs seismometers.
@author: davidnaranjo
"""

import numpy as np

import obspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.cross_correlation import correlate, xcorr_max

import os
import sys
import subprocess
from contextlib import contextmanager

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.linalg.lapack import dgelsd
import pandas as pd

mpl.style.use('bmh')

ext = os.path.dirname(os.path.abspath(__file__))
inv_ext = os.path.join(ext, 'inventory', 'inventory.xml')
#%%
# Some parameters that come as default but can be modified directly when 
# calling the functions.

program_dir = os.path.join(ext, 'recover_timing_errors-master/')
lf = 0.15 # Low freq. for the bandpass filter
hf = 0.3 # High freq. for the bandpass filter 
ref_vel = 2500 # m/s
dist_trh = 2.0 # Minimum station separation in terms of wavelength
snr_trh = 10 # Signal-to-noise ratio threshold
noise_st = 240 # start of the noise window.
apr_dt_st1 = 0. # A-priori estimate of station 1
apr_dt_st2 = 0. # A-priori estimate of station 1
dt_err = 0.004 # Sampling interval needs to be multiple of this value.
#%%
inv = obspy.read_inventory(inv_ext)

#%%
@contextmanager
def suppress_stdout():
    '''
    Function to hide the ouput in the terminal to make some of the processes 
    faster without the need of seing the terminal's output.

    Returns
    -------
    None.

    '''
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
            
def read_xcorrelations(station1, station2, directory=
                       '/Users/davidnaranjohernandez/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/'):
    '''
    Function to load the cross-correlations given the names of the station
    pair.

    Parameters
    ----------
    station1 : str
        Name of station 1.
    station2 : str
        Name of station 2.
    directory: str
        Directory where the cross-correlations are located.
    Returns
    -------
    None.

    '''
    if station1 == station2:
        print('No autocorrelations allowed')
        raise
    os.listdir(directory)

    # directory='/Users/davidnaranjohernandez/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/'
    
    files = os.listdir(directory)
    
    
    # Correlation files in a list to return.
    corr_st = obspy.Stream()
    corr_dirs = []
    for file in files:
        if station1 in file and station2 in file:
            correlations = os.listdir(directory + file)
            for corr_file in correlations:
                if '.SAC' in corr_file:
                    if 'total' in corr_file:
                        1+1
                        # total = obspy.read(os.path.join(directory, file,
                        # corr_file))
                    else:
                        corr_dir = os.path.join(directory, file, corr_file)
                        corr_tr = obspy.read(corr_dir)[0]
                        corr_dirs.append(corr_dir)
                        # The file header contains the day in the middle of the correlation
                        # with 50 days before and 50 days after the correlation starts.
                        center_date =  obspy.UTCDateTime(int(corr_file.split('_')[0]))
                        number_days = corr_file.split('_')[2].replace('.SAC', '')
                        number_days = float(number_days)*100
                        # start = (center_day - 50*24*60*60)
                        #end = (center_day + 50*24*60*60)
                        
                        # The file header contains the day in the middle of the correlation
                        # averaged over the availablle days
                        average_date =  obspy.UTCDateTime(int(corr_file.split('_')[1]))
                        corr_tr.stats.average_date = average_date
                        corr_tr.stats.number_of_days = number_days
                        corr_tr.stats.center_date = center_date
                        corr_tr.stats.corr_period_id = corr_file.split('_')[0]
                        corr_tr.stats.station_pair = station1 + '_' + station2
                        corr_st += corr_tr
    return corr_st, corr_dirs

def plotCCF(tr, maxlag, figurename=None):

    # define the time vector for the correlation (length of corr = corrwin + 1)
    limit = (len(tr.data) / 2.) * tr.stats.delta
    timevec = np.arange(-limit, limit, tr.stats.delta)

    plt.rcParams["figure.figsize"] = [10,5]
    plt.plot(timevec, tr.data, 'k')
    plt.xlim(-maxlag, maxlag)
    plt.xlabel('time [s]')
    plt.show()

def calculate_result_shift(data_dir, station_1, station_2, dir_name,
                           nsamples, lf, hf, ref_vel, dist_trh, snr_trh,
                           noise_st, apr_dt_st1, apr_dt_st2, dt_err,
                           inv=inv,
                           resp_details=False):
    '''
    Function that runs with fortran code located in program_dir to calculate
    the difference between the time of the negative and positive virtual
    source arrival. 
    The results will be saved in a directory called 
    '{program_dir}/temp/{station1}_{station2}_{dir_name}/'

    Parameters
    ----------
    data_dir : TYPE
        location of the correlation pairs.
    station_1 : TYPE
        name of station 1.
    station_2 : TYPE
        name of station 2.
    dir_name : TYPE
        The name of the directory where the results will be saved. I suggest
        to write here the date 
    nsamples : TYPE
        Number of samples in the trace.
    lf : TYPE
        low freq. for bandapss filter.
    hf : TYPE
        high freq. for bandapss filter..
    cpl_dist : TYPE
        Distance between station pair.
    ref_vel : TYPE
        Reference velocity.
    dist_trh : TYPE
        Minimum station separation.
    snr_trh : TYPE
        Signal2noise ratio threshold.
    noise_st : TYPE
        Beginning of the noise window.
    apr_dt_st1 : TYPE
        Apriori estimate of arrival time of station 1.
    apr_dt_st2 : TYPE
        Apriori estimate of arrival time of station 2.
    dt_err : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    sta1 = inv.select(station=station_1)[0][0]
    lat1, lon1 = sta1.latitude, sta1.longitude
    sta2 = inv.select(station=station_2)[0][0]
    lat2, lon2 = sta2.latitude, sta2.longitude
    # Great circle distance in m using WGS84 ellipsoid.
    cpl_dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)[0] 
    min_wl=ref_vel/hf
    if cpl_dist/min_wl < dist_trh:
        print("Station couple does not exceed minimum separation")
        return (np.nan, ['Station couple does not exceed minimum separation',
                         cpl_dist])
    
    params_dir = os.path.join(program_dir, 'params.txt')

    with open(params_dir, 'w') as file:
        file.write("# data_dir, station_1, station_2, dir_name, nsamples, lf, hf, "
        "cpl_dist, ref_vel, dist_trh, snr_trh, noise_st,"
        "apr_dt_st1, apr_dt_st2, dt_err, resp_details \n")

        for val in [data_dir, station_1, station_2, dir_name, nsamples,
                    lf, hf, cpl_dist, ref_vel, dist_trh, snr_trh,
                    noise_st, apr_dt_st1, apr_dt_st2, dt_err, resp_details]:
            file.write(str(val) + '\n')
    os.chdir(program_dir)
    # Uncomment the following line if you made changes to the fortran
    # routines. I don't recommend changing the executable.
    # subprocess.run(["make"], capture_output=True)
    result = subprocess.run(["./BIN/Timing_err_inv"], capture_output=True)
    errors = str(result.stderr).replace("b'","").split('\\n')
    output = str(result.stdout).replace("b'","").split('\\n')
    for a in errors:
        print(a)
    for a in output:
        print(a)
        if 'Result shift:' in a:
            shift = a.split(':')[1]
        if 'Results saved in folder:' in a:
            folder_dir = a.split(':')[1]
    if ('shift' in locals())==False:
            shift = np.nan
            folder_dir = output
    return float(shift), [folder_dir, cpl_dist]

#%%
def plot_signal(tr, maxlag=100, shift=0):
    # define the time vector for the correlation (length of corr = corrwin + 1)
    limit = (len(tr.data) / 2.) * tr.stats.delta
    timevec = np.arange(-limit, limit, tr.stats.delta)
    
    plt.rcParams["figure.figsize"] = [8,4]
    plt.subplot(211)
    plt.plot(timevec, tr.data, 'k', alpha=1, label='before')
    plt.xlim(-maxlag, maxlag)
    plt.legend()
    plt.subplot(212)
    plt.plot([x + shift for x in timevec], tr.data, label='Corrected signal')
    plt.xlim(-maxlag, maxlag)
    plt.xlabel('time [s]')
    plt.legend()
    plt.show()
# plot_signal(tr, maxlag=40)
def trim_correltation_trace(tr, min_t, max_t):
    '''
    Function to trim a cross-correlation trace where the zero is located in 
    the middle of the trace. If you want to trim a trace 5 seconds
    before the zero and 10 seconds after the zero then the fucntion should be
    used as: trim_correltation_trace(tr, min_t=-5, max_t=10)

    Parameters
    ----------
    tr : TYPE
        DESCRIPTION.
    min_t : TYPE
        DESCRIPTION.
    max_t : TYPE
        DESCRIPTION.

    Returns
    -------
    times2 : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.

    '''
    # limit = (len(tr.data) / 2.) * tr.stats.delta
    # timevec = np.arange(-limit, limit, tr.stats.delta)
    
    start = -tr.times()[-1] / 2.
    end = tr.times()[-1] / 2.
    times = np.linspace(start, end, tr.stats.npts)
    
    for i, (amp, t) in enumerate(zip(tr.data, times)):
        if t > min_t:
            low_index = i
            break
    for i, (amp, t) in enumerate(zip(tr.data, times)):
        if t < max_t:
            high_index = i
    tr2 = tr.copy().filter('bandpass', freqmin=0.15, freqmax= 0.3,
                           corners=4, zerophase=True)
    
    times2 = times[low_index: high_index]
    data = tr2.data[low_index: high_index]
    return times2, data

def calculate_apriori(st, t_before, t_after, plot=False):
    '''
    Calculate the aprior dt based on the earliest and latest causal wave.
    It computes cuts the trace between t_before and t_after, then computes 
    the cross-correlation between both traces and returns the shift in time
    between both traces. It assumes that the center of the traces is the zero 
    time.

    Parameters
    ----------
    st : obspy.Stream()
        Obspy stream with cross_correlation from the same station pair.
        It needs to have the attribute average date to know which one is the 
        earliest and which one is the last cross_correlation.
    t_before_zero : float
        can be negative or positive. Take into account that the middle of the
        trace is assumed to be the zero time.
    t_after_zero : float
        can be negative or positive. Take into account that the middle of the
        trace is assumed to be the zero time.
    plot: boolean

    Returns
    -------
    None.

    '''
    avg_dates = [tr.stats.average_date for tr in st]
    earliest_date = avg_dates[0]
    latest_date = avg_dates[0]
    earliest_i = 0
    latest_i = 0
    for i, t in enumerate(avg_dates):
        if t <= earliest_date:
            earliest_date = t
            earliest_i = i
        if t >= latest_date:
            latest_date = t
            latest_i = i
    
    # Now we calculate the cross-correlation between both traces an the shift.
    # This value will be used for the apriori estimates.    
    tr1 = st[earliest_i]
    tr2 = st[latest_i]
    t1, data1 = trim_correltation_trace(tr1, t_before, t_after)
    t2, data2 = trim_correltation_trace(tr2, t_before, t_after)
    
    
    cc = correlate(tr1, tr2, 1000)
    shift, value = xcorr_max(cc)
    
    if plot:
        f, (ax1, ax2) = plt.subplots(2, 1, sharey=True, 
                             figsize=(8,6))
        ax1.plot(t1, data1)
        ax1.plot(t2, data2)
        ax2.plot(t1, data1)
        ax2.plot(t2 + shift/tr1.stats.sampling_rate, data2)
        plt.show()
    # The shift is negative because we are shifting the firs trace,
    # but actually the trace that was shifted is the trace with older 
    # average date.
    time_shift = -shift/tr1.stats.sampling_rate
    
    delta_t = (latest_date - earliest_date)
    shift_rate = time_shift / delta_t
    
    apriori_dt = []
    for t in avg_dates:
        dt = (t - earliest_date) * shift_rate
        apriori_dt.append(dt)
    return apriori_dt

def calculate_apriori_using_a_and_b(a, b, st, reference_time):
    '''
    Calculate the aprior dt based on the earliest and latest causal wave.
    It computes cuts the trace between t_before and t_after, then computes 
    the cross-correlation between both traces and returns the shift in time
    between both traces. It assumes that the center of the traces is the zero 
    time.

    Parameters
    ----------
    st : obspy.Stream()
        Obspy stream with cross_correlation from the same station pair.
        It needs to have the attribute average date to know which one is the 
        earliest and which one is the last cross_correlation.
    t_before_zero : float
        can be negative or positive. Take into account that the middle of the
        trace is assumed to be the zero time.
    t_after_zero : float
        can be negative or positive. Take into account that the middle of the
        trace is assumed to be the zero time.
    plot: boolean

    Returns
    -------
    None.

    '''
    avg_dates = [tr.stats.average_date for tr in st]
    apriori_dt = []
    for t, tr in zip(avg_dates, st):
        # t_N_lps as in Weemstra et al., 2020.
        t_N_lps = (tr.stats.average_date - reference_time)/86400.0
        #ait(lps) + bi
        dt_ins = a * t_N_lps + b
        apriori_dt.append(dt_ins)
    return apriori_dt
    
class Corr_station():
    def __init__(self, station_name, index, corr_with, corr_with_indices,
                 **kwds):
        self._station_name = station_name
        self._index = index
        
        # commented on feb 24.2021
        # self._corr_with = corr_with
        # self._corr_with_indices = corr_with_indices
        if 'corr_directory' in kwds:
            self.corr_directory = (kwds['corr_directory'])
        if 'file_name' in kwds:
            self.file_name = (kwds['file_name'])
            
        # Because we assume that the land stations have no timing errors, 
        # we need to add the type of seismometer. These ones will be assumed to have
        # no timing errors. If you want to change which stations have no errors
        # then you can change this if condition.
        if station_name.startswith('O'):
            self.station_type = 'obs'
        else:
            self.station_type = 'land'

    # function that provides information when used with print statement
    def __repr__(self):
        # Idea: we construct a string with the information and return it:
        # (1) basic information:
        info_string = "\n Corr_station object\n"
        # (2) info on grid cells, spacing, extent:
        info_string += "Station: \t= (%s)\n" % str(self._station_name)
        # info_string += "Index\t= (%d)\n" % (self._index)
        if hasattr(self, 'correlated_with'):
            info_string += "Correlated with stations: " + str(self.correlated_with) + '\n'
            info_string += "Number of correlations: " + str(
                self.number_of_correlations) + '\n'
        if hasattr(self, 'a'):
            info_string += "a value : " + str(self.a)
        if hasattr(self, 'b'):
            info_string += "\nb value : " + str(self.b)
        return info_string
        
class Corr_period():
    
    def __init__(self, station_names, directory, date, file_name,
                 **kwds):
        self.set_station_names(station_names)
        self.set_directory(directory)
        self.set_date(date)
        
        # Key words which are contained in the file name to identify correctly
        # the correlation period.
        self.set_file_name(file_name) 
        
        self.set_stations(station_names)

    # function that provides information when used with print statement
    def __repr__(self):
        # Idea: we construct a string with the information and return it:
        # (1) basic information:
        info_string = "Corr_period object\n"
        # (2) info on grid cells, spacing, extent:
        info_string += "Correlation period with date: \t= (%s)\n" % str(self.date)
        info_string += "Number of stations\t= (%d)\n" % len(self.station_names)
        #info_string += "Files located in directory \t=\n" + str(self.directory)
        
        if all(hasattr(self, attr) for attr in ["matrix_A", "t_app",
                                                "station_pairs", 
                                                "reference_time"]):
            info_string += ("\nMatrix A and the apparent times were " + 
                            "calculated using the following reference time: \n" +
                            str(self.reference_time))
        return info_string
        # return "<{klass} @{id:x} {attrs}>".format(
        #         klass=self.__class__.__name__,
        #         id=id(self) & 0xFFFFFF,
        #         attrs=" ".join("{}={!r}".format(k, v) for k, v in self.__dict__.items()),
        #         )
    
    def __iter__(self):
        #return iter(self.__stations[1:]) #uncomment this if you wanted to 
        # skip the first element.
        return list(self._stations).__iter__()
    
    def __len__(self):
        """
        Return the number of correlation periods in the clock_drift object.
        
        """
        return len(self._stations)

    count = __len__

    def __getitem__(self, index):
        """
        Return the correlation period corresponding to the index.
        
        """
        return self._stations[index]
    
    def set_station_names(self,station_names):
        if isinstance(station_names, list) == False: 
            print("The name of the stations should be in a list of strings!")
            raise AttributeError
        self.station_names = station_names
        
    def set_directory(self,directory):
        if isinstance(directory, str) == False: 
            print("The directory should be a string!")
            raise AttributeError
        self.directory = directory
    
    def set_date(self,date):
        if isinstance(date, obspy.UTCDateTime) == False: 
            print("The date should be in Obspy.UtCDateTime format !")
            raise AttributeError
        self.date = date
        
    def set_file_name(self,file_name):
        if isinstance(file_name, str) == False: 
            print("The file_name should be a string !")
            raise AttributeError
        self.file_name = file_name
    

    def set_stations(self, station_names):
        # Get all the attributes for the stations, like possible correlations
        # with other stations, there indices.
        stations = []
        for i, station1 in enumerate(station_names):
            index = i
            corr_with_indices = []
            corr_with = []
            for j, station2 in enumerate(station_names):#[i+1:], start=i+1):
                if i != j:
                    corr_with_indices.append(j)
                    corr_with.append(station2)
            station = Corr_station(station_name=station1, index=index, 
                                    corr_with=corr_with, 
                                    corr_with_indices=corr_with_indices,
                                    corr_directory=self.directory,
                                    file_name=self.file_name
                                    )
            stations.append(station)
        self._stations = stations
    
    def calculate_shifts(self):
        for station in self._stations:
            shifts = []
            for sta2 in station._corr_with:
                sta1 = station._station_name
                
                st, corr_dirs = read_xcorrelations(sta1, sta2, 
                                                   directory=station.corr_directory,
                                                   )
            
                if any(station.file_name in file_ext for file_ext in corr_dirs):
                    for tr, file_ext in zip(st, corr_dirs):
                        if station.file_name in file_ext:
                            # Calculate the time shift between station 1 and 2.
                            with suppress_stdout():
                                shift, _ = calculate_result_shift( 
                                                               data_dir=file_ext, 
                                                               station_1=sta1, 
                                                               station_2=sta2,
                                                               dir_name=station.file_name,
                                                               nsamples=tr.stats.npts, 
                                                               lf=lf, hf=hf, 
                                                               ref_vel=ref_vel, 
                                                               dist_trh=dist_trh, 
                                                               snr_trh=snr_trh,
                                                               noise_st=noise_st, 
                                                               apr_dt_st1=apr_dt_st1,
                                                               apr_dt_st2=apr_dt_st2, 
                                                               dt_err=dt_err, 
                                                               inv=inv)
                else:
                    shift = np.nan
                
                shifts.append(shift)
                station.shifts = shifts

    def build_matrices(self,
                       lf = lf, # Low freq. for the bandpass filter
                       hf = hf, # High freq. for the bandpass filter 
                       ref_vel = ref_vel, # m/s
                       dist_trh = dist_trh, # Minimum station separation in terms of wavelength
                       snr_trh = snr_trh, # Signal-to-noise ratio threshold
                       noise_st = noise_st, # start of the noise window.
                       dt_err = dt_err, # Sampling interval needs to be multiple of this value.
                       inv=inv,
                       use_apriori = True, #
                       
                       # Reference time with whichwill be taken as the time zero.
                       reference_time = '2014-08-21T00:00:00.000000Z',
                       **kwds
                       # apr_dt_st1 = apr_dt_st1, # A-priori estimate of station 1
                       # apr_dt_st2 = apr_dt_st2, # A-priori estimate of station 2
                       ):
        
        # Matrix and vector with apparent arrival times.
        t_app = [] # vector with the apparent arrival times.
        A = [] # Matrix A shown in the paper Weemstra, Naranjo et al., (2021)
        station_pairs, average_date, cpl_distances, number_of_days = [], [], [], []
        center_date = []
        correlation_ids = []
        n = len(self._stations)*2 # There are two unkowns per station (a and b).
        reference_time = obspy.UTCDateTime(reference_time)
        for i, station1 in enumerate(self._stations):
            for station2 in self._stations[i+1:]: # To avoid repeating stations.
                if station1.station_type == 'obs' or station2.station_type == 'obs':
                    # Only stations that have at least one obs station.
                    sta1 = station1._station_name
                    sta2 = station2._station_name
                    st, corr_dirs = read_xcorrelations(sta1, sta2, 
                                           directory=station1.corr_directory,
                                           )
                    if len(st) == 0:
                        continue
                    if use_apriori:
                        # I included a function to calculate the apriori estimates.
                        apriori_dts_1 = calculate_apriori(st, t_before=-50, t_after=0,
                                                       plot=False)
                        # This method gives a zero shift to the apr of causal
                        # wave.
                        apriori_dts_2 = [0. for i in st]
                    # elif 'apr_dt_st1' not in kwds or 'apr_dt_st2' not in kwds:
                    #     print("If you don't use the apriori_method you need " +
                    #           "give the apriori estimates of station one and" + 
                    #           " station 2. This can be a float.")
                    #     raise AttributeError
                    else:
                        if station1.station_type == 'obs':
                            a_sta1, b_sta1 = float(station1.a), float(station1.b)
                            apriori_dts_1 = calculate_apriori_using_a_and_b(a_sta1,
                                                                            b_sta1,
                                                                            st,
                                                                            reference_time)
                        else:
                            apriori_dts_1 = [0. for i in st]
                        if station2.station_type == 'obs':
                            a_sta2, b_sta2 = float(station2.a), float(station2.b)
                            apriori_dts_2 = calculate_apriori_using_a_and_b(a_sta2,
                                                    b_sta2,
                                                    st,
                                                    reference_time)
                        else:
                            apriori_dts_2 = [0. for i in st]
                    # if sta1 == 'O15' and sta2 == 'KEF':
                    #     print('Aprioris O15:\use apriori method:' + str(use_apriori))
                    #     print(apriori_dts_1)
                    #     print('\n')
                        
                    for tr, f, dt_st1, dt_st2 in zip(st, 
                                                     corr_dirs,
                                                     apriori_dts_1, 
                                                     apriori_dts_2):
                        if self.file_name in f:
                            t_N_lps = (tr.stats.average_date - 
                                       reference_time)/86400.0
                            
                            # To avoid printing useless info on screen.
                            with suppress_stdout(): 
                            
                                shift, [_, cpl_dist] = calculate_result_shift(
                                    data_dir=f, station_1=sta1, 
                                    station_2=sta2, dir_name=station1.file_name,
                                    nsamples=tr.stats.npts, lf=lf, hf=hf, 
                                    ref_vel=ref_vel, dist_trh=dist_trh, 
                                    snr_trh=snr_trh, noise_st=noise_st, 
                                    apr_dt_st1=dt_st1,
                                    apr_dt_st2=dt_st2, 
                                    dt_err=dt_err, inv=inv)

                            if np.isnan(shift) == False:
                                t_app.append(shift)
                                # print(station1._index, station2._index, 
                                #       station1._station_name, 
                                #       station2._station_name,
                                #       self.file_name, shift)
                                a = np.zeros(n)
                                if station1.station_type == 'obs':
                                    a[station1._index*2:station1._index*2+2] = 2
                                    a[station1._index*2] = 2*t_N_lps
                                if station2.station_type == 'obs':
                                    a[station2._index*2:station2._index*2+2] = -2
                                    a[station2._index*2] = -2*t_N_lps
                                A.append(a)
                                station_pairs.append(sta1 + '_' + sta2)
                                #times.append(tr.stats.starttime)
                                cpl_distances.append(cpl_dist)
                                average_date.append(tr.stats.average_date)
                                center_date.append(tr.stats.center_date)
                                number_of_days.append(tr.stats.number_of_days)
                                correlation_ids.append(tr.stats.corr_period_id)
        A = np.asarray(A)
        t_app = np.asarray(t_app)
        self.matrix_A = A
        self.t_app = t_app
        self.station_pairs = station_pairs
        self.average_date = average_date
        self.center_date = center_date
        self.cpl_dist = cpl_distances
        self.number_of_days = number_of_days
        self.correlation_ids = correlation_ids
        self.reference_time = reference_time

#%%
class Clock_drift():
    
    def __init__(self, correlation_periods, 
                 station_names, directory, dates, file_names,
                 **kwds):
        
        self.set_correlation_periods(correlation_periods, 
                 station_names, directory, dates, file_names,
                 **kwds)
        
        self.set_station_names(station_names)
        self.set_directory(directory)
        self.iteration = 0
        info_string = "Clock_drift object initialized with " + str(len(self)) 
        info_string += " correlation periods and " +str(len(self.station_names))
        info_string += " stations. \n"
        info_string += "Number of stations\t= (%d)\n" % len(self.station_names)
        info_string += "Files are located in directory \t=\n" + str(self.directory)
        info_string += "\n \n To calculate the time shifts and build the matrix"
        info_string += "do: clock_drift.build_all_matrix(). \n"
        
    # function that provides information when used with print statement
    def __repr__(self):
        # Idea: we construct a string with the information and return it:
        # (1) basic information:
        info_string = "Clock_drift object\n"
        # (2) info on grid cells, spacing, extent:
        info_string += "There are " + str(len(self)) + " correlation periods"
        info_string += " in Clock_drift object. \n"
        info_string += "Number of stations\t= (%d)\n" % len(self.station_names)
        # info_string += "Files are located in directory \t=\n" + str(self.directory)
        # info_string += "\n \n To calculate the time shifts and build the "
        # info_string += "do: clock_drift.build_all_matrix(). \n"
        if self.iteration > 0:
            info_string += "The clock drift was calculated using " 
            info_string += str(self.iteration) + " iteration(s)."
        return info_string
        # return "{}({!r})".format(self.__class__.__name__, self.__dict__)
    
    def __iter__(self):
        #return iter(self.__stations[1:]) #uncomment this if you wanted to 
        # skip the first element.
        return list(self._corr_periods).__iter__()
    
    def __len__(self):
        """
        Return the number of correlation periods in the clock_drift object.
        
        """
        return len(self._corr_periods)

    count = __len__

    def __getitem__(self, index):
        """
        Return the correlation period corresponding to the index.
        
        """
        return self._corr_periods[index]
    
    def set_correlation_periods(self, correlation_periods, station_names, directory, 
                                dates, file_names):
        corr_periods = []
        for corr_p, date, file_name in zip(correlation_periods, dates, 
                                           file_names):
            t = Corr_period(station_names=station_names, directory=directory, 
                            date=date, file_name = file_name)
            corr_periods.append(t)
        self._corr_periods = corr_periods
        self.dates = dates
        
    def set_station_names(self,station_names):
        if isinstance(station_names, list) == False: 
            print("The name of the stations should be in a list of strings!")
            raise AttributeError
        self.station_names = station_names
    
    def set_directory(self,directory):
        if isinstance(directory, str) == False: 
            print("The directory should be a string!")
            raise AttributeError
        self.directory = directory
    
    # def build_all_matrix(self):
    #     for i, corr_period in enumerate(self._corr_periods):
    #         if i == 0:
    #             with suppress_stdout():
    #                 corr_period.build_matrix()
    #                 corr_period.calculate_shifts()
    #                 corr_period.build_t_app()
    #             A = corr_period._matrixA
    #             t_app = corr_period.t_app
    #         else:
    #             with suppress_stdout():
    #                 corr_period.build_matrix()
    #                 corr_period.calculate_shifts()
    #                 corr_period.build_t_app()
    #             a = corr_period._matrixA
    #             A = np.concatenate((A, a))
    #             t_temp = corr_period.t_app
    #             t_app = np.concatenate((t_app, t_temp))
    #     self.Matrix_A = A
    #     self.t_app = t_app
    
    def build_all_matrix(self, 
                         lf = lf, # Low freq. for the bandpass filter
                         hf = hf, # High freq. for the bandpass filter 
                         ref_vel = ref_vel, # m/s
                         dist_trh = dist_trh, # Minimum station separation in terms of wavelength
                         snr_trh = snr_trh, # Signal-to-noise ratio threshold
                         noise_st = noise_st, # start of the noise window.
                         apr_dt_st1 = apr_dt_st1, # A-priori estimate of station 1
                         apr_dt_st2 = apr_dt_st2, # A-priori estimate of station 1
                         dt_err = dt_err, # Sampling interval needs to be multiple of this value.
                         inv=inv):
        average_date = []
        center_date = []
        station_pairs = []
        cpl_distances = []
        corr_p_ids = []
        number_of_days = []
        A = []
        t_app = []
        for i, corr_period in enumerate(self._corr_periods):
            
            corr_period.build_matrices(lf = lf, hf = hf, ref_vel = ref_vel, 
                                       dist_trh = dist_trh, snr_trh = snr_trh, 
                                       noise_st = noise_st, apr_dt_st1 = apr_dt_st1, 
                                       apr_dt_st2 = apr_dt_st2, dt_err = dt_err, inv=inv)
            A_temp = corr_period.matrix_A
            t_temp = corr_period.t_app
            if len(t_temp) > 0:
                if len(A) == 0:
                    A = A_temp
                    t_app = t_temp
                else:
                    A = np.concatenate((A, A_temp))
                    t_temp = corr_period.t_app
                    t_app = np.concatenate((t_app, t_temp))
                average_date += corr_period.average_date
                center_date += corr_period.center_date
                station_pairs += corr_period.station_pairs
                cpl_distances += corr_period.cpl_dist
                corr_p_ids += corr_period.correlation_ids
                number_of_days += corr_period.number_of_days

        self.t_app = t_app
        self.average_date = average_date
        self.center_date = center_date
        self.station_pairs = station_pairs
        self.corr_period_ids = corr_p_ids
        df = pd.DataFrame(list(zip(station_pairs, average_date, center_date, 
                                   t_app, cpl_distances, 
                                   corr_p_ids, number_of_days)),
                          columns =['Station_pair', 'average_date', 'center_date',
                                    't_app', 'Station_separation[m]',
                                    'Corr_period_id', 'Number_of_days_correlated'])
        self.df = df
        
        column_names = []
        # for i, sta in enumerate(np.repeat(self.station_names, 2)):
        #     if sta.startswith('O'):
        #         if i%2==0:
        #             column_names.append('a*t_{N_lps} ('+sta+')')
        #         else:
        #             column_names.append('b ('+sta +')')
        # print(column_names)
        column_names = ['a*t_{N_lps} ('+sta+')' if i%2==0 else 'b ('+sta +')' for i,
                        sta in enumerate(np.repeat(self.station_names, 2))]
        matrix_A = pd.DataFrame(A, columns = column_names, 
                     index=station_pairs)
        matrix_A = matrix_A.loc[:, (matrix_A != 0).any(axis=0)]
        self.Matrix_A = matrix_A
        
    def solve_eq(self, method='lstsq'):
        
        A_dum = self.Matrix_A.copy()
        Tobs_dum = self.t_app.copy()
        
        if method == 'lstsq':
            x = np.linalg.lstsq(A_dum, Tobs_dum, rcond=None)[0]
            # print ('least squares solution')
            # print(x)
        
        
        elif method == 'lapack':
            nrows = A_dum.shape[0]
            ncols = A_dum.shape[1]
            iwork = 0
            rcond = 0.001
            smlsiz=50     # 25 is probably the approximate value. 50 is on the safe side.
            nlvl = max(0, int(np.log(np.real(min(nrows, ncols ))
                                     / ( smlsiz+1 ) )/np.log(2.)) + 1 )
            lwork=12*ncols + 2*ncols*smlsiz + 8*ncols*nlvl + ncols*1 + (smlsiz+1)**2
            x_lapack, s, rank, info = dgelsd(A_dum, Tobs_dum, lwork, iwork, 
                                      [rcond, 0, 0])
            x = x_lapack[:len(self.Matrix_A.columns)]
            print('lapack solution:')
            print(x)

        else:
            raise NotImplementedError
        
        column_names = [i.replace('*t_{N_lps}','') for i in self.Matrix_A.columns]
        sol = pd.DataFrame(columns = column_names)
        sol.loc['values'] = x
        self.solution = sol
        self.update_stations()
        
    def update_stations(self):
        '''
        Method to update the a and b of each station after using the method
        solve_eq()

        Returns
        -------
        None.

        '''
        
        if hasattr(self, 'solution'):
            for corr_p in self:
                for i, station in enumerate(corr_p._stations):
                    if station.station_type == 'obs':
                        sta_name = station._station_name
                        if 'a ('+sta_name+')' in self.solution.columns:
                            station.a = str(self.solution['a ('+sta_name+')'][0])
                            station.b = str(self.solution['b ('+sta_name+')'][0])
        if hasattr(self, 'Matrix_A'):
            station_pairs = self.station_pairs
            n_correlations_in_clock_drift = 0
            
            for corr_p in self:
                n_correlations_in_corr_p = 0
                
                for station in corr_p._stations:
                    correlated_with = []
                    station_name = station._station_name
                    
                    for pair in station_pairs:
                        if station_name in pair:
                            sta2 = pair.replace('_','').replace(station_name,'')
                            correlated_with.append(sta2)
                    station.number_of_correlations = len(correlated_with)
                    station.correlated_with = list(dict.fromkeys(correlated_with))
                corr_p.number_of_correlations = n_correlations_in_corr_p
            self.number_of_correlations = n_correlations_in_clock_drift
                
    def calculate_app_times(self,
                            use_apriori=True,
                            lf = lf, # Low freq. for the bandpass filter
                            hf = hf, # High freq. for the bandpass filter 
                            ref_vel = ref_vel, # m/s
                            dist_trh = dist_trh, # Minimum station separation in terms of wavelength
                            snr_trh = snr_trh, # Signal-to-noise ratio threshold
                            noise_st = noise_st, # start of the noise window.
                            #apr_dt_st1 = apr_dt_st1, # A-priori estimate of station 1
                            #apr_dt_st2 = apr_dt_st2, # A-priori estimate of station 1
                            dt_err = dt_err, # Sampling interval needs to be multiple of this value.
                            inv=inv,
                            **kwds):
        '''
        Parameters
        ----------
        use_apriori : boolean
            If it's the first iteration then the program will use the default
            method for estimating the apriori values. This method uses the shift
            of the first and last causal wavewithin the different correlation
            periods.
            DESCRIPTION. The default is True.
        program_dir : TYPE, optional
            DESCRIPTION. The default is program_dir.
        lf : TYPE, optional
            DESCRIPTION. The default is lf.
        # Low freq. for the bandpass filter                            hf : TYPE, optional
            DESCRIPTION. The default is hf.
        # High freq. for the bandpass filter                            ref_vel : TYPE, optional
            DESCRIPTION. The default is ref_vel.
        # m/s                            dist_trh : TYPE, optional
            DESCRIPTION. The default is dist_trh.
        # Minimum station separation in terms of wavelength                            snr_trh : TYPE, optional
            DESCRIPTION. The default is snr_trh.
        # Signal-to-noise ratio threshold                            noise_st : TYPE, optional
            DESCRIPTION. The default is noise_st.
        # start of the noise window.                            apr_dt_st1 : TYPE, optional
            DESCRIPTION. The default is apr_dt_st1.
        # A-priori estimate of station 1                            apr_dt_st2 : TYPE, optional
            DESCRIPTION. The default is apr_dt_st2.
        # A-priori estimate of station 1                            dt_err : TYPE, optional
            DESCRIPTION. The default is dt_err.
        # Sampling interval needs to be multiple of this value.                            inv : TYPE, optional
            DESCRIPTION. The default is inv.
        apriori_method : TYPE, optional
            DESCRIPTION. The default is 'standard'.
        **kwds : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        average_date = []
        center_date = []
        station_pairs = []
        cpl_distances = []
        corr_p_ids = []
        number_of_days = []
        A = []
        t_app = []
        for i, corr_period in enumerate(self._corr_periods):
            corr_period.build_matrices(lf = lf, hf = hf, ref_vel = ref_vel, 
                                       dist_trh = dist_trh, snr_trh = snr_trh, 
                                       noise_st = noise_st, 
                                       #apr_dt_st1 = apr_dt_st1, apr_dt_st2 = apr_dt_st2, 
                                       dt_err = dt_err, inv=inv,
                                       use_apriori=use_apriori)
            A_temp = corr_period.matrix_A
            t_temp = corr_period.t_app
            if len(t_temp) > 0:
                if len(A) == 0:
                    A = A_temp
                    t_app = t_temp
                else:
                    A = np.concatenate((A, A_temp))
                    t_temp = corr_period.t_app
                    t_app = np.concatenate((t_app, t_temp))
                average_date += corr_period.average_date
                center_date += corr_period.center_date
                station_pairs += corr_period.station_pairs
                cpl_distances += corr_period.cpl_dist
                corr_p_ids += corr_period.correlation_ids
                number_of_days += corr_period.number_of_days

        self.t_app = t_app
        self.average_date = average_date
        self.center_date = center_date
        self.station_pairs = station_pairs
        self.corr_period_ids = corr_p_ids
        df = pd.DataFrame(list(zip(station_pairs, average_date, center_date, 
                                   t_app, cpl_distances, 
                                   corr_p_ids, number_of_days)),
                          columns =['Station_pair', 'average_date', 'center_date',
                                    't_app', 'Station_separation[m]',
                                    'Corr_period_id', 'Number_of_days_correlated'])
        self.df = df
        
        column_names = []
        # for i, sta in enumerate(np.repeat(self.station_names, 2)):
        #     if sta.startswith('O'):
        #         if i%2==0:
        #             column_names.append('a*t_{N_lps} ('+sta+')')
        #         else:
        #             column_names.append('b ('+sta +')')
        # print(column_names)
        column_names = ['a*t_{N_lps} ('+sta+')' if i%2==0 else 'b ('+sta +')' for i,
                        sta in enumerate(np.repeat(self.station_names, 2))]
        matrix_A = pd.DataFrame(A, columns = column_names, 
                     index=station_pairs)
        matrix_A = matrix_A.loc[:, (matrix_A != 0).any(axis=0)]
        self.Matrix_A = matrix_A
    

    def calculate_clock_drift(self, max_no_iter=2):
        for i in range(0, max_no_iter+1):
            print("Iteration: " + str(self.iteration))
            if self.iteration == 0:
                self.calculate_app_times(use_apriori=True)
                self.solve_eq()
                self.iteration += 1
            else:
                self.calculate_app_times(use_apriori=False)
                self.solve_eq()
                self.iteration += 1

    #%%
    def plot_cross_corr(self, i, maxlag=100):
        '''
        Plot cross-correlation with index i.

        Returns
        -------
        None.

        '''
        
        if hasattr(self, 'df') == False:
            raise AttributeError("First you need to build the DataFrame!")
        
        if isinstance(i, int) == False:
            raise AttributeError("The index must be an integer!")
        
        sta1, sta2 = self.station_pairs[i].split('_')
        corr_period = self.corr_period_ids[i]
        st, corr_dirs = read_xcorrelations(sta1, sta2,
                                           directory=self.directory
                                           )
        for tr_, f in zip(st, corr_dirs):
            if corr_period in f:
                tr = tr_.copy()
        tr = tr.filter('bandpass', freqmin=lf, freqmax=hf, corners=4, 
                       zerophase=True)
        plot_signal(tr, maxlag, shift=-self.t_app[i])

    def plot_solution(self, station, savefig_as=False, title=None):
        last_day = 374
        a_val = float(self.solution['a (' + station + ')'])
        b_val = float(self.solution['b (' + station + ')'])
        dates_from_ref = [(i - self[0].reference_time)/86400.0 for i in self.dates]
        dates_from_ref.insert(0, 0)
        dates_from_ref.append(last_day)
        
        plt.figure(figsize=(8,4))
        plt.plot(dates_from_ref, np.asarray(dates_from_ref)*a_val+b_val, 
                 label='a: ' + str(a_val)[:6] + '\n b: ' + str(b_val)[:6])
        plt.scatter(max(dates_from_ref), max(dates_from_ref)*a_val+b_val)
        plt.annotate(str(max(dates_from_ref)*a_val+b_val)[:5], 
                     (max(dates_from_ref), max(dates_from_ref)*a_val+b_val))
        
        plt.scatter(min(dates_from_ref), min(dates_from_ref)*a_val+b_val)
        plt.annotate(str(min(dates_from_ref)*a_val+b_val)[:5], 
                     (min(dates_from_ref), min(dates_from_ref)*a_val+b_val))
        plt.legend()
        if title==None:
            plt.title('Station ' + station)
        else:
            plt.title(title)
        
        plt.xlabel('Days from ' + str(self[0].reference_time)[:10] +
                   '\n until ' + str(self[0].reference_time + 86400*last_day)[:10])

        plt.ylabel('Clock differences [s]')
        if savefig_as == False:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(savefig_as)
            plt.close()
#%%
# # ################################################################
# directory = '/Users/davidnaranjohernandez/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/'

# station_names = ['O15','KEF']#, 'O22']
# # Setting the Clock  drift object,
# correlation_periods = ["1430827200",
#                         "1426507200",
#                         "1422187200",
#                         "1417867200",
#                         "1413547200"]
# dates = [obspy.UTCDateTime(int(t)) for t in correlation_periods]
# file_names = correlation_periods
# clock_drift = Clock_drift(correlation_periods, station_names, directory, 
#                           dates, file_names)
# clock_drift.build_all_matrix()
# clock_drift.solve_eq(method='lapack')
# print(clock_drift.solution)
# clock_drift.plot_solution(station='O15')

#%%
################################################################
# # Setting the Clock  drift object,
# station_names = ['O15', 'KEF','O22','HOS']
# station_names = [sta.code for sta in inv[0]]
# correlation_periods = ["1430827200",
#                        "1426507200",
#                        "1422187200",
#                        "1417867200",
#                        "1413547200"]
# dates = [obspy.UTCDateTime(int(t)) for t in correlation_periods]
# file_names = correlation_periods
# clock_drift = Clock_drift(correlation_periods, station_names, directory, 
#                           dates, file_names)
# clock_drift.build_all_matrix()
# clock_drift.solve_eq(method='lstsq')

# clock_drift.plot_solution(station='O15')
# clock_drift.plot_solution(station='O22')
