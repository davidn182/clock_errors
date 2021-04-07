#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:58:16 2021

@author: davidnaranjo
"""
from contextlib import contextmanager
import obspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.cross_correlation import correlate, xcorr_max
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import sys

# Public functions.
def read_xcorrelations(station1_code, station2_code, path2data_dir):
    '''
    Function to load all the available cross-correlation for a given station
    pair. All the correlation files should be in the same directory and should 
    have the following structure:
    f= station1_station2_averageDateOfCrosscorrelations_numberOfDatesCorrelated

    Parameters
    ----------
    station1_code : string
        name of the station.
    station2_code : string
        name of the station.
    directory : string
        path to data directory

    Returns
    -------
    corr_st : obspy.Stream()
    corr_dirs : list. directory paths of the files

    '''
    if station1_code == station2_code:
        print('No autocorrelations allowed')
        raise
    files = os.listdir(path2data_dir)
    
    # Correlation files in a list to return.
    correlation_stream = obspy.Stream()
    correlation_paths = []
    for file in files:
        if station1_code not in file:
            continue
        if station2_code not in file:
            continue
        if '.sac' not in file:
            continue
        
        correlation_dir = os.path.join(path2data_dir, file)
        correlation_tr = obspy.read(correlation_dir)[0]
        
        # The file header contains the day in the middle of the correlation
        # averaged over the available days
        average_date =  obspy.UTCDateTime(int(file.split('_')[2]))
        number_days = float(file.split('_')[-1].replace('.sac', ''))
    
        correlation_tr.stats.average_date = average_date
        correlation_tr.stats.number_of_days = number_days
        correlation_tr.stats.station_pair = station1_code + '_' + station2_code
        correlation_stream += correlation_tr
        correlation_paths.append(correlation_dir)
    return correlation_stream, correlation_paths

# Public functions.
def read_correlation_file(path2file):
    '''
    Function to load correlation file using the path to the file.

    Parameters
    ----------
    path2file : string
        path to the stream file.

    Returns
    -------
    trace : obspy.Trace()

    '''
    if os.path.isfile(path2file) == False:
        msg = "The file does not exist."
        raise Exception(msg)

    file = os.path.basename(path2file)
    # The file header contains the day in the middle of the correlation
    # averaged over the available days
    splitted_file = file.split("_")
    average_date =  obspy.UTCDateTime(int(splitted_file[2]))
    number_days = float(splitted_file[-1].replace('.sac', ''))
    station1_code, station2_code = str(splitted_file [0]), str(splitted_file [1])
    
    correlation_tr = obspy.read(path2file)[0]
    correlation_tr.stats.average_date = average_date
    correlation_tr.stats.number_of_days = number_days
    correlation_tr.stats.station_pair = station1_code + '_' + station2_code
    
    return correlation_tr

def trim_correlation_trace(tr, min_t, max_t):
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

def calculate_first_apriori_dt(clock_drift_object, correlations, plot=False,
                               **kwargs):
    '''
    Calculates de apriori estimate of given several correlation files of the 
    same station pair, given that the correlation was perform in the same order
    for all the files (meaning station 1 is the same and station 2 is the 
                       same)

    Parameters
    ----------
    clock_drift_object : Clock_drift()
        DESCRIPTION.
    correlations : list
      list of Correlations object. You can use the following function
      to retrieve all the correlations for a given station pair:
      correlations = Clock_drift.get_correlations_of_stationpair(station1_code,
                                                            station2_code)
    if plot is set to tru provide a min_t and t_max to trim the correlation 
    in the times you want to check

    Returns
    -------
    None.

    '''
    freqmin = clock_drift_object.processing_parameters.freqmin
    freqmax = clock_drift_object.processing_parameters.freqmax
    if len(correlations) < 2:
        msg = "There should be at least two correlations to use this method"
        raise Exception(msg)

    sta1 = list(set([correlation.station1_code for correlation in correlations]))
    sta2 = list(set([correlation.station1_code for correlation in correlations]))
    if len(sta1)!= 1 or len(sta2) != 1:
        msg = "The first and second station in the correlations are not the "
        msg += "same for all the correlations."
        raise Exception(msg)

    avg_dates = [correlation.average_date for correlation in correlations]
    
    # Read the correlation of the earliest date
    earliest_date = min(avg_dates)
    earliest_index = avg_dates.index(earliest_date)
    earliest_path2file = correlations[earliest_index].file_path
    earliest_tr = read_correlation_file(path2file=earliest_path2file)
    earliest_tr = earliest_tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                       corners=4, zerophase=True)
    
    # Read the correlation with the latest date.
    latest_date = max(avg_dates)
    latest_index = avg_dates.index(latest_date)
    latest_path2file = correlations[latest_index].file_path
    latest_tr = read_correlation_file(path2file=latest_path2file)
    latest_tr = latest_tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                           corners=4, zerophase=True)
        
    cc = correlate(earliest_tr.data, latest_tr.data, 1000)
    shift, value = xcorr_max(cc)
    time_shift = shift/earliest_tr.stats.sampling_rate
    
    delta_t = (latest_date - earliest_date)
    shift_rate = time_shift / delta_t
    
    for correlation in correlations:
        t = correlation.average_date
        dt = (t - earliest_date) * shift_rate
        if clock_drift_object.get_station(correlation.station1_code).needs_correction ==True:
            correlation.first_apriori_dt1 = dt
            correlation.first_apriori_dt2 = 0
        elif clock_drift_object.get_station(correlation.station2_code).needs_correction ==True:
            correlation.first_apriori_dt1 = 0
            correlation.first_apriori_dt2 = dt
        else:
            raise
        
    if plot:
        min_t = kwargs['min_t']
        max_t = kwargs['max_t']
        t1, data1 = trim_correlation_trace(earliest_tr, min_t, max_t)
        t2, data2 = trim_correlation_trace(latest_tr, min_t, max_t)
        f, (ax1, ax2) = plt.subplots(2, 1, sharey=True, figsize=(8,6))
        ax1.set_title('Before correction ' + earliest_tr.stats.station_pair)
        ax2.set_title('After correction '+ earliest_tr.stats.station_pair)
        ax1.plot(t1, data1, label=earliest_tr.stats.average_date)
        ax1.plot(t2, data2, label=latest_tr.stats.average_date)
        ax2.plot(t1, data1, label=earliest_tr.stats.average_date)
        ax2.plot(t2 + time_shift, data2, label=latest_tr.stats.average_date)
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel('Amplitudes')
        ax1.set_ylabel('Amplitudes')
        plt.tight_layout()
        ax1.legend(loc='best')
        ax2.legend(loc='best')
        plt.show()

def correlations_of_station_exist(station_code, path2data_dir):
    '''
    Function that returns True if there are correlation files for station.
    Remember that the file must contain the station code in the name.
    '''
    for file in os.listdir(path2data_dir):
        if station_code in file:
            return True
    return False

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

# Definitions of the classes ###############################################
class Processing_parameters():
    '''
    '''
    def __init__(self, 
                 freqmin = 0.15, # Low freq. for the bandpass filter
                 freqmax = 0.3, # High freq. for the bandpass filter 
                 ref_vel = 2500, # m/s
                 dist_trh = 2.0, # Minimum station separation in terms of wavelength
                 snr_trh = 10, # Signal-to-noise ratio threshold
                 noise_st = 240, # start of the noise window.
                 dt_err = 0.004, # Sampling interval needs to be multiple of this value.
                 resp_details = False):
        self.freqmin = float(freqmin)
        self.freqmax = float(freqmax)
        self.ref_vel = float(ref_vel)
        self.dist_trh = float(dist_trh)
        self.snr_trh = 10 # Signal-to-noise ratio threshold
        self.noise_st = 240 # start of the noise window.
        self.dt_err = 0.004 # Sampling interval needs to be multiple of this value.
    def __repr__(self):
        msg =  "Processing Parameters Object:"
        msg += "\nMinimum frequency: " + str(self.freqmin)
        msg += "\nMaximum frequency: " + str(self.freqmax)
        msg += "\nReference surface wave velocity: " + str(self.ref_vel)
        msg += "\nMinimum station separation in terms of wavelength: "
        msg += str(self.dist_trh)
        msg += "\nSignal-to-noise ratio threshold: " + str(self.snr_trh)
        msg += "\nStart of the noise window: " + str(self.noise_st)
        return msg

class Station():
    '''
    
    '''
    def __init__(self, code, index, needs_correction, latitude, longitude,
                 elevation, sensor_type, project):
        '''
        Initialize the station object.

        Parameters
        ----------
        code : str
        needs_correction : bool
        latitude : float
        longitude : float
        elevation : float
        sensor_type : str
        project : str
        '''
        # TODO: Check if there is a correlation file with this station.
        self.code = str(code)
        self.sensor_type = str(sensor_type)
        self.project = str(project)
        self.index = index
        
        if needs_correction == 'True':
            self.needs_correction = True
        elif needs_correction == 'False':
            self.needs_correction = False
        else:
            msg = "Error with station " + self.code
            msg += "\nThe value for needs_correction should be True or False."
            raise Exception(msg)
        
        # Only stations that need correction will have a function describing
        # their clock drift. This function is f(t) = at + b.
        # We initialize both variables as lists because we will iteratively
        # calculate the a's and b's until reaching a stable solution.
        if needs_correction: 
            self.a = []
            self.b = []
        try:
            self.latitude = float(latitude)
            self.longitude = float(longitude)
            self.elevation = float(elevation)
        except:
            msg = "Error with station " + self.code
            msg += "\nCheck that the lat, lon or elevation values are real "
            msg += "numbers."
            raise Exception(msg)
        # self.path2corrFile = path2corrFile
    
    def __repr__(self):
        info_string = ("\n Station object\n Code: " + self.code + 
                       "\n Index: " + str(self.index) + 
                       "\n Project: " + self.project + 
                       "\n Sensor type: " + self.sensor_type +
                       "\n Needs correction: " + str(self.needs_correction) +
                       "\n Latitude: " + str(self.latitude) +
                       "\n Longitude: " + str(self.longitude) +
                       "\n Elevation: " + str(self.elevation) + 
                       "\n a values: " + str(self.a) + 
                       "\n b values: " + str(self.b))
        return info_string

class Correlation():
    '''
    '''
    
    def __init__(self, station1_code, station2_code, average_date, number_days,
                 file_path, npts, sampling_rate, length_of_file_s, delta, 
                 cpl_dist):
        self.station1_code = str(station1_code)
        self.station2_code = str(station2_code)
        self.average_date = obspy.UTCDateTime(average_date)
        self.number_days = float(number_days)
        if os.path.exists(file_path) == False:
            msg = "Directory with the correlation file couldn't be found: "
            raise Exception(msg)
        self.file_path = file_path
        self.npts = int(npts)
        self.sampling_rate = float(sampling_rate)
        self.length_of_file_s = float(length_of_file_s)
        self.delta = float(delta)
        self.cpl_dist = cpl_dist
        
        self.t_app = "Not calculated yet."
        self.first_apriori_dt1 = "Not calculated yet."
        self.first_apriori_dt2 = "Not calculated yet."
        self.apriori_station1 = "Not calculated yet."
        self.apriori_station2 = "Not calculated yet."
        
    def __repr__(self):
        info_string = ("\n Correlation object" +
                       "\n Station 1: " + self.station1_code + 
                       "\n Station 2: " + self.station2_code + 
                       "\n Average date of CC: " + str(self.average_date) +
                       "\n Number of days: " + str(self.number_days) +
                       "\n Path: " + self.file_path + 
                       "\n first_apriori_dt1: " + str(self.first_apriori_dt1) +
                       "\n first_apriori_dt2: " + str(self.first_apriori_dt2) +
                       "\n apriori_station1: " + str(self.apriori_station1) +
                       "\n apriori_station2: " + str(self.apriori_station2) +
                       "\n t_app: " + str(self.t_app) + "\n")
        return info_string
    
    def calculate_t_app(self, freqmin, freqmax, ref_vel, dist_trh, snr_trh,
                               noise_st, dt_err, resp_details=False):
        '''
        It uses the last apriori estimate calculated for the station1 and 
        station 2.
    
        Parameters
        ----------
        freqmin : TYPE
            DESCRIPTION.
        freqmax : TYPE
            DESCRIPTION.
        ref_vel : TYPE
            DESCRIPTION.
        dist_trh : TYPE
            DESCRIPTION.
        snr_trh : TYPE
            DESCRIPTION.
        noise_st : TYPE
            DESCRIPTION.
        dt_err : TYPE
            DESCRIPTION.
        resp_details : TYPE, optional
            DESCRIPTION. The default is False.
    
        '''
        
        # Set variables for locating the fortran codes and the parameter file.
        current_path = os.path.abspath(__file__)
        clock_errors_py_dir =(Path(current_path).parents[0])
        program_dir = os.path.join(clock_errors_py_dir,
                                   'recover_timing_errors-master')
        params_dir = os.path.join(program_dir, 'params.txt')
        xcorr_path = self.file_path
        
        results_dir_name = os.path.basename(self.file_path)
        results_dir_name.replace('.sac','')
        results_dir_name.replace(self.station1_code+'-'+self.station2_code,'')
        
        # Station separation. Great circle distance in m using WGS84 ellipsoid.
        cpl_dist = self.cpl_dist
        min_wl=ref_vel/freqmax # Minimum wavelength separation.
        
        if cpl_dist/min_wl < dist_trh:
            self.resp_details = 'Station couple does not exceed minimum separation'
            if type(self.t_app) is not list:
                self.t_app = [np.nan]
            else:
                self.t_app.append(np.nan)
            return(np.nan, ['Station couple does not exceed minimum separation',
                              cpl_dist])
        try:
            # Apriori estimates of both stations.
            apr_dt_st1 = float(self.apriori_station1[-1])
            apr_dt_st2 = float(self.apriori_station2[-1])
        except:
            msg = "No apriori estimate found for station 1 and station 2."
            raise Exception(msg)
        
        with open(params_dir, 'w') as file:
            file.write("# data_dir, station_1, station_2, results_dir_name, "
                       "nsamples, freqmin, freqmax, cpl_dist, ref_vel, dist_trh, snr_trh, "
                       "noise_st, apr_dt_st1, apr_dt_st2, dt_err, resp_details \n")
        
            for val in [xcorr_path, self.station1_code, self.station2_code,
                        results_dir_name, self.npts, freqmin, freqmax, cpl_dist, 
                        ref_vel, dist_trh, snr_trh, noise_st, apr_dt_st1,
                        apr_dt_st2, dt_err, resp_details]:
                file.write(str(val) + '\n')
        os.chdir(program_dir)
        
        result = subprocess.run(["./BIN/Timing_err_inv"], capture_output=True)
        errors = str(result.stderr).replace("b'","").split('\\n')
        output = str(result.stdout).replace("b'","").split('\\n')
        for a in errors:
            print(a)
        for a in output:
            print(a)
            if 'Result shift:' in a:
                shift = float(a.split(':')[1])
            if 'Results saved in folder:' in a:
                folder_dir = a.split(':')[1]
        if ('shift' in locals())==False:
                shift = np.nan
                folder_dir = output
        if resp_details:
            self.resp_details = folder_dir
        if type(self.t_app) is not list:
            self.t_app = [shift]
        else:
            self.t_app.append(shift)
        self.station_separation = cpl_dist
        # TODO: Calculate signal to noise raito.

class Clock_drift():
    '''
    '''
    def __init__(self, station_file, path2data_dir, reference_time, processing_parameters):
        '''
        Parameters
        ----------
        station_file : str file_path
            Path to location of the station file that contains all the
            information of all the stations.
        path2data_dir : str folder_path
            Path to the folder that contains all the cross correlations.
        reference_time : str
            reference time or the date considered the zero time. 
        processing_parameters : Processing_Parameters object
            DESCRIPTION.

        '''
        self.reference_time = obspy.UTCDateTime(reference_time)
        self.set_stations(station_file, path2data_dir)
        self.set_correlations(path2data_dir)
        self.processing_parameters = processing_parameters
    def __repr__(self):
        info_string = "Clock_drift object\n"
        info_string += "There are " + str(len(self.stations)) + " stations"
        info_string += " in Clock_drift object. \n"
        return info_string
        
    def set_correlations(self, path2data_dir):
        if os.path.exists(path2data_dir) == False:
            msg = "Directory with the correlation files couldn't be found: "
            msg += path2data_dir
            raise Exception(msg)
        correlations = []
        for file in sorted(os.listdir(path2data_dir)):
            # The file header contains the day in the middle of the correlation
            # averaged over the available days
            if '.sac' not in file:
                continue
            attributes = file.replace('.sac', '').split('_')
            station1_code = attributes[0]
            station2_code = attributes[1]
            try:
                station1 = self.get_station(station1_code)
                station2 = self.get_station(station2_code)
            except:
                msg = "Couldn´t find the stations of the file: " 
                msg += str(path2data_dir)
                msg += " in the inventory of the stations provided."
                print(msg)
                continue
            
            file_path = os.path.join(path2data_dir, file)
            try:
                tr = read_correlation_file(file_path)
            except:
                msg = "Couldn´t open the file: " 
                msg += str(file_path)
                msg += "\n Remember that it needs to have the format:"
                msg += "station1code_station2code_averageDate_noCorrelatedDays"
                msg += ".sac"
                print(msg)
                continue
            # print(type(tr.stats.average_date))
            # Station separation. Great circle distance in m using WGS84 ellipsoid.
            cpl_dist = gps2dist_azimuth(station1.latitude, station1.longitude,
                                        station2.latitude, station2.longitude)[0]
            average_date =  tr.stats.average_date
            number_days = tr.stats.number_of_days
            npts = tr.stats.npts
            delta = tr.stats.sac.delta
            sampling_rate = tr.stats.sampling_rate
            length_of_file_s = tr.stats.endtime-tr.stats.starttime
            correlation = Correlation(station1_code, station2_code, 
                                      average_date, number_days, file_path,
                                      npts, sampling_rate, length_of_file_s,
                                      delta, cpl_dist)
            t_N_lps = (average_date - self.reference_time)/86400.0
            correlation.t_N_lps = t_N_lps
            correlations.append(correlation)
        self.correlations = correlations
        self.path2data_dir = path2data_dir
        
    def set_stations(self, station_file, path2data_dir):
        if os.path.exists(station_file) == False:
            msg = "Station file couldn't be found: " + station_file
            raise Exception(msg)
        
        stations = []
        with open(station_file) as file:
            rows = file.readlines()[1:]
            index = 0
            for row in rows:
                columns = row.split()
                project = columns[0]
                sensor_code = columns[1]
                needs_correction = columns[2]
                latitude = columns[3]
                longitude = columns[4]
                elevation = columns[5]
                elevation = 0 if elevation == '-' else elevation
                sensor_type = str(row.split()[6])
                if correlations_of_station_exist(sensor_code, path2data_dir):
                    sta = Station(code=sensor_code, index=index,
                                  needs_correction=needs_correction,
                                  latitude=latitude, longitude=longitude,
                                  elevation=elevation, sensor_type=sensor_type,
                                  project=project)
                    stations.append(sta)
                    index += 1
                else:
                    print("No correlation files found for station:" +
                          sensor_code)
        self.stations = stations
        self.station_names = [sta.code for sta in stations]
    
    def get_station(self, station_code):
        for station in self.stations:
            if station.code == station_code:
                return station
        
        msg = "Station not found"
        raise Exception(msg)
    
    def get_correlations_of_stationpair(self, station1_code, station2_code):
        get_correlations = []
        if station1_code == station2_code:
            msg = "You have to choose two different stations."
            raise Exception(msg)
        for correlation in self.correlations:
            if (station1_code != correlation.station1_code and 
                station1_code != correlation.station2_code):
                continue
            if (station2_code != correlation.station1_code and
                station2_code != correlation.station2_code):
                continue
            get_correlations.append(correlation)
        return get_correlations
    
    def calculate_first_aprioridt_4_allcorrelations(self):
        '''
        Function that calculates the first apriori dt for all correlation files. 
        Given the all the stations contained in the clock drift object, the 
        function calculates all the possible station-pair combinations
        and then calculates the apriori estimate for each of the correlations.
    
        Returns
        -------
        None.
    
        '''
        print("Calculating the first apriori estimates for each stationpair")
        for i, station1 in enumerate(self.stations):
            for station2 in self.stations[i+1:]: # To avoid repeating stations.
                correlations = self.get_correlations_of_stationpair(station1.code,
                                                                  station2.code)
                # If there are no corelations for station pair... skip
                if len(correlations) == 0:
                    continue
                
                # If there is only one corelation assume the first estimate as
                # 0 time shift.
                if len(correlations) == 1:
                    correlations[0].first_apriori_dt1 = 0
                    correlations[0].first_apriori_dt2 = 0
                    continue
                # Else calculate the first apriori estimate.
                calculate_first_apriori_dt(self, correlations)

    def calculate_appriori_estimates(self):
        '''
        Calculate the apriori estimates for each of the correlation files.
        Each pariori estimate has a value for the a-causal wave (station 1)
        and a value for the causal wave (station 2).
        
        To calculate the apriori estimates we use the a and b value of both
        stations.
        
        First, we check if the first apriori estimate exists, otherwise it
        gets calculated. 
        Then, we retrieve all the cross-correlations present in the clock
        drift object. 
        After, if there are previously calculated a and b values for the
        stations we calculate the t_N_lps (look at Naranjo et al.,
        2021).
        Otherwise, we use the first_apriori estimates.
        If the station doesn't need correction, the apriori estimate will be 0.
        
        The apriori estimates will be save in a list so that we can check how
        they evolve with each iteration.
        
        Parameters
        ----------

        Returns
        -------
        None.

        '''
        # Check if the first apriori estimate was already calculated. Except:
        #  calculate it.
        try:
            [float(c.first_apriori_dt1) for c in self.correlations]
            print("Calculating the apriori estimates for each stationpair.")
        except:
            self.calculate_first_aprioridt_4_allcorrelations()
            
        for correlation in self.correlations:
            station1 = self.get_station(correlation.station1_code)
            station2 = self.get_station(correlation.station2_code)
            average_date = obspy.UTCDateTime(correlation.average_date)
            t_N_lps = (average_date - self.reference_time)/86400.0
            if station1.needs_correction:
                # If there is an a value or b value for the station
                # use them for calculating the apriori estimates.
                if len(station1.a) > 0 and len(station1.b) > 0:
                    a_val = station1.a[-1] # Take the last calculated a_value
                    b_val = station1.b[-1] # Take the last calculated b_value
                    apriori_station1 = a_val * t_N_lps + b_val
                    correlation.apriori_station1.append(apriori_station1)
                # If there are no a or b values then use the first apriori
                # estimate as the t_app.
                else:
                    correlation.apriori_station1 = [correlation.first_apriori_dt1]
            # If the station doesn't need correction then the a=b=0
            else:
                if len(station1.a) > 0:
                    correlation.apriori_station1.append(0)
                else:
                    correlation.apriori_station1 = [0]
            if station2.needs_correction:
                if len(station2.a) > 0:
                    a_val = station2.a[-1] # Take the last calculated a_value
                    b_val = station2.b[-1] # Take the last calculated b_value
                    apriori_station2 = a_val * t_N_lps + b_val
                    correlation.apriori_station2.append(apriori_station2)
                else:
                    correlation.apriori_station2 = [correlation.first_apriori_dt2]
            else:
                if len(station2.a) > 0:
                    correlation.apriori_station2.append(0)
                else:
                    correlation.apriori_station2 = [0]
                
    def calculate_tapp_4_allcorrelations(self):
        '''
        

        Returns
        -------
        None.

        '''
        freqmin = self.processing_parameters.freqmin
        freqmax = self.processing_parameters.freqmax
        ref_vel = self.processing_parameters.ref_vel
        dist_trh = self.processing_parameters.dist_trh
        snr_trh = self.processing_parameters.snr_trh
        noise_st = self.processing_parameters.noise_st
        dt_err = self.processing_parameters.dt_err
        try:
            for correlation in self.correlations:
                float(correlation.apriori_station1)
                float(correlation.apriori_station1)
        except:
            self.calculate_appriori_estimates()
        
        print("Calculating the t_app for each stationpair.")
        with suppress_stdout():
            for correlation in self.correlations:
                correlation.calculate_t_app(freqmin, freqmax, ref_vel, 
                                            dist_trh, snr_trh, noise_st, dt_err)
                
    def build_matrices(self):
        '''
        Method to build the matrices and vector mentioned in the paper 
        Naranjo et al., (2021)
    
        Returns
        -------
        None.
    
        '''
        # Check first if the apparent times were calculated for all the correlations
        try:
            [float(correlation.t_app[-1]) for correlation in self.correlations]
        except:
            self.calculate_tapp_4_allcorrelations()
        A = [] # Matrix A shown in the paper Weemstra, Naranjo et al., (2021)
        vector_tapp = []
        station1_codes, station2_codes = [], []
        average_date, cpl_distances, number_of_days = [], [], []
        
        # There are two unkowns per station (a and b).
        n = len(self.stations)*2 
        for correlation in self.correlations:
            
            # Get the stations of the cross-correlation file.
            station1 = self.get_station(correlation.station1_code)
            station2 = self.get_station(correlation.station2_code)
            t_N_lps = correlation.t_N_lps
            t_app = correlation.t_app[-1]
            
            # If the apparent time is nan the skip.
            if np.isnan(t_app) == True:
                continue
            
            # Make the row to append matrix A.
            a = np.zeros(n)
            if station1.needs_correction:
                a[station1.index*2:station1.index*2+2] = 2
                a[station1.index*2] = 2*t_N_lps
            if station2.needs_correction:
                a[station2.index*2:station2.index*2+2] = -2
                a[station2.index*2] = -2*t_N_lps
            
            # Store the information of the correlation with shift.
            A.append(a)
            vector_tapp.append(t_app)
            station1_codes.append(correlation.station1_code)
            station2_codes.append(correlation.station2_code)
            average_date.append(correlation.average_date)
            cpl_distances.append(correlation.cpl_dist)
            number_of_days.append(correlation.number_days)
        
        # Make the matrix and vector asarray.
        A = np.asarray(A)
        vector_tapp = np.asarray(vector_tapp)
        
        # We can create a dataframe with all the information.
        df = pd.DataFrame(list(zip(station1_codes, station2_codes, average_date,
                                   vector_tapp, cpl_distances, number_of_days)),
                          columns =['Station1', 'Station2', 'average_date',
                                    't_app[s]', 'Station_separation[m]', 
                                    'Number_of_days_correlated'])
        
        # Then we save matrix A as a dataframe with clear headers.
        column_names = []
        
        # Station pairs is a list with the station pairs of each correlation
        # pair that was used to build the matrix. This will be used in the 
        # dataframe.
        station_pairs = [sta1 + "_" + sta2 for sta1, sta2 in zip(station1_codes,
                                                                 station2_codes)]
        
        # Make the header for the matrix A to make it easier to check.
        for i in range(n):
            for station in self.stations:
                if station.index == i:
                    sta_code = station.code
                    column_names.append("a*t_{N_lps} (" + sta_code + ")")
                    column_names.append("b (" + sta_code + ")")
                    break
        matrix_A = pd.DataFrame(A, columns = column_names, index=station_pairs)
        #  We remove columns that only contain zeros
        matrix_A = matrix_A.loc[:, (matrix_A != 0).any(axis=0)]
        self.matrix_A = matrix_A
        self.df = df
        # Check if there are stations without apparent times.
        for sta in self.stations:
            if sta.needs_correction == True:
                found = False
                for name in matrix_A.columns:
                    if sta.code in name:
                        found = True
                        break
                if found == False:
                    print("No t_app found for Station: " + sta.code)
    def solve_eq(self):
        '''
        
    
        Returns
        -------
        None.
    
        '''
        try:
            A_dum = self.matrix_A.copy()
            Tobs_dum = self.df['t_app[s]'].copy()
        except:
            self.build_matrices()
            A_dum = self.matrix_A.copy()
            Tobs_dum = self.df['t_app[s]'].copy()
        
        print("Inverting the matrix and calculating a and b for each station.")
        x = np.linalg.lstsq(A_dum, Tobs_dum, rcond=None)[0]
        
        column_names = [i.replace('*t_{N_lps}','') for i in self.matrix_A.columns]
        sol = pd.DataFrame(columns = column_names)
        sol.loc['values'] = x
        for value, header in zip(x, column_names):
            if "a" in header:
                station_code = header.replace("a (","").replace(")","")
                station = self.get_station(station_code)
                station.a.append(value)
                continue
            if "b" in header:
                station_code = header.replace("b (","").replace(")","")
                station = self.get_station(station_code)
                station.b.append(value)
        self.solution = sol
    
    ##########################################################################3
    # TODO: Add next step iteration to double check that steps are not repeated.
    ##########################################################################
    
    def plot_before_n_after_first_apriori_estimation(self, 
                                                 station1_code,
                                                 station2_code,
                                                 min_t = -40,
                                                 max_t = 30):
        '''
        Function to generate plot of the cross-correlations before and after 
        applying the correction using the first apriori estimate function.
    
        Parameters
        ----------
        station1_code : TYPE
            DESCRIPTION.
        station2_code : TYPE
            DESCRIPTION.
    
        Returns
        -------
        None.
    
        '''
        correlations = self.get_correlations_of_stationpair(station1_code, 
                                                            station2_code)
        f, (ax1, ax2) = plt.subplots(2, 1, sharey=True, figsize=(8,6))
        f.suptitle('Before and after first apriori estimation')
        for correlation in correlations:
            tr = read_correlation_file(correlation.file_path)
            try:
                first_apriori_dt1 = float(correlation.first_apriori_dt1)
                first_apriori_dt2 = float(correlation.first_apriori_dt2)
                time_shift = first_apriori_dt1 + first_apriori_dt2
            except:
                msg = "You need to calculate the first apriori estimates "
                msg += "before running this function."
                raise Exception(msg)
            t1, data = trim_correlation_trace(tr, min_t, max_t)
            ax1.plot(t1, data, label=str(tr.stats.average_date)[:10])
            ax2.plot(t1 + time_shift, data, label=str(tr.stats.average_date)[:10])
        
        ax1.set_title('Before correction ' + tr.stats.station_pair)
        ax2.set_title('After correction '+ tr.stats.station_pair)
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel('Amplitudes')
        ax1.set_ylabel('Amplitudes')
        ax1.legend(loc=2)
        ax2.legend(loc=2)
        plt.tight_layout()
        plt.show()
        
# if __name__ == "__main__":
#     main()
#%%
# ext = '/Users/localadmin/Dropbox/GitHub/'
# # Parameters
# station_file = ext + "station_info"
# datadir = ext + "data"
# cd = Clock_drift(station_file, datadir, 0)
# cd.calculate_appriori_estimates()

# #%%
# station1_code = 'O01'
# station2_code = 'O20'
# cd.plot_before_n_after_first_apriori_estimation(station1_code, station2_code,
#                                                 min_t = -40, max_t = 30)
# #%%
# freqmin = 0.15 # Low freq. for the bandpass filter
# freqmax = 0.3 # High freq. for the bandpass filter 
# ref_vel = 2500 # m/s
# dist_trh = 2.0 # Minimum station separation in terms of wavelength
# snr_trh = 10 # Signal-to-noise ratio threshold
# noise_st = 240 # start of the noise window.
# dt_err = 0.004 # Sampling interval needs to be multiple of this value.
# resp_details = False
# a = cd.correlations[20]
# #########workflow###########3
# cd.calculate_appriori_estimates()
# cd.calculate_tapp_4_allcorrelations()
# cd.build_matrix
# a.calculate_t_app(freqmin, freqmax, ref_vel, dist_trh, snr_trh, noise_st, dt_err)