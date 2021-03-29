#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:58:16 2021

@author: davidnaranjo
"""
import obspy
from obspy.signal.cross_correlation import correlate, xcorr_max
import os
import matplotlib.pyplot as plt
import numpy as np

# Public functions.
def read_xcorrelations(station1, station2, path2data_dir):
    '''
    Function to load all the available cross-correlation for a given station
    pair. All the correlation files should be in the same directory and should 
    have the following structure:
    f= station1_station2_averageDateOfCrosscorrelations_numberOfDatesCorrelated

    Parameters
    ----------
    station1 : string
        name of the station.
    station2 : string
        name of the station.
    directory : string
        path to data directory

    Returns
    -------
    corr_st : obspy.Stream()
    corr_dirs : list. directory paths of the files

    '''
    if station1 == station2:
        print('No autocorrelations allowed')
        raise
    files = os.listdir(path2data_dir)
    
    # Correlation files in a list to return.
    correlation_stream = obspy.Stream()
    correlation_paths = []
    for file in files:
        if station1 not in file:
            continue
        if station2 not in file:
            continue
        if '.mseed' not in file:
            continue
        
        correlation_dir = os.path.join(path2data_dir, file)
        correlation_tr = obspy.read(correlation_dir)[0]
        
        # The file header contains the day in the middle of the correlation
        # averaged over the available days
        average_date =  obspy.UTCDateTime(int(file.split('_')[2]))
        number_days = float(file.split('_')[-1].replace('.mseed', ''))
    
        correlation_tr.stats.average_date = average_date
        correlation_tr.stats.number_of_days = number_days
        correlation_tr.stats.station_pair = station1 + '_' + station2
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
    number_days = float(splitted_file[-1].replace('.mseed', ''))
    station1, station2 = str(splitted_file [0]), str(splitted_file [1])
    
    correlation_tr = obspy.read(path2file)[0]
    correlation_tr.stats.average_date = average_date
    correlation_tr.stats.number_of_days = number_days
    correlation_tr.stats.station_pair = station1 + '_' + station2
    
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
      correlations = Clock_drift.get_correlations_of_stationpair(station1.code,
                                                            station2.code)

    Returns
    -------
    None.

    '''
    if len(correlations) < 2:
        msg = "There should be at least two correlations to use this method"
        raise Exception(msg)

    sta1 = list(set([correlation.station1 for correlation in correlations]))
    sta2 = list(set([correlation.station1 for correlation in correlations]))
    if len(sta1)!= 1 or len(sta2) != 1:
        msg = "The first and second station in the correlations are not the "
        msg += "same for all the correlations."
        raise Exception(msg)

    avg_dates = [correlation.average_date for correlation in correlations]
    
    earliest_date = min(avg_dates)
    earliest_index = avg_dates.index(earliest_date)
    earliest_path2file = correlations[earliest_index].file_path
    earliest_tr = read_correlation_file(path2file=earliest_path2file)
    
    latest_date = max(avg_dates)
    latest_index = avg_dates.index(latest_date)
    latest_path2file = correlations[latest_index].file_path
    latest_tr = read_correlation_file(path2file=latest_path2file)
    
    cc = correlate(earliest_tr.data, latest_tr.data, 1000)
    shift, value = xcorr_max(cc)
    time_shift = shift/earliest_tr.stats.sampling_rate
    
    delta_t = (latest_date - earliest_date)
    shift_rate = time_shift / delta_t
    
    for correlation in correlations:
        t = correlation.average_date
        dt = (t - earliest_date) * shift_rate
        if clock_drift_object.get_station(correlation.station1).needs_correction ==True:
            correlation.first_apriori_dt1 = dt
            correlation.first_apriori_dt2 = 0
        elif clock_drift_object.get_station(correlation.station2).needs_correction ==True:
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
                       "\n Elevation: " + str(self.elevation) + '\n')
        return info_string

class Correlation():
    def __init__(self, station1, station2, average_date, number_days, 
                 file_path):
        self.station1 = str(station1)
        self.station2 = str(station2)
        self.average_date = obspy.UTCDateTime(average_date)
        self.number_days = float(number_days)
        if os.path.exists(file_path) == False:
            msg = "Directory with the correlation file couldn't be found: "
            raise Exception(msg)
        self.file_path = file_path
        self.t_app = "Not calculated yet."
        self.first_apriori_dt1 = "Not calculated yet."
        self.first_apriori_dt2 = "Not calculated yet."
        self.apriori_station1 = "Not calculated yet."
        self.apriori_station2 = "Not calculated yet."
    def __repr__(self):
        info_string = ("\n Correlation object" +
                       "\n Station 1: " + self.station1 + 
                       "\n Station 2: " + self.station2 + 
                       "\n Average date of CC: " + str(self.average_date) +
                       "\n Number of days: " + str(self.number_days) +
                       "\n Path: " + self.file_path + 
                       "\n first_apriori_dt1: " + str(self.first_apriori_dt1) +
                       "\n first_apriori_dt2: " + str(self.first_apriori_dt2) +
                       "\n aprioris_dt1: " + str(self.apriori_station1) +
                       "\n aprioris_dt2: " + str(self.apriori_station2) +
                       "\n t_app: " + str(self.t_app) + "\n")
        return info_string
        
class Clock_drift():
    '''
    
    '''
    def __init__(self, station_file, path2data_dir, reference_time):
        self.set_stations(station_file)
        self.set_correlations(datadir)
        self.reference_time = obspy.UTCDateTime(reference_time)
    def __repr__(self):
        info_string = "Clock_drift object\n"
        info_string += "There are " + str(len(self.stations)) + " stations"
        info_string += " in Clock_drift object. \n"
        return info_string
        
    def set_correlations(self, path2data_dir):
        if os.path.exists(datadir) == False:
            msg = "Directory with the correlation files couldn't be found: "
            msg += datadir
            raise Exception(msg)
        correlations = []
        for file in sorted(os.listdir(path2data_dir)):
            # The file header contains the day in the middle of the correlation
            # averaged over the available days
            if '.mseed' not in file:
                continue
            attributes = file.replace('.mseed', '').split('_')
            station1 = attributes[0]
            station2 = attributes[1]
            average_date =  obspy.UTCDateTime(int(attributes[2]))
            number_days = float(attributes[-1])
            file_path = os.path.join(path2data_dir, file)
            correlation = Correlation(station1, station2, average_date,
                                      number_days, file_path)
            correlations.append(correlation)
        self.correlations = correlations
        self.path2data_dir = path2data_dir
        
    def set_stations(self, station_file):
        if os.path.exists(station_file) == False:
            msg = "Station file couldn't be found: " + station_file
            raise Exception(msg)
        
        stations = []
        with open(station_file) as file:
            rows = file.readlines()[1:]
            for i, row in enumerate(rows):
                columns = row.split()
                project = columns[0]
                sensor_code = columns[1]
                needs_correction = columns[2]
                latitude = columns[3]
                longitude = columns[4]
                elevation = columns[5]
                elevation = 0 if elevation == '-' else elevation
                sensor_type = str(row.split()[6])
                
                sta = Station(code=sensor_code, index=i,
                              needs_correction=needs_correction,
                              latitude=latitude, longitude=longitude,
                              elevation=elevation, sensor_type=sensor_type,
                              project=project)
                stations.append(sta)
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
            if (station1_code != correlation.station1 and 
                station1_code != correlation.station2):
                continue
            if (station2_code != correlation.station1 and
                station2_code != correlation.station2):
                continue
            get_correlations.append(correlation)
        return get_correlations
        
    def calculate_appriori_estimates(self):
        '''
        Calculate the apriori estimates for each of the correlation files.
        Each pariori estimate has a value for the a-causal wave (station 1)
        and a value for the causal wave (station 2).
        
        To calculate the apriori estimates we use the a and b value of both
        stations.

        Parameters
        ----------
        stations : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        try:
            [float(c.first_apriori_dt1) for c in self.correlations]
        except:
            self.calculate_first_aprioridt_4_allcorrelations()
        for correlation in self.correlations:
            station1 = self.get_station(correlation.station1)
            station2 = self.get_station(correlation.station2)
            average_date = obspy.UTCDateTime(correlation.average_date)
            t_N_lps = (average_date - self.reference_time)/86400.0
            if station1.needs_correction:
                if len(station1.a) > 0:
                    a_val = station1.a[-1] # Take the last calculated a_value
                    b_val = station1.b[-1] # Take the last calculated b_value
                    apriori_station1 = a_val * t_N_lps + b_val
                    correlation.apriori_station1.append(apriori_station1)
                else:
                    correlation.apriori_station1 = [correlation.first_apriori_dt1]
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
        for i, station1 in enumerate(self.stations):
            for station2 in self.stations[i+1:]: # To avoid repeating stations.
                correlations = self.get_correlations_of_stationpair(station1.code,
                                                                  station2.code)
                if len(correlations) == 0:
                    continue
                
                if len(correlations) == 1:
                    correlations[0].first_apriori_dt1 = 0
                    correlations[0].first_apriori_dt2 = 0
                    continue
                
                calculate_first_apriori_dt(self, correlations)
    
#%%
ext = '/Users/localadmin/Dropbox/GitHub/'
# Parameters
station_file = ext + "station_info"
datadir = ext + "data"
cd = Clock_drift(station_file, datadir, 0)

cd.calculate_appriori_estimates()

print(len(cd.get_correlations_of_stationpair('KEF', 'RAR')))