#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 17:27:30 2021

@author: davidnaranjo
"""

import numpy as np
import obspy
import matplotlib.style
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.style.use('bmh')
import sys
sys.path.insert(0, './')
from main import read_xcorrelations, calculate_result_shift, suppress_stdout
from main import Corr_station, Corr_period, Clock_drift
import main_adapted

#%%
# Parameters
directory = '/Users/localadmin/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/'
inv = obspy.read_inventory('/Users/localadmin/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/inventory.xml')
reference_date = obspy.UTCDateTime('2014-08-21T00:00:00.000000Z')

################################################################
# Setting the Clock  drift object,
station_names = ['O16', 'KEF']
correlation_periods = ["1430827200",
                       "1426507200",
                       "1422187200",
                       "1417867200",
                       "1413547200"]
dates = [obspy.UTCDateTime(int(t)) for t in correlation_periods]
file_names = correlation_periods
clock_drift = Clock_drift(correlation_periods, station_names, directory, 
                          dates, file_names)
clock_drift.calculate_app_times()#program_dir='/Users/localadmin/Documents/GitHub/OBS_clock_errors/clock_errors_py/recover_timing_errors-master')

#%%
################################################################
# Setting the Clock  drift object,
station_names = [sta.code for sta in inv[0]]
station_names.remove('ONG')
correlation_periods = ["1430827200",
                       "1426507200",
                       "1422187200",
                       "1417867200",
                       "1413547200"]
dates = [obspy.UTCDateTime(int(t)) for t in correlation_periods]
file_names = correlation_periods
clock_drift2 = Clock_drift(correlation_periods, station_names, directory, 
                          dates, file_names)
#%%
clock_drift2.build_all_matrix(snr_trh=60)
#%%
clock_drift2.solve_eq(method='lstsq')