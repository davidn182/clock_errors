#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 21:03:52 2021

@author: davidnaranjo
"""

from main_adapted import Processing_parameters, Clock_drift
import time
ext = '/Users/localadmin/Dropbox/GitHub/'
# Parameters
station_file = ext + "station_info"
datadir = ext + "data"


params = Processing_parameters(
                 freqmin = 0.15, # Low freq. for the bandpass filter
                 freqmax = 0.3, # High freq. for the bandpass filter 
                 ref_vel = 2500, # m/s
                 dist_trh = 2.0, # Minimum station separation in terms of wavelength
                 snr_trh = 10, # Signal-to-noise ratio threshold
                 noise_st = 240, # start of the noise window.
                 dt_err = 0.004, # Sampling interval needs to be multiple of this value.
                 resp_details = False)
#%%

start_time = time.time()
cd = Clock_drift(station_file, datadir, reference_time = '2014-08-21T00:00:00.000000Z', processing_parameters=params)
print("---Initializing the Clock_drift object takes %s seconds ---" % (time.time() - start_time))

#%%
start_time = time.time()
cd.calculate_appriori_estimates()
print("---Calculating apriori times takes --- %s seconds ---" % (time.time() - start_time))
#%%
start_time = time.time()
cd.calculate_tapp_4_allcorrelations()
print("---Calculating apparent times takes --- %s seconds ---" % (time.time() - start_time))

#%%
start_time = time.time()
cd.build_matrices()
print("---Building the matrix takes --- %s seconds ---" % (time.time() - start_time))

#%%
start_time = time.time()
cd.solve_eq()
print("---Inverting the matrix takes --- %s seconds ---" % (time.time() - start_time))
#%%
# ###################3 Some plotting functions ###############################
station1_code = 'O01'
station2_code = 'O20'
cd.plot_before_n_after_first_apriori_estimation(station1_code, station2_code,
                                                min_t = -40, max_t = 30)

#%%
start_time = time.time()
cd = Clock_drift(station_file, datadir, reference_time = '2014-08-21T00:00:00.000000Z', processing_parameters=params)
cd.calculate_appriori_estimates()
cd.calculate_tapp_4_allcorrelations()
cd.build_matrices()
cd.solve_eq()
print("---One iteration takes --- %s seconds ---" % (time.time() - start_time))

for i in range(10):
    start_time = time.time()
    cd.calculate_appriori_estimates()
    cd.calculate_tapp_4_allcorrelations()
    cd.build_matrices()
    cd.solve_eq()
    print("---One iteration takes --- %s seconds ---" % (time.time() - start_time))
#%%
# ########## TO save yur progress ########################################
import pickle
dump_dir = ext + "clock_errors/clock_errors_py/dump/"
dump_file = dump_dir + 'clock_drift.obj'
file_pi = open(dump_file, 'wb') 
pickle.dump(cd, file_pi)

with open(dump_file, 'rb') as f:
        cd2 = pickle.load(f)
        
"/Users/localadmin/Dropbox/GitHub/clock_errors/clock_errors_py/jupyter-tutorials/data_test/",
"/Users/localadmin/Dropbox/GitHub/clock_errors/clock_errors_py/jupyter-tutorials/station_info"
# Parameters for locating the files where the correlation files and station 
# information is contained.