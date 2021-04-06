#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:44:00 2021

@author: davidnaranjo
"""
from main_adapted import Processing_parameters, Clock_drift
import numpy as np
import pandas as pd
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


cd = Clock_drift(station_file, datadir, reference_time = '2014-08-21T00:00:00.000000Z', processing_parameters=params)
#%%
# cd.calculate_appriori_estimates()
# cd.calculate_tapp_4_allcorrelations()
cd.build_matrices()

cd.solve_eq()
cd.calculate_appriori_estimates()
cd.build_matrices()
cd.solve_eq()
cd.

