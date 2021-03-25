#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:52:13 2021

@author: davidnaranjo
"""
import obspy
import os

data_path='/Users/localadmin/Dropbox/David/PhDTUDelft/data/'
output_dir = '/Users/localadmin/Documents/GitHub/data/'
minutes = 45
for folder in os.listdir(data_path):
    if '.DS_Store' in folder:
        continue
    corrpair_path = os.path.join(data_path, folder)
    
    for corr_file in os.listdir(corrpair_path):
        if 'total' in corr_file:
            continue
        corr_file_path = os.path.join(corrpair_path, corr_file)
        st = obspy.read(corr_file_path, debug_headers=True)
        st.trim(st[0].stats.starttime + 60*minutes,
            st[0].stats.endtime - 60*minutes)
        output_file = os.path.join(output_dir, folder + '_' + corr_file)
    
        st.write(output_file + '.mseed', format='MSEED')