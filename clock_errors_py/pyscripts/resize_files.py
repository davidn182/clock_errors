#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:52:13 2021

@author: davidnaranjo
"""
import obspy
import os

<<<<<<< HEAD
data_path='/Users/localadmin/Dropbox/David/PhDTUDelft/data/'
output_dir = '/Users/localadmin/Documents/GitHub/OBS_clock_errors/clock_errors_py/data/'
minutes = 55
for folder in os.listdir(data_path):
    if '.DS_Store' in folder:
        continue
    corrpair_path = os.path.join(data_path, folder)
=======
data_path='/Users/localadmin/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data'


for path in os.listdir(data_path):
    corrpair_path = os.path.join(data_path, path)
>>>>>>> parent of 7d7a291 (Adding the data)
    for corr_file in os.listdir(corrpair_path):
        if 'total' in corr_file:
            continue
        corr_file_path = os.path.join(corrpair_path, corr_file)
<<<<<<< HEAD
        st = obspy.read(corr_file_path, debug_headers=True)
        st.trim(st[0].stats.starttime + 60*minutes,
            st[0].stats.endtime - 60*minutes)
        output_file = os.path.join(output_dir, folder + '_' + corr_file)
    
        st.write(output_file + '.mseed', format='MSEED')
        
    
#%%

st = obspy.read('/Users/localadmin/Downloads/O23_ONG/1417867200_1417882899_0.99.SAC', 
                debug_headers=False)
#%%
sac = obspy.io.sac.SACTrace.read(corr_file_path)
=======
        print(corr_file_path)
    break
        
#%%

st = obspy.read(data_path + 'ARN_O01/1417867200_1418651674_0.82.SAC')
tr = st[0]
print(tr.stats)
minutes = 55
tr.trim(tr.stats.starttime + 60*minutes, tr.stats.endtime - 60*minutes)
print(tr.stats)
tr.write('/Users/localadmin/Downloads/Ttr.mseed',
         format='MSEED')

>>>>>>> parent of 7d7a291 (Adding the data)

# st.write("example.mseed", format="MSEED")  