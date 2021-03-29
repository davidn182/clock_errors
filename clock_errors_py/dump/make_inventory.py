#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:20:15 2021
@author: davidnaranjo


Script showing how to make inventory files from station info  in text file.
Modify the station_file in the current directory. The description of each 
of the columns in this file can be seen in the header. The header is skipped
when reading the file.

The only important parameters are the sensor type, code, lat, lon, elev.

If you want to correct the sensors' clock drift 
"""

# Import all the necessary packages.
from obspy.core.inventory import Inventory, Network, Station
import os
# Variables that you can modify.
source = 'dn_IMAGE-project'
input_file = 'station_info' 
# The input file should be in the same directory as this code and have the
# format : header
# PROJECT	SENSORCODE needs_correction(y/n) LATITUDE	LONGITUDE ELEVATION(m) SENSORTYPE
# ...
# We'll first create all the various objects. These strongly follow the
# hierarchy of StationXML files.
inv = Inventory(
                # We'll add networks later.
                networks=[],
                # The source should be the id whoever create the file.
                source=source)

net = Network(
              # This is the network code according to the SEED standard.
              code=".",
              # A list of stations. We'll add one later.
              stations=[],
              description = ("This inventory file is used for correcting the " +
                             "OBS clock errors. We will ignore the different" +
                             " networks present in our data and assume that " +
                             "they all belong to a single network. This " +
                             "facilitates the processing."))

with open(os.path.dirname(__file__) + '/' + input_file) as file:
    rows = file.readlines()[1:]
    for row in rows:
        project = str(row.split()[0])
        sensor_code = str(row.split()[1])
        needs_correction = str(row.split()[2])
        latitude = str(row.split()[3])
        longitude = str(row.split()[4])
        elevation = str(row.split()[5])
        if elevation == '-':
            elevation = 0
        sensor_type = str(row.split()[6])
        
        sta = Station(
                      # This is the station code according to the SEED standard.
                      code=sensor_code,
                      latitude=latitude,
                      longitude=longitude,
                      elevation=elevation)
        sta.sensor_type = sensor_type
        sta.needs_correction = needs_correction
        net.stations.append(sta)
inv.networks.append(net)
print(inv)



# And finally write it to a StationXML file. We also force a validation against
# the StationXML schema to ensure it produces a valid StationXML file.
#
# Note that it is also possible to serialize to any of the other inventory
# output formats ObsPy supports.
inv.plot(projection='local')
#%%
inv.write("inventory.xml", format="stationxml", validate=True)
