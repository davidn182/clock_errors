#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:58:16 2021

@author: davidnaranjo
"""
import os

class Station():
    '''
    
    '''
    def __init__(self, code, needs_correction, latitude, longitude,
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
        
        try:
            self.needs_correction = bool(needs_correction)
        except:
            raise Exception("The value for 'needs_correction' should be a" +
                            " boolean")
        try:
            self.latitude = float(latitude)
            self.longitude = float(longitude)
            self.elevation = float(elevation)
        except:
            raise Exception("Check the lat, lon or elevation of station " +
                            self.code + " the values should be numbers.")
        # self.path2corrFile = path2corrFile
    
    def __repr__(self):
        info_string = ("\n Station object\n Code: " + self.code + 
                       "\n Project: " + self.project + 
                       "\n Sensor type: " + self.sensor_type +
                       "\n Needs correction: " + str(self.needs_correction) +
                       "\n Latitude: " + str(self.latitude) +
                       "\n Longitude: " + str(self.longitude) +
                       "\n Elevation: " + str(self.elevation) + '\n')
        return info_string

class Station_couples():
    pass

class Clock_drift():
    '''
    
    '''
    def __init__(self, station_file, datadir):
        self.set_stations(station_file)
        self.set_datadir(datadir)
       
    def __repr__(self):
        info_string = "Clock_drift object\n"
        info_string += "There are " + str(len(self.stations)) + " stations"
        info_string += " in Clock_drift object. \n"
        # if self.iteration > 0:
        #     info_string += "The clock drift was calculated using " 
        #     info_string += str(self.iteration) + " iteration(s)."
        return info_string
        
    def set_datadir(self, datadir):
        if os.path.exists(datadir) == False:
            raise Exception("Directory with the correlation files couldn't " 
                            + "be found")
    def set_stations(self, station_file):
        if os.path.exists(station_file) == False:
            print("Station file couldn't be found")
            raise Exception("Station file couldn't be found")
        
        stations = []
        with open(station_file) as file:
            rows = file.readlines()[1:]
            for row in rows:
                columns = row.split()
                project = columns[0]
                sensor_code = columns[1]
                needs_correction = columns[2]
                latitude = columns[3]
                longitude = columns[4]
                elevation = columns[5]
                if elevation == '-':
                    elevation = 0
                sensor_type = str(row.split()[6])
                
                sta = Station(code=sensor_code, 
                              needs_correction=needs_correction,
                              latitude=latitude, longitude=longitude,
                              elevation=elevation, sensor_type=sensor_type,
                              project=project)
                stations.append(sta)
        self.stations = stations
        self.station_names = [sta.code for sta in stations]
    def calculate_apparent_times(stations):
        pass
#%%
station_file = "/Users/localadmin/Documents/GitHub/OBS_clock_errors/clock_errors_py/inventory/station_info"
datadir = "/Users/localadmin/Dropbox/David/KAUST/Jupyter-notebooks/reykjanes/data/"
cd = Clock_drift(station_file, datadir)
print(len(cd.stations))


# if __name__== "__main__":
#    main()