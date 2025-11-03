#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 23:20:32 2025

@author: decide
"""
%pylab

import pandas as pd
import matplotlib.pylab as plt
import os
import geopandas as gpd
import numpy as np
from dateutil.tz import tzoffset
import matplotlib as mpl

mpl.rcParams.update({"axes.grid" : True})

subfold_name = os.path.join("Tests_Cognin", "Test_tuyau_moutbl")

filename_temp = "time_constant_v_0m_s"
filename_gpx = "trace_gps"
filename_sig = "trace_and_temp"

# Define the sensors
sensors = {"2800028D10000093": "moustiquaire", "28B2968A100000D4": "tuyau"}

# Read temperature
filepath_temp = os.path.join(os.path.abspath(os.path.curdir), "Data", subfold_name, filename_temp + ".csv")
df = pd.read_csv(filepath_temp, header = 0, index_col = 0, parse_dates = True)
# Set timezone to UTC+2 using tzoffset
df.index = df.index.tz_localize(tzoffset(None, 2 * 3600))

# Read gpx
filepath_gpx = os.path.join(os.path.abspath(os.path.curdir), "Data", subfold_name, filename_gpx + ".gpx")
gdf_gpx = gpd.read_file(filepath_gpx, layer = "track_points")

# Calculate distance and time steps and thus speed and filters inconsistent high speed
gdf_gpx = gdf_gpx.reindex(gdf_gpx.time.drop_duplicates().index)
gdf_gpx["distance_step"] = gdf_gpx.to_crs(2154).distance(gdf_gpx.to_crs(2154).shift(-1))
gdf_gpx["time_step"] = gdf_gpx.time.diff()/np.timedelta64(1, 's')
gdf_gpx["speed"] = gdf_gpx["distance_step"].divide(gdf_gpx["time_step"])
gdf_gpx["speed"] = gdf_gpx[gdf_gpx["speed"] <= 10]["speed"]
gdf_speed = gdf_gpx.set_index("time")[["geometry","speed"]]

# Put all data at the same times and rename sensors
union_index = df.index.union(gdf_speed.index)
df = df.reindex(union_index).interpolate(method = "time", limit = 8).reindex(gdf_speed.index)
gdf_speed["speed"] = gdf_speed["speed"].interpolate(method  ="time", limit = 5)\
    .rolling(5).mean()
df.rename(sensors, axis = 1, inplace = True)

# Plot temporal figure
fig, ax = plt.subplots(nrows=2, sharex = True)
df.plot(ax = ax[0], marker = "o")
gdf_speed["speed"].plot(ax = ax[1])
fig.savefig(os.path.join(os.path.abspath(os.path.curdir), "Figures", subfold_name.replace(os.sep, "_") + ".png"))

# Create a GIS file containing both temperature, speed and position
filepath_sig = os.path.join(os.path.abspath(os.path.curdir), "Data", subfold_name, filename_sig + ".fgb")
gdf_all = gpd.GeoDataFrame(pd.concat([df, gdf_speed], axis = 1).reset_index())
gdf_all.rename(sensors, axis = 1, inplace = True)
gdf_all.to_file(filepath_sig)