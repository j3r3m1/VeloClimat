"""
Created on Sun May 11 23:20:32 2025

Averages the difference of air temperature between each cell of a 
grid and one (or several) reference cells of the same grid.
If no grid is given as input, a regular square grid is created. 
To account for shelter + sensor thermal inertia, two methods can be used
for averaging:
    - "lag": data are considered in a cell only after a given duration in the grid 
    (it can be 0 for testing the sensitivity of this method),
    - "moving_avg": a moving average (along time) is used but the value of the window
    is attributed to the beginning of the window


@author: Jérémy Bernard, EDYTEM and Lab-STICC, CNRS
"""

# %pylab

import geopandas as gpd
from sqlalchemy import create_engine
import os
from util import create_grid
from glob import glob
import matplotlib.pyplot as plt

from pathlib import Path
import os
import pandas as pd
import datetime
import numpy as np
from shapely.geometry import Point
from datetime import timedelta, timezone
from dateutil.tz import tzoffset




input_data_path = os.path.abspath(os.path.join(os.path.curdir, "Data", "Inputs"))
output_data_path = os.path.abspath(os.path.join(os.path.curdir, "Data", "Outputs"))


input_data = glob(os.path.join(input_data_path, "*.fgb"))

gdf = {}
for f in input_data:
    name = Path(f).stem
    gdf[name] = gpd.read_file(f)
    
    gdf[name].set_index("timestamp")
    
    # Set timezone to UTC+2 using tzoffset
    if name == "labsticc_sensor0":
        gdf[name].timestamp = pd.DatetimeIndex(gdf[name].timestamp).tz_localize(tzoffset(None, 2 * 3600))
    elif name == "veloclimatmeter0":
        gdf[name].timestamp = pd.DatetimeIndex(gdf[name].timestamp).tz_localize(tzoffset(None, 0 * 3600)).tz_convert('UTC')
    
    gdf[name].drop_duplicates(subset = ["timestamp"], inplace = True)
    gdf[name].set_index("timestamp", inplace = True)
    gdf[name].sort_index(inplace = True)
    
    if name == "veloclimatmeter0":
        gdf[name].index = gdf[name].index + pd.offsets.Second(36)
    
# Put all data at the same times and rename sensors
union_index = gdf["labsticc_sensor0"].index.union(gdf["veloclimatmeter0"].index)
df = pd.concat([gdf["labsticc_sensor0"].temperature.reindex(union_index)\
                .interpolate(method = "time", limit = 8)\
                    .reindex(gdf["veloclimatmeter0"].index).rename("labsticc"),
                gdf["veloclimatmeter0"].temperature.rename("lecea")], axis = 1)

# Plot temporal figure
fig, ax = plt.subplots(nrows=2, sharex = True)
df.plot(ax = ax[0], marker = "o")
gdf["veloclimatmeter0"]["vitesse"].plot(ax = ax[1])
fig.savefig(os.path.join(os.path.abspath(os.path.curdir), "Figures", subfold_name.replace(os.sep, "_") + ".png"))

# Create a GIS file containing both temperature, speed and position
filepath_sig = os.path.join(os.path.abspath(os.path.curdir), "Data", subfold_name, filename_sig + ".fgb")
gdf_all = gpd.GeoDataFrame(pd.concat([df, gdf_speed], axis = 1).reset_index())
gdf_all.rename(sensors, axis = 1, inplace = True)
gdf_all.to_file(filepath_sig)