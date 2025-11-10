#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
from pathlib import Path

import os
import pandas as pd
import datetime
import numpy as np
from shapely.geometry import Point
from datetime import timedelta, timezone
import matplotlib.pyplot as plt
# import pandas as pd
# import matplotlib.pylab as plt
# import os
# import geopandas as gpd
# import numpy as np
# from dateutil.tz import tzoffset
# import matplotlib as mpl

# mpl.rcParams.update({"axes.grid" : True})

# File paths











    


###############################################################################
############### 1. DECLARE VARIABLES ##########################################
###############################################################################
# Name of the dataset (thermoparty) of interest
# Possible solutions
#   - "NULL"
#   - "Saint-Jean La Poterie"
#   - "Laval"
#   - "Alençon"
#   - "Amiens"
#   - "Lille"
thermo_name = "Gand"

# Select steps that need to be run
step1 = True        # Download the data
step2 = True        # Interpolates the data along time and calculate UHI from ref
step3 = True        # Calculates average values within grid cells

# Reference sensor(s) - it is a list since it can be the mean of several sensors
ref_sens = []

# Whether or not we interpolate linearly the ref if no data
interp_ref = True

# Name of the variable of interest ("deltaT", "T")
var = "deltaT"

# Start and end time of the measurement campaign in UTC
campaign_start = datetime.datetime(2025,7,3,18,35,0, tzinfo = datetime.timezone.utc)
campaign_end = datetime.datetime(2023,6,12,21,33,39, tzinfo = datetime.timezone.utc)

# Path of the file to use as grid / else size of the gri used to average the data (m)
grid_info = 0
grid_srid = 2154

# Averaging method chosen ("lag" or "moving_avg")
avg_meth = "moving_avg"

# Time spent in a grid cell before taking into account the data / 
# windows size for moving average (s)
time_lag = 10

# Extension for GIS files
gis_ext = ".fgb"

# Percentiles to use for min and max values per grid cell
pval_max = 0.9
pval_min = 0.1

###############################################################################
############### 2. ALREADY SET VARIABLES ######################################
###############################################################################
# Database config
jdbc_driver_jar = "/home/decide/.m2/repository/org/postgresql/postgresql/42.6.0/postgresql-42.6.0.jar"
db_infos = pd.read_csv("/home/decide/Code/DB_informations.csv",
                       header = 0, 
                       index_col = 0)
host = db_infos.loc["veloclimat", "host"]
dbname = db_infos.loc["veloclimat", "dbname"]
username = db_infos.loc["veloclimat", "username"]
password = db_infos.loc["veloclimat", "password"]


# Schema and table names in the database
schema_met = "veloclimat"
schema_geo = "Mapuce"
labsticc_sensor_tab = "labsticc_sensor"
veloclimatmeter_tab = "veloclimatmeter"
lcz_tab = "grid_indicators_100"
lcz_tab = "rsu_lcz"

# Column names
sensor_name_col = "sensor_name"
temperature_col = "temperature"
thermo_name_col = "thermo_name"
timestamp_col = "timestamp"
geom_col = "the_geom"
lcz_id = "ID_ZONE"

# Reference
reference_name = "Reference"

###############################################################################
############### 3. DECLARE PATHS ##############################################
###############################################################################
input_data_path = os.path.abspath(os.path.join(os.path.curdir, "Data", "Inputs"))
output_data_path = os.path.abspath(os.path.join(os.path.curdir, "Data", "Outputs"))
figures_path = os.path.abspath(os.path.join(os.path.curdir, "Figures"))

labsticc_dir_in = os.path.join(input_data_path, thermo_name)
labsticc_dir_out = os.path.join(output_data_path, thermo_name)

veloclimat_raw_data_path = os.path.join(input_data_path, "{0}" + f"{gis_ext}")
labsticc_raw_data_path = os.path.join(labsticc_dir_in, "{0}" + f"{gis_ext}")
uhi_data_path = os.path.join(labsticc_dir_out, "{0}" + f"{gis_ext}")
input_rsu_path = os.path.join(labsticc_dir_in, "{0}" + f"{gis_ext}")
output_rsu_path = os.path.join(labsticc_dir_out, "{0}" + f"_fin_{avg_meth}{gis_ext}")

uhi_filename = "{0}" + f"_ref_{'&'.join(ref_sens)}"
uhi_lcz_figure = os.path.join(figures_path, "{0}_uhi_lcz.png")

###############################################################################
############### 4. SCRIPT BEGINS #################################################
###############################################################################
# If not exists, create a folder for the data and figures
if not os.path.exists(labsticc_dir_in):
    os.mkdir(labsticc_dir_in)
if not os.path.exists(labsticc_dir_out):
    os.mkdir(labsticc_dir_out)
if not os.path.exists(figures_path):
    os.mkdir(figures_path)

#------------------------------------------------------------------
# 1. SAVE DATA LOCALLY
#------------------------------------------------------------------
if step1:
    # Connect to database
    engine = create_engine(f"postgresql://{username}:{password}@{host}:5432/{dbname}")
    
    # # Read data via SQL query, download locally and save locally
    # # VeloClimameter data
    # sql_velo = f"SELECT * FROM {schema_met}.{veloclimatmeter_tab}"
    # gdf_lab = gpd.read_postgis(sql_velo, engine, geom_col=geom_col)
    # gdf_lab.to_file(veloclimat_raw_data_path.format(veloclimatmeter_tab))
    
    # # Lab-STICC data
    # if thermo_name == "NULL":
    #     sql_lab = f"SELECT * FROM {schema_met}.{labsticc_sensor_tab} WHERE THERMO_NAME IS {thermo_name}"
    # else:
    #     sql_lab = f"""SELECT * FROM {schema_met}.{labsticc_sensor_tab} 
    #                     WHERE THERMO_NAME = '{thermo_name}' OR THERMO_NAME = '{thermo_name} - {reference_name}'"""
    # gdf_lab = gpd.read_postgis(sql_lab, engine, geom_col=geom_col)
    # if not ref_sens:
    #     ref_sens = gdf_lab[gdf_lab[thermo_name_col] == f'{thermo_name} - {reference_name}'][sensor_name_col].unique().tolist()
    # gdf_lab[[temperature_col, geom_col, timestamp_col, sensor_name_col]].drop_duplicates()\
    #     .to_file(labsticc_raw_data_path.format(labsticc_sensor_tab))

    # RSU_LCZ data
    sql_lcz = f"""SELECT * FROM {schema_geo}.{lcz_tab}
                    WHERE {lcz_id} = '{thermo_name}'"""
    gdf_lcz = gpd.read_postgis(sql_lcz, engine, geom_col=geom_col).to_crs(grid_srid)
    gdf_lcz.drop_duplicates().to_file(labsticc_raw_data_path.format(lcz_tab))                    
    
#------------------------------------------------------------------
# 2. INTERPOLATES ALONG TIME AND CALCULATE UHI FROM REF
#------------------------------------------------------------------
# Load the data and remove nan
gdf_lab = gpd.read_file(labsticc_raw_data_path.format(labsticc_sensor_tab)).dropna(subset = "temperature")

#gdf_lab = gdf_lab[pd.DatetimeIndex(gdf_lab["timestamp"]).hour >= campaign_start.hour]
# Save one file for each sensor
list_sensors = gdf_lab[sensor_name_col].unique()
list_mobile_sensors = [x for x in list_sensors if x not in ref_sens]
if step2:
    gdf_lab_dic = {}
    for s in list_sensors:
        gdf_lab_dic[s] = gdf_lab[gdf_lab[sensor_name_col] == s]
        gdf_lab_dic[s].rename({temperature_col: f"{temperature_col}{s}"}, inplace = True, axis = 1)
        gdf_lab_dic[s][[timestamp_col, "geometry", f"{temperature_col}{s}"]]
        gdf_lab_dic[s].set_index(timestamp_col, inplace = True)
        # Drop duplicated time steps (keep only first)
        gdf_lab_dic[s] = gdf_lab_dic[s][~gdf_lab_dic[s].index.duplicated(keep='first')].sort_index()
    
    # Temporary line to create training dataset
    # gdf_lab_dic[s] = gdf_lab_dic[s].shift(freq = -pd.offsets.Day(5) + pd.offsets.Minute(13))
    
    # Read meteo data
    df_all_meteo = pd.concat([gdf_lab_dic[s][[f"{temperature_col}{s}"]] 
                              for s in list_sensors], 
                              axis = 1).sort_index()
    
    # Reindex reference air temperature to all indexes, interpolate linearly with time
    # and average the reference stations
    for s in list_sensors:
        if interp_ref:
            # Reference meteorological condition
            df_meteo_ref = df_all_meteo[[f"{temperature_col}{s_i}" for s_i in ref_sens]]\
                .interpolate(method = "time").mean(axis = 1)
        else:
            # Reference meteorological condition
            df_meteo_ref = df_all_meteo[[f"{temperature_col}{s_i}" for s_i in ref_sens]]\
                .interpolate(method = "time", limit = len(list_sensors) + 1)\
                    .mean(axis = 1)

    
    gdf_uhi_dic = {s: gdf_lab_dic[s].copy(deep = True) for s in gdf_lab_dic.keys()}
    for s in list_mobile_sensors:
        if var == "deltaT":
            gdf_uhi_dic[s][f"{temperature_col}{s}"] = \
                gdf_lab_dic[s][f"{temperature_col}{s}"].subtract(df_meteo_ref)\
                    .dropna()
        elif var == "T":
            gdf_uhi_dic[s][f"{temperature_col}{s}"] = \
                gdf_lab_dic[s][f"{temperature_col}{s}"].copy(deep=True)
        gdf_uhi_dic[s].rename({f"{temperature_col}{s}" : var}, 
                              axis = 1).dropna().to_file(uhi_data_path.format(uhi_filename.format(s)))
    
#------------------------------------------------------------------
# 2. AVERAGE THE VARIABLES WITHIN GRID CELL
#------------------------------------------------------------------
if step2:
    print("Step 2 started")
    # If a grid file is given as input
    if grid_info == 0:
        lcz_file = labsticc_raw_data_path.format(lcz_tab)
        if os.path.exists(lcz_file):
            rsu = Path(lcz_file).stem
            gdf_grid = gpd.read_file(lcz_file)
            # Set the index as ID
            gdf_grid["ID"] = gdf_grid.index
        
    # Else create a grid file using square grid cells
    else:
        rsu = "grid"
        xmin, ymin, xmax, ymax = \
            pd.concat({s: gpd.read_file(uhi_data_path.format(uhi_filename.format(s)))
                             for s in list_mobile_sensors}, 
                            ignore_index = True).sort_values(axis = 0,
                                                             by = timestamp_col).total_bounds
        s_bbox_corners = gpd.GeoSeries([Point(xmin, ymin), Point(xmax, ymax)],
                                     crs = 4326)
        gdf_grid = create_grid(s_bbox_corners = s_bbox_corners,
                               cell_size = grid_info,
                               srid = grid_srid)
        gdf_grid.to_file(input_rsu_path.format(rsu))
        
    # Identify measurement points intersecting each grid cell
    # Do it for each station independantly to implement rules for the averaging
    join_fin = {}
    for s in list_mobile_sensors:
        s_obs = gpd.read_file(uhi_data_path.format(uhi_filename.format(s))).to_crs(grid_srid)
        
        if avg_meth == "lag" or avg_meth == "moving_avg":
            # Identify the intersection between station location and grid cell
            # and sort in ascending datetime
            join = gpd.sjoin(left_df=s_obs[[timestamp_col,"geometry", var]],
                             right_df=gdf_grid[["ID","geometry"]], 
                             how="left", 
                             predicate="intersects").sort_values(timestamp_col)        
            if avg_meth == "lag":
                # Add a column with the time where the station get into a new grid cell
                join["in_time"] = join[timestamp_col][join["ID"].diff() != 0].reindex(join.index)
                # Calculate the time spent in the grid since the entrance
                join["in_time"] = join[timestamp_col].subtract(pd.Series(join.set_index(timestamp_col)["in_time"].ffill().values,
                                                                         index = join.index).dt.tz_localize(timezone(timedelta(hours = 0))))
                # Remove points being in a cell since less than 'time_lag'
                join_fin[s] = join[join["in_time"] >= pd.Timedelta(pd.offsets.Second(time_lag))]\
                    .drop(["in_time"], axis = 1)
                
            elif avg_meth == "moving_avg":
                # First move the air temperature further in the future
                join = join.set_index(timestamp_col)
                join_T = join[var].shift(periods = 1,
                                         freq = pd.offsets.Second(-time_lag))
                # Then calculate a rolling average, else there is no more geometries for time some values
                join[var] = join_T.reindex(join_T.index.union(join.index))\
                    .rolling(window = pd.offsets.Second(time_lag), 
                             min_periods = int(time_lag / 4)).mean()\
                        .reindex(join.index)
                join_fin[s] = join.copy(deep = True)
        else:
            print(f"The averaging method '{avg_meth}' does not exists. Please select a valid one.")
        
    # Gather all measurements in a single dataframe
    grid_fin = pd.concat(join_fin, ignore_index = True)
    
    # Aggregate values by grid cell using some statistics
    gdf_grid_fin = gdf_grid.set_index("ID")
    gdf_grid_fin[f"{var}_mean"] = grid_fin.groupby("ID")[var].mean()
    gdf_grid_fin[f"{var}_q{pval_max}"] = grid_fin.groupby("ID")[var].quantile(pval_max)
    gdf_grid_fin[f"{var}_q{pval_min}"] = grid_fin.groupby("ID")[var].quantile(pval_min)
    gdf_grid_fin[f"{var}_median"] = grid_fin.groupby("ID")[var].median()
    gdf_grid_fin[f"{var}_sqrt"] = grid_fin.groupby("ID")[var].std()
    
    # Save the resulting file
    gdf_grid_fin.to_file(output_rsu_path.format(rsu))
    
    # Plot the results by LCZ type
    fig, ax = plt.subplots(figsize = (10, 7))
    gdf_grid_groupbylcz = gdf_grid_fin[["lcz_primary", "deltaT_mean"]].groupby("lcz_primary")
    gdf_grid_groupbylcz.boxplot(column = "deltaT_mean", subplots = False,
                                rot = 35, ax = ax)
    ax.set_xticklabels(pd.Series(gdf_grid_groupbylcz.groups.keys(), dtype = int))
    ax.set_xlabel("LCZ type")
    ax.set_ylabel(u"Différence de température / référence (°C)")
    fig.savefig(uhi_lcz_figure.format(thermo_name))