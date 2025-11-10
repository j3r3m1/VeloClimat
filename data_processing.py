#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 23:20:32 2025

Main script of the repo to download, pre-process and process the data
acquired during the VeloClimat experiment


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
   
from GlobalVariables import MEAS_NAME, INTERP_REF, VAR, AVG_METH, TIME_LAG, PVAL_MAX,\
    PVAL_MIN, GRID_INFO, GRID_SRID, JDBC_DRIVER_VAR, DB_INFOS, HOST, DBNAME, \
    USERNAME, PASSWORD, SCHEMA_MET, SCHEMA_GEO, LABSTICC_SENSOR_TAB, VELOCLIMATMETER_TAB, \
    LCZ_TAB, SENSOR_NAME_COL, TEMPERATURE_COL, THERMO_NAME_COL, TIMESTAMP_COL,\
    GEOM_COL, LCZ_ID, REFERENCE_NAME, VELOCLIMAMETER, VELOCLIMAP, GIS_EXT,\
    INPUT_DATA_PATH, OUTPUT_DATA_PATH, FIGURES_PATH, LABSTICC_DIR_IN, LABSTICC_DIR_OUT,\
    VELOCLIMAT_RAW_DATA_PATH, LABSTICC_RAW_DATA_PATH, UHI_DATA_PATH, INPUT_RSU_PATH,\
    OUTPUT_RSU_PATH, UHI_LCZ_FIGURE

###############################################################################
############### 1. DECLARE VARIABLES ##########################################
###############################################################################
# -----------------------------------------------------------------------------
# 0. SELECT STEPS THAT NEED TO BE RUN -----------------------------------------
# -----------------------------------------------------------------------------
download = True
cleanAndFormat = True
gather = True
createRef = True
tempDiff = True
investigations = True

###############################################################################
############### 4. SCRIPT BEGINS #################################################
###############################################################################
# If the meas name is a date, download the VeloClimameter data 
if type(MEAS_NAME) == pd._libs.tslibs.timestamps.Timestamp:
    data_type = VELOCLIMAMETER
else:
    data_type = VELOCLIMAP

#------------------------------------------------------------------------------
# 1. DOWNLOAD AND SAVE DATA LOCALLY -------------------------------------------
#------------------------------------------------------------------------------

if download:
    # If not exists, create a folder for the data
    labsticc_dir_in = os.path.join(INPUT_DATA_PATH, MEAS_NAME)
    if not os.path.exists(labsticc_dir_in):
        os.mkdir(labsticc_dir_in)
        
    # Connect to database
    engine = create_engine(f"postgresql://{USERNAME}:{PASSWORD}@{HOST}:5432/{DBNAME}")
    
    # If the meas name is a date, download the VeloClimameter data 
    if data_type == VELOCLIMAMETER:
        sql_velo = f"SELECT * FROM {SCHEMA_MET}.{VELOCLIMATMETER_TAB}"
        gdf_lab = gpd.read_postgis(sql_velo, engine, geom_col=GEOM_COL)
        gdf_lab.to_file(VELOCLIMAT_RAW_DATA_PATH.format(VELOCLIMATMETER_TAB))
    
    # Else download the Lab-STICC data
    elif data_type == VELOCLIMAP:
        if MEAS_NAME == "NULL":
            sql_lab = f"""SELECT * FROM {SCHEMA_MET}.{LABSTICC_SENSOR_TAB} 
                        WHERE THERMO_NAME IS {MEAS_NAME}"""
        else:
            sql_lab = f"""SELECT * FROM {SCHEMA_MET}.{LABSTICC_SENSOR_TAB} 
                          WHERE THERMO_NAME = '{MEAS_NAME}' OR THERMO_NAME = '{MEAS_NAME} - {REFERENCE_NAME}'"""
        gdf_lab = gpd.read_postgis(sql_lab, engine, geom_col=GEOM_COL)
        gdf_lab[[TEMPERATURE_COL, GEOM_COL, TIMESTAMP_COL, SENSOR_NAME_COL]].drop_duplicates()\
            .to_file(LABSTICC_RAW_DATA_PATH.format(LABSTICC_SENSOR_TAB))

    # RSU_LCZ data (need to be change to take the extent of the gdf_lab data
    # as zones to download)
    sql_lcz = f"""SELECT * FROM {SCHEMA_GEO}.{LCZ_TAB}
                  WHERE {LCZ_ID} = '{MEAS_NAME}'"""
    gdf_lcz = gpd.read_postgis(sql_lcz, engine, geom_col=GEOM_COL).to_crs(GRID_SRID)
    gdf_lcz.drop_duplicates().to_file(LABSTICC_RAW_DATA_PATH.format(LCZ_TAB))                    
    

#------------------------------------------------------------------------------
# 1. DOWNLOAD AND SAVE DATA LOCALLY -------------------------------------------
#------------------------------------------------------------------------------
# Clean and format all stations at the same time step to gather all data 
# in one dataFrame, calculate delta from a reference station if needed
if cleanAndFormat:
    # Load the data and remove nan
    gdf_lab = gpd.read_file(LABSTICC_RAW_DATA_PATH.format(LABSTICC_SENSOR_TAB))\
        .dropna(subset = "temperature")
    
    # Identify ref sensors if any
    if data_type == veloclimap:
        ref_sens = gdf_lab[gdf_lab[thermo_name_col] == f'{meas_name} - {reference_name}'][sensor_name_col].unique().tolist()
        uhi_filename = "{0}" + f"_ref_{'&'.join(ref_sens)}"
    else:
        ref_sens = ""

    cleanAndFormat(gdf_lab, ref_sens)


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
                               srid = GRID_SRID)
        gdf_grid.to_file(input_rsu_path.format(rsu))
        
    # Identify measurement points intersecting each grid cell
    # Do it for each station independantly to implement rules for the averaging
    join_fin = {}
    for s in list_mobile_sensors:
        s_obs = gpd.read_file(uhi_data_path.format(uhi_filename.format(s))).to_crs(GRID_SRID)
        
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
    fig.savefig(uhi_lcz_figure.format(meas_name))