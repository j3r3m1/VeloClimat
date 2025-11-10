#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 00:02:59 2025

@author: decide
"""

from digitemp.master import UART_Adapter
from digitemp.device import DS18B20


import os
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from GlobalVariables import MEAS_NAME, INTERP_REF, VAR, AVG_METH, TIME_LAG, PVAL_MAX,\
    PVAL_MIN, GRID_INFO, GRID_SRID, JDBC_DRIVER_VAR, DB_INFOS, HOST, DBNAME, \
    USERNAME, PASSWORD, SCHEMA_MET, SCHEMA_GEO,  LABSTICC_SENSOR_TAB, VELOCLIMATMETER_TAB, \
    LCZ_TAB, SENSOR_NAME_COL, TEMPERATURE_COL, THERMO_NAME_COL, TIMESTAMP_COL,\
    GEOM_COL, LCZ_ID, REFERENCE_NAME, VELOCLIMAMETER, VELOCLIMAP, GIS_EXT,\
    INPUT_DATA_PATH, OUTPUT_DATA_PATH, FIGURES_PATH, LABSTICC_DIR_IN, LABSTICC_DIR_OUT,\
    VELOCLIMAT_RAW_DATA_PATH, LABSTICC_RAW_DATA_PATH, UHI_DATA_PATH, INPUT_RSU_PATH,\
    OUTPUT_RSU_PATH, UHI_LCZ_FIGURE

def read_temp(adress):
    try:
        sensor = DS18B20(UART_Adapter(adress))
        # bus = sensor.bus.name
        temp = sensor.get_temperature()
        rom_id = sensor.rom_code.hex().upper()
        sensor = None
    except:
        temp = None
        rom_id = None
    
    return rom_id, temp


def create_grid(s_bbox_corners, cell_size, srid):
    # Convert to local CRS
    bbox_corner_srid = s_bbox_corners.to_crs(srid)
    
    # Create a small extent to have cells starting a bit further than at 
    # the extreme location of the measurement
    x_ll = bbox_corner_srid.x[0] - cell_size / 2
    y_ll = bbox_corner_srid.y[0] - cell_size / 2
    x_ur = bbox_corner_srid.x[1] + cell_size / 2
    y_ur = bbox_corner_srid.y[1] + cell_size / 2
    
    # Calculate the number of cells along x and y axis
    nx = int(np.trunc((x_ur - x_ll) / cell_size) + 1)
    ny = int(np.trunc((y_ur - y_ll) / cell_size) + 1)
    
    gdf = gpd.GeoDataFrame({"geometry" : [Polygon([(x_ll + i * cell_size, y_ll + j * cell_size),
                                         (x_ll + (i + 1) * cell_size, y_ll + j * cell_size), 
                                         (x_ll + (i + 1) * cell_size, y_ll + (j + 1) * cell_size),
                                         (x_ll + i * cell_size, y_ll + (j + 1) * cell_size),
                                         (x_ll + i * cell_size, y_ll + j * cell_size)])
                            for i in range(nx) for j in range(ny)]},
                           crs = srid)
    gdf["ID"] = [i + 1 for i in range(gdf.index.size)]
    
    return gdf

def cleanAndFormat(gdf_lab, ref_sens, uhi_filename):
    # Save one file for each sensor
    list_sensors = gdf_lab[SENSOR_NAME_COL].unique()
    list_mobile_sensors = [x for x in list_sensors if x not in ref_sens]
    # If not exists, create a folder for the data
    if not os.path.exists(LABSTICC_DIR_OUT):
        os.mkdir(LABSTICC_DIR_OUT)
        
    gdf_lab_dic = {}
    for s in list_sensors:
        gdf_lab_dic[s] = gdf_lab[gdf_lab[SENSOR_NAME_COL] == s]
        gdf_lab_dic[s].rename({TEMPERATURE_COL: f"{TEMPERATURE_COL}{s}"}, inplace = True, axis = 1)
        gdf_lab_dic[s][[TIMESTAMP_COL, "geometry", f"{TEMPERATURE_COL}{s}"]]
        gdf_lab_dic[s].set_index(TIMESTAMP_COL, inplace = True)
        
        # Drop duplicated time steps (keep only first)
        gdf_lab_dic[s] = gdf_lab_dic[s][~gdf_lab_dic[s].index.duplicated(keep='first')].sort_index()
    

    # Read meteo data
    df_all_meteo = pd.concat([gdf_lab_dic[s][[f"{TEMPERATURE_COL}{s}"]] 
                              for s in list_sensors], 
                              axis = 1).sort_index()
    
    # Reindex reference air temperature to all indexes, interpolate linearly with time
    # and average the reference stations
    for s in list_sensors:
        if INTERP_REF:
            # Reference meteorological condition
            df_meteo_ref = df_all_meteo[[f"{TEMPERATURE_COL}{s_i}" for s_i in ref_sens]]\
                .interpolate(method = "time").mean(axis = 1)
        else:
            # Reference meteorological condition
            df_meteo_ref = df_all_meteo[[f"{TEMPERATURE_COL}{s_i}" for s_i in ref_sens]]\
                .interpolate(method = "time", limit = len(list_sensors) + 1)\
                    .mean(axis = 1)

    # Calculate the temperature difference from the reference sensor if needed
    gdf_uhi_dic = {s: gdf_lab_dic[s].copy(deep = True) for s in gdf_lab_dic.keys()}
    for s in list_mobile_sensors:
        if VAR == "deltaT":
            gdf_uhi_dic[s][f"{TEMPERATURE_COL}{s}"] = \
                gdf_lab_dic[s][f"{TEMPERATURE_COL}{s}"].subtract(df_meteo_ref)\
                    .dropna()
        elif VAR == "T":
            gdf_uhi_dic[s][f"{TEMPERATURE_COL}{s}"] = \
                gdf_lab_dic[s][f"{TEMPERATURE_COL}{s}"].copy(deep=True)
        gdf_uhi_dic[s].rename({f"{TEMPERATURE_COL}{s}" : VAR}, 
                              axis = 1).dropna().to_file(UHI_DATA_PATH.format(uhi_filename.format(s)))
        
    return 