#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 16:18:48 2025

@author: Jérémy Bernard, EDYTEM and Lab-STICC, CNRS
"""
import pandas as pd
import os

###############################################################################
############### 0. ONLY VARIABLES NEEDED TO BE SET ############################
###############################################################################
# For VeloCliMap thermoparties: city name. Possible solutions:
#   - "NULL"
#   - "Saint-Jean La Poterie"
#   - "Laval"
#   - "Alençon"
#   - "Amiens"
#   - "Lille"
# For VeloClimameter measurements: day dates. Possible solutions: 
#   - "2025-06-27"
#   - ...
#   - "2025-07-06"
MEAS_NAME = "Gand"

# # Start and end time of the measurement campaign in UTC
# campaign_start = datetime.datetime(2025,7,3,18,35,0, tzinfo = datetime.timezone.utc)
# campaign_end = datetime.datetime(2023,6,12,21,33,39, tzinfo = datetime.timezone.utc)


###############################################################################
############### 1. VARIABLES FOR SENSITIVITY TESTS ############################
###############################################################################
# Whether or not we interpolate linearly the ref if no data
INTERP_REF = True

# Name of the variable of interest ("deltaT", "T")
VAR = "deltaT"

# Averaging method chosen ("lag" or "moving_avg")
AVG_METH = "moving_avg"

# Time spent in a grid cell before taking into account the data / 
# windows size for moving average (s)
TIME_LAG = 10

# Percentiles to use for min and max values per grid cell
PVAL_MAX = 0.9
PVAL_MIN = 0.1

# Path of the file to use as grid / else size of the grid cells used to average the data (m)
GRID_INFO = 0

# SRID code of the grid used
GRID_SRID = 2154

###############################################################################
############### 2. ALREADY SET VARIABLES ######################################
###############################################################################
# Database config
JDBC_DRIVER_VAR = "/home/decide/.m2/repository/org/postgresql/postgresql/42.6.0/postgresql-42.6.0.jar"
DB_INFOS = pd.read_csv("/home/decide/Code/DB_informations.csv",
                       header = 0, 
                       index_col = 0)
HOST = DB_INFOS.loc["veloclimat", "host"]
DBNAME = DB_INFOS.loc["veloclimat", "dbname"]
USERNAME = DB_INFOS.loc["veloclimat", "username"]
PASSWORD = DB_INFOS.loc["veloclimat", "password"]

# Schema and table names in the database
SCHEMA_MET = "veloclimat"
SCHEMA_GEO = "Mapuce"
LABSTICC_SENSOR_TAB = "labsticc_sensor"
VELOCLIMATMETER_TAB = "veloclimatmeter"
LCZ_TAB = "grid_indicators_100" # (or "rsu_indicators")

# Column names
SENSOR_NAME_COL = "sensor_name"
TEMPERATURE_COL = "temperature"
THERMO_NAME_COL = "thermo_name"
TIMESTAMP_COL = "timestamp"
GEOM_COL = "the_geom"
LCZ_ID = "ID_ZONE"

# Reference
REFERENCE_NAME = "Reference"

VELOCLIMAMETER = "veloclimameter"
VELOCLIMAP = "veloclimap"

# Extension for GIS files
GIS_EXT = ".fgb"

###############################################################################
############### 3. DECLARE PATHS ##############################################
###############################################################################
INPUT_DATA_PATH = os.path.abspath(os.path.join(os.path.curdir, "Data", "Inputs"))
OUTPUT_DATA_PATH = os.path.abspath(os.path.join(os.path.curdir, "Data", "Outputs"))
FIGURES_PATH = os.path.abspath(os.path.join(os.path.curdir, "Figures"))

LABSTICC_DIR_IN = os.path.join(INPUT_DATA_PATH, MEAS_NAME)
LABSTICC_DIR_OUT = os.path.join(OUTPUT_DATA_PATH, MEAS_NAME)

VELOCLIMAT_RAW_DATA_PATH = os.path.join(INPUT_DATA_PATH, "{0}" + f"{GIS_EXT}")
LABSTICC_RAW_DATA_PATH = os.path.join(LABSTICC_DIR_IN, "{0}" + f"{GIS_EXT}")
UHI_DATA_PATH = os.path.join(LABSTICC_DIR_OUT, "{0}" + f"{GIS_EXT}")
INPUT_RSU_PATH = os.path.join(LABSTICC_DIR_IN, "{0}" + f"{GIS_EXT}")
OUTPUT_RSU_PATH = os.path.join(LABSTICC_DIR_OUT, "{0}" + f"_fin_{AVG_METH}{GIS_EXT}")

UHI_LCZ_FIGURE = os.path.join(FIGURES_PATH, "{0}_uhi_lcz.png")


###############################################################################
###############. 4. MODIFIED VARIABLES ########################################
###############################################################################
# Identify whether VeloClimameter or VeloCliMap data are manipulated 
try:
    MEAS_NAME = pd.to_datetime(MEAS_NAME)
except:
    pass