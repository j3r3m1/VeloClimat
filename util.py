#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 00:02:59 2025

@author: decide
"""

from digitemp.master import UART_Adapter
from digitemp.device import DS18B20
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon

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