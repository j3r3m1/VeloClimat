#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 00:02:59 2025

@author: decide
"""

from digitemp.master import UART_Adapter
from digitemp.device import DS18B20

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