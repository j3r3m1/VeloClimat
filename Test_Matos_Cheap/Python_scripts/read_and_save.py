#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 23:20:32 2025

@author: decide
"""

import time
import pandas as pd
import util
import threading
import os
import datetime

# File name
filename = "time_constant_v_0m_s.csv"

# Time step in seconds
dt = 3

# Adresses where may be one of the 2 sensors
adress = "/dev/ttyUSB{0}"

# Sensor ROM ID
rom_ids = ["2800028D10000093", "28B2968A100000D4"]

def poll():
    df = pd.DataFrame(columns = rom_ids)
    start_time = time.time()
    i = 0
    while not stop_flag.is_set():        
        t = i * dt
        s = 0
        sread = 0
        while not ((not s < 10) or (not sread < 2)):
            rom_id, temp = util.read_temp(adress.format(s))
            if temp:
                sread += 1
                df.loc[datetime.datetime.now(), rom_id] = temp
            s += 1
        print(t)
        # Go to next step once dt seconds are elapsed
        elapsed = time.time() - start_time
        time.sleep(max(0, t + dt - elapsed))
        i += 1
    
    df.to_csv(os.path.join(os.path.abspath(os.path.curdir), "Data", filename))
    
stop_flag = threading.Event()
t = threading.Thread(target=poll)
t.start()
        
        
try:
    while True:
        time.sleep(1)

except KeyboardInterrupt:
    print("Stopping...")
    stop_flag.set()
    t.join()
    print("Stopped.")


