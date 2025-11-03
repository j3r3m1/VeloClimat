# VeloClimat data processing
## Data

This repository contains (or can access to) data produced or
 used within the VeloClimat project.

Some data directly derive from the project:
 - VeloClim_App: Air temperature observed during collaborative mapping events
 - VeloClimameter: several micro-meteorological and physical
    parameters observed all along the itinerary
 - Physio: Physiological variables observed for 4 participants

Other data is used to analyze the previous ones:
 - Meteo: Meteorological variables observed at Meteo-France stations
located along the VeloClimat itinerary
 - GeoClimate: indicators describing the geography and urban canopy 
(e.g. Local Climate Zones - LCZ) of the territory
where the measurement have been carried out  intersected by the itinerary
(calculated using the GeoClimate software and the OpenStreetMap - OSM - data)
 - Terrain: Numerical Terrain Model to have the elevation of each
measurement
     


## Research assumptions 
Assumptions about the GeoClimatic relations (GCA):
 - GCA1. LCZ type may have an impact on air temperature and
IR radiation (up and down)
 - GCA2. Air temperature differences between LCZ may vary
within the day


Assumption about the Sensors (SA):
 - SA1. Wt low speed, radiation might have an effect on air temperature 
and IR radiation measurement 
 - SA2. A lag of several seconds may happen on air temperature measurements


Assumptions about effects of climate variables on physiology (C2P):
 - C2P1. Climate variables (air temperature, RH and IR temperature) have a
a significant impact on heart rate or body temperature


## Method of investigation
To account for the potential effect of radiation, lag 
and LCZ effect (spatial AND temporal), there is a need 
to investigate one parameter making sure
the others have very little effect. Below is described
how to minimize the effect of each assumption while analyzing
the effect of an other:
 - GCA1. Gathering results by LCZ type (and using MF stations as references)
 - GCA2. Gathering results by time periods (and using MF stations as references)
 - SA1. Keeping values only for wind speed higher than a certain threshold
 - SA2. Keeping values only when there is a minimum number of seconds spent in a same LCZ

    
The investigations need to be performed in the following order:
 - SA1: For each set of times of the day and LCZ types (and a high number
of seconds spent in a same LCZ), display temp difference
from the MF station for several ranges of speed
   => Set minimum speed
 - SA2: For each set of times of the day and LCZ types and above the speed
previously identified (SA1), display temp difference from the MF station
for several ranges of minimum seconds spent in a same LCZ
   => Set minimum time spent in a same LCZ
 - GCA1 & GCA2: For each time of the day, above the speed identified in SA1 and
when the time spent in a same LCZ is > the threshold set in SA2, display
the temperature difference and compare all LCZ

The C2P1 assumption can be investigated at the end, testing several times to
average physiological and climatic conditions (to account for potential
inertial - lag - effects). Each participant can be analyzed individually. 
This assumption is tested only for VeloClimameter data.



### Code
Several scripts have been developped in order to download
and analyze the VeloClimat data.

It consists in a main script that calls functions (for VeloClimameter or VeloClim'App):
 1. download: To download (and store) the data,
 2. clean: To clean the data (remove duplicate points)
 3. gather: 
    - VeloClimameter: with physiological record with RSU indicators and
elevation
    - VeloClim'App: with RSU indicators and elevation and calculate wind speed
 5. createRef: To create a reference air temperature based on a bilinear
interpolation of the MF stations and an altitude correction for day-time
VeloClimameter data.
 6. tempDifferences: Calculate air temperature differences from the reference
 7. investigations: investigate the assumptions described above in the order 
they are listed in the method section.

