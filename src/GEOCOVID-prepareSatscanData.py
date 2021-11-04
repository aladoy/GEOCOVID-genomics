# GEOCOVID-prepareSatscanData.py

'''
This code prepares the three text files that are required for the SaTScan
analysis (Poisson model), namely the case file containing COVID-19 tests (.cas)
, the population file (.pop), and the geographic file (.geo) containing the
points coordinates.
'''

import os
import pandas as pd

# DIRECTORIES
path_dir: str = (r"/mnt/data/GEOSAN/RESEARCH PROJECTS/GEOCOVID @ CHUV/"
                 "GEOCOVID-phase2")

covid_data_filename = os.sep.join(
    [path_dir, "processed_data/s2_filtered_tests_geo_14d.feather"])

reli_data_filename = os.sep.join([path_dir, "data/NPA/PLZO_CSV_LV95.csv"])

pop_data_filename = os.sep.join(
    [path_dir, "processed_data/s2_notfiltered_tests_geo.gpkg"])

output_data_dir = os.sep.join([path_dir, "processed_data/rsatscan/"])

# PARAMETERS - Define the software that will be used
# (rsatscan requires that population data contain dates while satscan does not)
usr_soft = input(
    ('Choose the software that will be used to run SaTScan analyses:,'
     '(rsatscan or satscan)'))
if usr_soft in ('rsatscan', 'satscan'):
    soft = usr_soft
else:
    raise ValueError("Invalid input. Please use another value.")

# LOAD FILES
ha =



if soft == 'rsatscan':
    # CROSS JOIN BETWEEN RELI + COVID (FOR RSATSCAN)
    ha_j = pd.DataFrame(ha)
    ha_j['key'] = 1
    date = pd.Series(covid.date.unique()).to_frame(name='date')
    date['key']=1
    ha_date = pd.merge(ha_j, date, on ='key').drop("key", 1)
    # #Add covid data
    satscan=pd.merge(ha_date,covid,on=['reli','date'],how='left')
    satscan['ncov_tests']=satscan.ncov_tests.fillna(value=0)
    satscan['ncov_cases']=satscan.ncov_cases.fillna(value=0)
    satscan=satscan.convert_dtypes() #Convert to float to int
    satscan['month']=satscan.date.dt.month #Add month
    satscan['week']=satscan.date.dt.isocalendar().week
    satscan=gpd.GeoDataFrame(satscan,crs="EPSG:4326", geometry='geometry')
    satscan['x']=satscan.geometry.x
    satscan['y']=satscan.geometry.y

    #Add month / week information to covid data
    covid['month']=covid.date.dt.month #Add month
    covid['week']=covid.date.dt.isocalendar().week #Add week

    #Add lat / lon to ha
    ha['x']=ha.geometry.x
    ha['y']=ha.geometry.y

    #Split date between cases,at-risk population and geo
    cases=satscan[satscan.ncov_cases>0][['reli','date','ncov_cases','month','week']]
    population=satscan[['reli','date','b19btot']]
    geo=ha[['reli','x','y']]

else:
    #Add month / week information to covid data
    covid['month']=covid.date.dt.month #Add month
    covid['week']=covid.date.dt.isocalendar().week #Add week

    #Add lat / lon to ha
    ha['x']=ha.geometry.x
    ha['y']=ha.geometry.y

    #Split date between cases,at-risk population and geo
    cases=covid[covid.ncov_cases>0][['reli','date','ncov_cases','month','week']]
    population=ha[['reli','b19btot']]
    geo=ha[['reli','x','y']]
