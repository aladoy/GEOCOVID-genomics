# GEOCOVID-prepareSatscanData.py

'''
This code prepares the three text files that are required for the SaTScan
analysis (Poisson model), namely the case file containing COVID-19 tests (.cas)
, the population file (.pop), and the geographic file (.geo) containing the
points coordinates.
'''

# LIBRARIES
import os
import pandas as pd
import geopandas as gpd
import sys
import importlib

# Import functions from GIRAPH-functions repository
sys.path.append(r'/mnt/data/GEOSAN/FUNCTIONS/GIRAPH-functions/')
try:
    import basic_utils as b
except FileNotFoundError:
    print("Wrong file or file path")

# DIRECTORIES
project_dir: str = (r"/mnt/data/GEOSAN/RESEARCH PROJECTS/GEOCOVID @ CHUV/"
                    "GEOCOVID-genomics")

geosan_db_dir: str = (r"/mnt/data/GEOSAN/GEOSAN DB/data")

output_dir = os.sep.join([project_dir, 'processed_data/satscan_models/model2'])


# IMPORT DATA
# Hectometric grid
reli = gpd.read_file(os.sep.join(
    [geosan_db_dir, "STATPOP 2020/processed_data/statpopVD_poly.gpkg"]),
                     crs=2056, geometry='geometry')
# Lower column names
reli.columns = reli.columns.str.lower()
# Extract centroids and total population only
reli_centroids = reli[['reli', 'b20btot', 'geometry']]
reli_centroids['geometry'] = reli_centroids.geometry.centroid


# RT-PCR tests
tests = gpd.read_file(os.sep.join(
    [project_dir, "processed_data/s2_geocoded_data",
     "s2_notfiltered_tests_geo.gpkg"]),  # change this line to use another file
                      crs=2056, geometry='geometry')
tests['date_reception'] = tests.date_reception.map(pd.to_datetime)
# filter for study 1 only
tests = tests[tests.date_reception < '2020-07-01']

# MATCH TESTS WITH RELI (closest RELI if outside)
tests_reli = tests.sjoin(reli[['reli', 'geometry']],
                         how='left', op='intersects')
print('Number of tests: ', tests.shape[0])
print('Number of tests after the join: ', tests_reli.shape[0])
print('Number of tests that are located exactly between two RELIs: ',
      tests_reli[tests_reli.duplicated(
          subset='id_demande_study2', keep=False)].shape[0])
print('Number of tests that fall outside RELIs: ',
      tests_reli[tests_reli.reli.isna()].shape[0])

# for location in the border line of two relis, just keep the first occurence
tests_reli.drop_duplicates(subset='id_demande_study2',
                           keep='first', inplace=True)

# use nearest neighbor for outside locations
importlib.reload(b)
# tests_reli.loc[tests_reli.reli.isna()].to_file('outside_RELI.gpkg', driver='GPKG')
tests_reli.loc[tests_reli.reli.isna(), 'reli'] = tests_reli[tests_reli.reli.isna()].apply(
        lambda row: b.find_nearest_reli(row, reli_centroids, radius=500),
        axis=1)
# remove tests that are > 500 meters away from the nearest reli
print('Nb of tests that were removed because too far away from the next RELI:',
      tests_reli.loc[tests_reli.reli.isna()].shape[0])
tests_reli = tests_reli[~tests_reli.reli.isna()]
tests_reli['reli'] = tests_reli.reli.astype('int')

# COMPUTE DAILY NUMBER OF CASES PER RELI
# extract only cases
cases = tests_reli[tests_reli.res_cov == 1]

# compute daily number of cases
cases_agg = cases.groupby(by=['date_reception', 'reli']).size()
cases_agg = cases.groupby(
    by=['date_reception', 'reli'])['id_demande_study2'].agg(
        'count').reset_index(drop=False)
cases_agg.rename(columns={'id_demande_study2': 'nb_cases',
                          'date_reception': 'date'}, inplace=True)
# verify we have the same number of cases
print(sum(cases_agg.nb_cases) == cases.shape[0])

# PARAMETERS - Define the software that will be used
# (rsatscan requires that population data contain dates while satscan does not)
usr_soft = input(
    ('Choose the software that will be used to run SaTScan analyses:,'
     '(rsatscan or satscan)'))
if usr_soft in ('rsatscan', 'satscan'):
    soft = usr_soft
else:
    raise ValueError("Invalid input. Please use another value.")


if soft == 'rsatscan':

    # cross join between reli + cases (for rsatscan)
    ha_j = pd.DataFrame(reli_centroids)
    ha_j['key'] = 1
    date = pd.Series(cases_agg.date.unique()).to_frame(name='date')
    date['key'] = 1
    ha_date = pd.merge(ha_j, date, on='key').drop("key", 1)

    # add month / week information to covid data
    cases_agg['month'] = cases_agg.date.dt.month  # add month
    cases_agg['week'] = cases_agg.date.dt.isocalendar().week  # add week

    # add lat / lon to ha
    reli_centroids['x'] = reli_centroids.to_crs("epsg:4326").geometry.x
    reli_centroids['y'] = reli_centroids.to_crs("epsg:4326").geometry.y

    # split date between cases,at-risk population and geo
    population = ha_date[['reli', 'date', 'b20btot']]
    geo = reli_centroids[['reli', 'x', 'y']]

    # # add covid data
    # satscan = pd.merge(ha_date, cases_agg, on=['reli', 'date'], how='left')
    # satscan['nb_cases'] = satscan.nb_cases.fillna(value=0)
    # satscan = satscan.convert_dtypes()  # convert to float to int
    # satscan['month'] = satscan.date.dt.month  # add month
    # satscan['week'] = satscan.date.dt.isocalendar().week
    # satscan = gpd.GeoDataFrame(satscan, crs="EPSG:4326", geometry='geometry')
    # satscan['x'] = satscan.geometry.x
    # satscan['y'] = satscan.geometry.y

    # # add month / week information to covid data
    # cases_agg['month'] = cases_agg.date.dt.month  # add month
    # cases_agg['week'] = cases_agg.date.dt.isocalendar().week  # add week

    # # add lat / lon to ha
    # reli_centroids['x'] = reli_centroids.geometry.x

    # # split date between cases,at-risk population and geo
    # cases_agg = satscan[
    #     satscan.nb_cases > 0][['reli', 'date', 'nb_cases', 'month', 'week']]
    # population = satscan[['reli', 'date', 'b20btot']]
    # geo = reli_centroids[['reli', 'x', 'y']]

else:
    # add month / week information to covid data
    cases_agg['month'] = cases_agg.date.dt.month  # add month
    cases_agg['week'] = cases_agg.date.dt.isocalendar().week  # add week

    # add lat / lon to ha
    reli_centroids['x'] = reli_centroids.to_crs("epsg:4326").geometry.x
    reli_centroids['y'] = reli_centroids.to_crs("epsg:4326").geometry.y

    # split date between cases,at-risk population and geo
    cases_agg = cases_agg[['reli', 'date', 'nb_cases', 'month', 'week']]
    population = reli_centroids[['reli', 'b20btot']]
    geo = reli_centroids[['reli', 'x', 'y']]


# SAVE FILES
# save geographical coordinates (same for all the days)
geo.to_csv(os.sep.join([output_dir, 'geo_reli_' + soft + '.csv']), index=False)
# cases
cases_agg.to_csv(os.sep.join(
    [output_dir,  'cases_' + soft + '.csv']), index=False)
# population
population.to_csv(os.sep.join(
    [output_dir,  'pop_' + soft + '.csv']), index=False)
