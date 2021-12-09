# GEOCOVID-AssignTestsToRELI.py

'''
This code assign RT-PCR tests to VD hectares (RELI). For tests that do not
intersect a RELI, the code finds the nearest reli in a 500-meters distance.
The code also creates a CSV file with the unique combination of RELI-date for
the entire period of analysis.
'''

# LIBRARIES
import os
import geopandas as gpd
import pandas as pd
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

output_dir = os.sep.join([project_dir, 'processed_data/s2_geocoded_data'])


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


# EXPORT CSV WITH HA
# cross join between reli + cases (for rsatscan)
ha_j = pd.DataFrame(reli_centroids)
ha_j['key'] = 1
date = pd.Series(tests.date_reception.unique()).to_frame(name='date').sort_values('date', ascending=True)
date['key'] = 1
ha_date = pd.merge(ha_j, date, on='key').drop("key", 1)
ha_date.to_csv(os.sep.join(project_dir, 'processed_data/s2_reli_date.csv'))


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

# EXPORT TO GEOPACKAGE FILE
b.save_gdf(output_dir, 's2_notfiltered_tests_geo_wreli.gpkg', tests_reli, 'GPKG')
