# GEOCOVID-distance Matrix.py
"""
Compute the euclidean distance between sequenced individuals.
Return a spatial distance matrix with individuals' ID (alias) as index/headers.
"""

import pandas as pd
from scipy.spatial.distance import squareform, pdist
import os
import geopandas as gpd

# DIRECTORIES
path_dir: str = r"/mnt/data/Thesis/GEOSAN/GEOCOVID @ CHUV/Study1_Genomics/"
data_filename = os.sep.join([path_dir, "data/210820_geocovid_sequenced.csv"])
outputs_filename = os.sep.join([path_dir, "outputs/dist_matrix_sequenced.csv"])

# LOAD DATASET
ind_seq = pd.read_csv(data_filename, sep=',', index_col='alias')
print(ind_seq.shape)
print(ind_seq.dtypes)

# CONVERT LAT/LON TO X/Y
"""
Lat / Lon coordinates must be converted to projected coordinates before
computing euclidean distance matrix.
Otherwise, results will be in degrees and not in meters.
"""
obs = gpd.GeoDataFrame(ind_seq[['lat', 'lon']],
                       geometry=gpd.points_from_xy(ind_seq.lon, ind_seq.lat))
obs.crs = "EPSG:4326"  # Specify CRS (WGS 84)
obs.to_crs("EPSG:2056", inplace=True)  # Convert to Swiss CRS (CH1903+/LV95)
obs['x'], obs['y'] = obs.geometry.x, obs.geometry.y
obs = pd.DataFrame(obs.drop(columns=['lat', 'lon', 'geometry']))

# COMPUTE DISTANCE MATRIX
# Create an array for observations (alias, x, y)
# cdist = squareform(pdist) but pdist function is faster
dist_matrix = squareform(pdist(obs, metric='euclidean'))
# Convert array to DataFrame
dist_df = pd.DataFrame(dist_matrix, index=obs.index, columns=obs.index)

# SAVE RESULTS
dist_df.to_csv(outputs_filename)
