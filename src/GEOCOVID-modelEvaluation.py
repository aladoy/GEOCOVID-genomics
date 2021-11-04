# LIBRAIRIES
# Basic
import pandas as pd
import os
import re
from sklearn.metrics import rand_score, adjusted_rand_score

# DIRECTORIES
project_dir: str = (r"/mnt/data/GEOSAN/RESEARCH PROJECTS/GEOCOVID @ CHUV/"
                    "GEOCOVID-phase1-genomics")

geosan_db_dir: str = (r"/mnt/data/GEOSAN/GEOSAN DB/data")


# IMPORT DATA
gen_clusters = pd.read_csv(
    os.sep.join([project_dir,
                 'processed_data/genomic_clusters/gen_clusters_02112021.csv']))

gen_clusters.shape


# EXTRACT GEOGRAPHIC CLUSTERS THAT WERE DETECTED IN S1
# s1_clusters_dir = os.sep.join(
#     [project_dir, 'processed_data/cluster_selection',
#      'dossier_damien_18_03/individus_clusters_csv'])
#
# filelist = os.listdir(s1_clusters_dir)
#
# # create an empty list to store geo clusters from study 1
# s1_geoclusters_list = []
#
# for file in filelist:
#     # import file
#     df = pd.read_csv(os.sep.join([s1_clusters_dir, file]))
#     # use regex expression to find cluster id in filename
#     cluster_id = re.search(r'(?<=cases_).*(?=.csv)', file).group(0)
#     # extract IDs only
#     df = df[['id_demande', 'id_patient', 'charge_virale']]
#     # add new column with cluster_id
#     df['geo_cluster_s1'] = cluster_id
#     # append cluster to list
#     s1_geoclusters_list.append(df)
#
# # concatenate all the dataframes in the list
# s1_geoclusters = pd.concat(s1_geoclusters_list)

# FIND THE CORRESPONDANCE OF IDs
sequenced_data = pd.read_csv(
    os.sep.join([project_dir, 'data/210820_geocovid_sequenced.csv']))
sequenced_data = sequenced_data[
    ['alias', 'AllClusters', 'First3cases', 'id_demande',
     'ClusterID_all', 'ClusterID1_first3']]


# MERGE DATAFRAMES
sequenced = pd.merge(sequenced_data, gen_clusters, on='alias')

sequenced
# Add geo_cluster for Model 1 (First study)
sequenced['geo_cluster_M1'] = sequenced['ClusterID_all'].fillna(value=sequenced.ClusterID1_first3)
# Extract only ID for genomic clusters
sequenced['gen_cluster'] = sequenced['gen_cluster'].str.extract(r'gen_cluster(.*)')




rand_score(sequenced.gen_cluster, sequenced.geo_cluster_M1)
# Adjusted rand score to use a metric that is independent of the number of clusters and samples
adjusted_rand_score(sequenced.gen_cluster, sequenced.geo_cluster_M1)
sequenced

# s2_geocoded_data = gpd.read_file(
#     os.sep.join([project_dir, 'processed_data',
#                  's2_geocoded_data/s2_notfiltered_tests_geo.gpkg']))
# s2_geocoded_data['date_reception'] = s2_geocoded_data.date_reception.map(
#     pd.to_datetime)
# # Filter for the s1 time period
# matching_keys_cases = s2_geocoded_data.loc[
#     (s2_geocoded_data.date_reception < '2020-07-01')
#     & (s2_geocoded_data.res_cov == 1),
#     ['id_demande_study2', 'id_patient_study2', 'id_demande_study1',
#      'id_patient_study1']]
#
# a=pd.merge(s1_geoclusters, matching_keys_cases, left_on='id_demande', right_on='id_demande_study1')
# pd.merge(a, gen_clusters, left_on='id_demande_study2', right_on='alias')
