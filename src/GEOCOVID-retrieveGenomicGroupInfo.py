'''This code is retrieving information including residential address, age, sex,
viral load from the individuals constituting the genomic groups.
Input files: genomic groups as indicated in the Map Spanning Tree (file was
sent by Yangji Choi)
Output files: individuals within geonomic groups (sequenced.gpkg) geocoded
based on the first the first (layer s1) and the second (layer s2) geocoding
algorithm.
'''

# LIBRAIRIES
import pandas as pd
import os
# from sklearn.metrics import rand_score, adjusted_rand_score, f1_score
import sys
import geopandas as gpd

# DIRECTORIES
project_dir: str = (r"/mnt/data/GEOSAN/RESEARCH PROJECTS/GEOCOVID @ CHUV/"
                    "GEOCOVID-genomics")

geosan_db_dir: str = (r"/mnt/data/GEOSAN/GEOSAN DB/data")


# CONNECT TO GEOSAN DB
sys.path.append(r'/mnt/data/GEOSAN/FUNCTIONS/GIRAPH-functions/')
try:
    import db_utils as db
except FileNotFoundError:
    print("Wrong file or file path")
engine, conn, cursor = db.connect_db('geosan', 'aladoy')

# IMPORT DATA
# Genomic groups for the whole cluster design
gen_clusters = pd.read_csv(
    os.sep.join([project_dir,
                 'processed_data/genomic_clusters/Geocovid_genclusters.csv']))
#   remove individuals that have unique genomic sequence
gen_clusters = gen_clusters[gen_clusters.Gen_Group != 'UNIQUE']
#   delete the second column (keep the new index for genomic group)
gen_clusters.drop('gen_cluster', axis=1, inplace=True)
print('Size of the input file: ', gen_clusters.shape)

# # Genomic clusters
# gen_clusters = pd.read_csv(
#     os.sep.join([project_dir,
#                  'processed_data/genomic_clusters/gen_clusters_02112021.csv']))
# gen_clusters.shape
#
# # Removed the ones that do not appear in the next CSV file
# gen_clusters_all = pd.read_csv(
#     os.sep.join([project_dir,
#                  'processed_data/genomic_clusters/Geocovid_alias_220202.csv']))
# gen_clusters = gen_clusters[gen_clusters.alias.isin(gen_clusters_all.alias)]
#
# gen_clusters.shape
# # Compute number of individuals per genomic group
# nb_per_group = pd.DataFrame(gen_clusters.groupby(
#     'gen_cluster').size().reset_index(drop=False))
# nb_per_group.columns = ['gen_cluster', 'size']
# # Merge with initial df
# gen_clusters = gen_clusters.merge(nb_per_group, how='inner', on='gen_cluster')
#
# # Remove genomic cluster having only
# gen_clusters = gen_clusters[gen_clusters['size'] != 1]

# FIND THE CORRESPONDANCE OF IDs
sequenced_data = pd.read_csv(
    os.sep.join([project_dir, 'data/210820_geocovid_sequenced.csv']))
sequenced_data = sequenced_data[
    ['alias', 'AllClusters', 'First3cases', 'id_demande',
     'ClusterID_all', 'ClusterID1_first3']]

# Extract RT-PCR tests data (for the first study)
sql = "SELECT * \
FROM geocovid.s1_notfiltered_tests_geo \
WHERE res_cov=1 "
# Extract data from database
cases_s1 = gpd.read_postgis(sql, conn, geom_col='geometry')
cases_s1.shape

# Extract RT-PCR tests data (for the second study but in the first period)
sql = "SELECT * \
FROM geocovid.s2_notfiltered_tests_geo_wreli_p1 \
WHERE res_cov=1 "
# Extract data from database
cases_s2 = gpd.read_postgis(sql, conn, geom_col='geometry')
cases_s2.shape

# MERGE DATAFRAMES
sequenced = pd.merge(sequenced_data, gen_clusters, on='alias')
sequenced.shape

# Add geo_cluster for Model 1 (First study)
sequenced['geo_cluster_M1'] = sequenced['ClusterID_all'].fillna(
    value=sequenced.ClusterID1_first3)
# Extract the genomic group index
sequenced['Gen_Group'] = sequenced['Gen_Group'].str.extract(
    r'oup(.*)', expand=False).str.lstrip(to_strip='0')
# sequenced['gen_cluster'] = sequenced['gen_cluster'].str.extract(
#     r'gen_cluster(.*)').astype(int)

# Create a file for the location of sequenced inviduals (related to S1)
seq_s1 = pd.merge(sequenced, cases_s1, how='inner', on='id_demande')
seq_s1.shape

# Create a file for the location of sequenced inviduals (related to S2)
seq_s2 = pd.merge(sequenced, cases_s2, how='inner', left_on='id_demande',
                  right_on='id_demande_study1')
seq_s2.shape
# Fewer rows than seq_s1 because geocoding has been improved and some (wrongly)
# geocoded individuals have been removed.

# convert to geodataframe
seq_s1 = gpd.GeoDataFrame(seq_s1, geometry='geometry', crs='epsg:2056')
seq_s2 = gpd.GeoDataFrame(seq_s2, geometry='geometry', crs='epsg:2056')

# save files
seq_s1.to_file(os.sep.join([project_dir, 'results/sequenced.gpkg']),
               layer='s1', driver='GPKG')
seq_s2.to_file(os.sep.join([project_dir, 'results/sequenced.gpkg']),
               layer='s2', driver='GPKG')


# FIND PERCENTAGE OF SEQUENCED BY CLUSTER
# Import file
gen_clusters_alias = pd.read_csv(
    os.sep.join([project_dir,
                 'processed_data/genomic_clusters/Geocovid_alias_220202.csv']))
# Retrieve id_demande
gen_clusters_alias = pd.merge(gen_clusters_alias,
                              sequenced_data, how='left', on='alias')
# Retrieve individuals information from Study1
gen_clusters_alias = pd.merge(
    gen_clusters_alias, cases_s1, how='inner', on='id_demande')
# convert to geodataframe
gen_clusters_alias = gpd.GeoDataFrame(gen_clusters_alias, geometry='geometry',
                                      crs='epsg:2056')
# save file
gen_clusters_alias.to_file(os.sep.join(
    [project_dir, 'results/sequenced_whole_clusters.gpkg']), driver='GPKG')
# compute number per geographic cluster
gen_clusters_alias.groupby('ClusterID_all').size()


gen_clusters_alias
# FIND THE CORRESPONDANCE OF IDs
sequenced_data = pd.read_csv(
    os.sep.join([project_dir, 'data/210820_geocovid_sequenced.csv']))
sequenced_data = sequenced_data[
    ['alias', 'AllClusters', 'First3cases', 'id_demande',
     'ClusterID_all', 'ClusterID1_first3']]

# Extract RT-PCR tests data (for the first study)
sql = "SELECT * \
FROM geocovid.s1_notfiltered_tests_geo \
WHERE res_cov=1 "
# Extract data from database
cases_s1 = gpd.read_postgis(sql, conn, geom_col='geometry')
cases_s1.shape

#
# # WHOLE CLUSTER DESIGN SAMPLES WITH CORRESPONDING GENOMIC GROUP (TABLE)
#
# # Import file
# gen_clusters_indiv = pd.read_csv(
#     os.sep.join([project_dir,
#                  'processed_data/genomic_clusters/Geocovid_genclusters.csv']))
# gen_clusters_indiv.shape
#
# # Remove individuals with unique genomic sequence
# gen_clusters_indiv = gen_clusters_indiv[
#     gen_clusters_indiv.Gen_Group != 'UNIQUE']
#
# gen_clusters_indiv.shape
#
# gen_clusters_indiv
#
#
# a=seq_s1[['id_demande', 'geometry']].merge(seq_s2[['id_demande', 'geometry']],how='inner',on='id_demande')
# a[(a.gkode_x == a.gkode_y) & (a.gkodn_x == a.gkodn_y)].shape
#
# a[(a.geometry_y.st_buffer]
#
# a
#
# seq_s2
# seq_s1
#
# # EXTRACT GEOGRAPHIC CLUSTERS THAT WERE DETECTED IN S1
# # s1_clusters_dir = os.sep.join(
# #     [project_dir, 'processed_data/cluster_selection',
# #      'dossier_damien_18_03/individus_clusters_csv'])
# #
# # filelist = os.listdir(s1_clusters_dir)
# #
# # # create an empty list to store geo clusters from study 1
# # s1_geoclusters_list = []
# #
# # for file in filelist:
# #     # import file
# #     df = pd.read_csv(os.sep.join([s1_clusters_dir, file]))
# #     # use regex expression to find cluster id in filename
# #     cluster_id = re.search(r'(?<=cases_).*(?=.csv)', file).group(0)
# #     # extract IDs only
# #     df = df[['id_demande', 'id_patient', 'charge_virale']]
# #     # add new column with cluster_id
# #     df['geo_cluster_s1'] = cluster_id
# #     # append cluster to list
# #     s1_geoclusters_list.append(df)
# #
# # # concatenate all the dataframes in the list
# # s1_geoclusters = pd.concat(s1_geoclusters_list)
#
#
# # data = {'gen_cluster':[1,2,2,1,2], 'geo_cluster_M1':[4,4,6,2,6], 'geo_cluster_M2':[1,2,2,2,3], 'geo_cluster_M3':[1,2,2,1,2]}
# # dataset=pd.DataFrame(data)
#
# # f1_score(dataset.gen_cluster, dataset.geo_cluster_M1, average='samples')
# # rand_score(dataset.gen_cluster, dataset.geo_cluster_M1)
# #
# # dataset
#
# rand_score(sequenced.gen_cluster, sequenced.geo_cluster_M1)
# # Adjusted rand score to use a metric that is independent of the number of clusters and samples
# adjusted_rand_score(sequenced.gen_cluster, sequenced.geo_cluster_M1)
# sequenced
#
# # s2_geocoded_data = gpd.read_file(
# #     os.sep.join([project_dir, 'processed_data',
# #                  's2_geocoded_data/s2_notfiltered_tests_geo.gpkg']))
# # s2_geocoded_data['date_reception'] = s2_geocoded_data.date_reception.map(
# #     pd.to_datetime)
# # # Filter for the s1 time period
# # matching_keys_cases = s2_geocoded_data.loc[
# #     (s2_geocoded_data.date_reception < '2020-07-01')
# #     & (s2_geocoded_data.res_cov == 1),
# #     ['id_demande_study2', 'id_patient_study2', 'id_demande_study1',
# #      'id_patient_study1']]
# #
# # a=pd.merge(s1_geoclusters, matching_keys_cases, left_on='id_demande', right_on='id_demande_study1')
# # pd.merge(a, gen_clusters, left_on='id_demande_study2', right_on='alias')
