#GEOCOVID-selectionClusters.py

import pandas as pd
import geopandas as gpd
import numpy as np
import getpass
from sqlalchemy import create_engine
import psycopg2 as ps


outdir='outputs/'


#CONNECT TO DB WITH user=aladoy
pw=getpass.getpass() #Ask for user password
engine=create_engine("postgresql+psycopg2://aladoy:{}@localhost/geocovid".format(pw)) #Create SQLAlchemy engine
conn=ps.connect("dbname='geocovid' user='aladoy' host='localhost' password='{}'".format(pw)) #Create a connection object
cursor=conn.cursor() #Create a cursor object


#IMPORT DATA
#Genomes that were sequenced (cf. mail Damien Jacquot 16/02/2021)
sequenced=pd.read_excel('data/df_merged_GEOCOVID-2.xlsx', engine='openpyxl')
sequenced=sequenced.drop_duplicates('id_demande')[['clusterID','startDate','endDate','observedCases','id_demande']]
sequenced=tuple(sequenced['id_demande'].tolist())
#Extract info from database (location of individuals +++)
sequenced=gpd.read_postgis("SELECT id_demande, date_reception, charge_virale, cat_charge_virale,geometry FROM covid_tests_vd WHERE id_demande IN {}".format(sequenced),conn,geom_col='geometry')
#Save to geojson
sequenced.to_file(outdir+'sequenced.geojson',driver='GeoJSON')


#FUNCTION THAT RETURNS CASES WITHIN A SPECIFIC SIGNIFICANT CLUSTER THAT WERE INCLUDED IN SATSCAN ANALYSIS
#Arguments:
#   cluster_id: cluser of interest (must be significant since we didn't save backup SaTScan cases for non significant clusters)
#   save: if we want to save the results in geojson
#Outputs:
#   clusterCasesAnalysis: list of cases that were included in SaTScan analysis for the given cluster
def casesIncludedSaTScan(cluster_id, save=True):
    #Read file
    casesAnalysis=gpd.read_file('../Part1_Description/outputs/COVID_satscan/ByRELI/ByDays/cluster_cases.gpkg')
    casesAnalysis['end_date']=pd.to_datetime(casesAnalysis.end_date)
    #Add unique ID of clusters from the database
    clustIDs=pd.read_sql_query("SELECT id, cluster, end_date FROM satscan_clusters_firstvague WHERE significant=True", conn)
    #Dataframe with significant clusters + within id_demande
    casesAnalysis=pd.merge(casesAnalysis,clustIDs, how='inner', on=['cluster','end_date'])

    #Extract cases within the cluster of interest included in SaTScan analysis
    clusterCasesAnalysis=tuple(casesAnalysis[casesAnalysis.id==cluster_id]['id_demande'].tolist())
    clusterCasesAnalysis=gpd.read_postgis("SELECT * FROM covid_tests_vd WHERE id_demande IN {}".format(clusterCasesAnalysis),conn,geom_col='geometry')

    if save==True:
        clusterCasesAnalysis.to_file(outdir+'casesAnalysis'+str(cluster_id)+'.geojson',driver='GeoJSON')

    return clusterCasesAnalysis


#FUNCTION THAT RETURNS CLUSTERS FOR A GIVEN AREA
#Arguments:
#   loc: desired area (name of the municipality)
#   seq: list of all sequenced individuals
#Outputs:
#   optionsClust: clusters from locClust that contains at least one individual that have been already sequenced + their caracteristics (in order to select the best match)
#   casesClust: cluster and their cases (geometric intersection so not necessarly the cases included in SaTScan analyses)

def clustersLocation(loc, seq):

    #Clusters that intersect the location
    sql="SELECT cl.id, cl.start_date, cl.end_date, cl.significant, cl.observed, cl.duration, cl.geometry, c.id_demande, c.note_geocoding \
        FROM (SELECT id_demande, date_reception, note_geocoding, geometry FROM covid_tests_vd WHERE res_cov=1) as c INNER JOIN \
        (SELECT s.* FROM satscan_clusters_firstvague s, municipalities m WHERE m.name='{}' AND st_intersects(s.geometry, m.geometry)) cl \
        ON st_intersects(cl.geometry, c.geometry) WHERE (c.date_reception>=cl.start_date AND c.date_reception<=cl.end_date)".format(loc)
    casesClust=gpd.read_postgis(sql,conn,geom_col='geometry')

    #Add dummy var for merge
    sequenced['sequenced']=1
    casesClust=pd.merge(casesClust,sequenced[['id_demande','sequenced']], how='left', on='id_demande')
    casesClust['geo_building']=np.where(casesClust.note_geocoding.str.startswith('Geocoded at building.'),1,np.nan)

    #For each cluster, compute the number of cases considered in SaTScan analysis (observed), the number of cases within cluster if we are doing a geometric intersection with residential coordinates (nb_cases)
    # including the number that are geocoded at building (nb_geoBuilding) and the number that have been already sequenced (nb_sequenced)
    optionsClust=casesClust.groupby(['id','duration','significant','observed']).agg({'id_demande':'count','geo_building':'count','sequenced':'count'}).reset_index()
    optionsClust.columns=['id','duration','significant','observed','nb_cases','nb_geoBuilding','nb_sequenced']
    optionsClust=optionsClust.sort_values('nb_sequenced',ascending=False)

    print('Number of potential clusters we can select: ' + str(optionsClust.shape[0]))
    print('Number of sequenced genomes in this area: ' + str(indSequenced.shape[0]))

    # #Clusters that intersect the location
    # locClust=gpd.read_postgis("SELECT s.* FROM satscan_clusters_firstvague s, municipalities m WHERE m.name='{}' AND st_intersects(s.geometry, m.geometry)".format(loc),conn,geom_col='geometry')
    # print('Number of clusters in the first vague within ' + loc + ': ' + str(casesClust.drop_duplicates('id_demande').sequenced.sum()))
    # #Clusters and individuals sequenced
    # seqClust=gpd.sjoin(locClust, sequenced, how="inner", op='contains')
    # #Keep only individuals that are tested during cluster duration
    # seqClust=seqClust[(seqClust.date_reception>=seqClust.start_date) & (seqClust.date_reception<=seqClust.end_date)]
    # #seqClust=pd.merge(seq[['id', 'id_demande']], locClust, how='inner',on='id')
    # optionsClust=seqClust.groupby(['id','duration','significant','observed'])['id_demande'].agg('count').reset_index().sort_values('id_demande',ascending=False)
    # print('Number of potential clusters we can select: ' + str(optionsClust.shape[0]))
    # #Individuals sequenced
    # #indSequenced=seq[seq.id.isin(locClust.id)].drop_duplicates('id_demande')
    # indSequenced=gpd.sjoin(sequenced, locClust[['geometry']], how="inner", op='within').drop_duplicates('id_demande')
    # #Are there individuals already sequenced in these clusters?
    # print('Number of sequenced genomes in this area: ' + str(indSequenced.shape[0]))

    return optionsClust, casesClust



#FUNCTION THAT RETURNS INDIVIDUALS THAT WE MUST SEQUENCED
#Arguments:
#   optionsClust: list of clusters in the desired area that already contain sequenced individuals
#   selected_cluster: id of cluster we selected in the list cluster_options
#   indSequenced: list individuals that were sequenced in the loc clust
#Outputs:
#   ind: list of individuals we must sequenced
#   selectedClust: caracteristics of selectedClust

def ind2sequenced(optionsClust, selected_cluster, indSequenced): #selected_cluster=id of the cluster we keep
    #Caracteristics of the cluster selected
    selectedClust=gpd.read_postgis("SELECT s.* FROM satscan_clusters_firstvague s WHERE s.id={}".format(selected_cluster), conn, geom_col='geometry')

    #Extract nb of individuals to sequenced
    sql="select c.* FROM covid_tests_vd c, satscan_clusters_firstvague s WHERE s.id={} AND (c.date_reception BETWEEN s.start_date AND s.end_date) AND st_within(c.geometry,s.geometry) AND res_cov=1".format(selected_cluster)
    ind=gpd.read_postgis(sql, conn, geom_col='geometry')


    if (optionsClust.empty | indSequenced.empty):
        nbIndSequenced=0
        print('Number of new individuals to sequenced: ' + str(ind.shape[0]))
    else:
        #Number of sequenced individuals contained in the cluster selected
        nbIndSequenced=optionsClust.loc[optionsClust.id==selected_cluster,'id_demande'].values[0]
        print('Number of new individuals to sequenced: ' + str(ind.shape[0]-nbIndSequenced))
        print('Number of sequenced individuals left: ' + str(indSequenced.shape[0]-nbIndSequenced))

    return ind, selectedClust


#NYON
optionsClustNyon, casesClustNyon= clustersLocation('Nyon', sequenced)



optionsClust[optionsClust.observed==optionsClust.nb_cases]

optionsClust.head(20)

optionsClustNyon

#Select a subset of clusters (long with a lot of individuals sequenced) in order to contain all individuals sequenced
#seqNyon[seqNyon.id.isin([58, 794, 336, 797, 736])].drop_duplicates('id_demande')
#Once it's done, save into QGIS for better choice
#locNyon[locNyon.id.isin([58, 794, 336, 797, 736])].to_file(outdir+'nyon.geojson',driver='GeoJSON')

#Keep 1 between cluster n° 736 & n°794: choose n°794 because it contains more individuals sequenced (3 instead of 2)
indNyon,selectedClustNyon=ind2sequenced(optionsClustNyon,794,indSequencedNyon)

selectedClustNyon


#VALLEE DE JOUX (LE CHENIT)
locClustJoux, optionsClustJoux, indSequencedJoux= clustersLocation('Le Chenit', seqAnalysis)

optionsClustJoux
#Cluster choice:  n°514 because max number of already sequenced genomes, longer cluster before it starts to grow a lot
indJoux,selectedClustJoux=ind2sequenced(optionsClustJoux,514,indSequencedJoux)


indJoux


casesAnalysis[casesAnalysis.id==514].shape

selectedClustJoux


casesAnalysis
casesAnalysis[casesAnalysis.id==514].shape[0]
indVJ=casesAnalysis[casesAnalysis.id==514]

indVJ.head()


indVJ['already_sequenced']=np.where(indVJ.id_demande.isin(['D02225','D01614','D01761']), True, False)
indVJ[indVJ.already_sequenced==True]

indVJ=indVJ[['id','id_demande','date_reception','ct','charge_virale','already_sequenced']]

indVJ.shape[0]

indVJ.to_csv(outdir+'cluster_valleeJoux.csv',index=False)

indVJ



#MORGES
locClustMorges, optionsClustMorges, indSequencedMorges= clustersLocation('Morges', sequenced)

optionsClustMorges

selectedClustMorges

#Select 532 because it's the longest at the same location before it becomes not significant
indMorges,selectedClustMorges=ind2sequenced(optionsClustMorges,532,indSequencedMorges)

16+32+13+4

#PULLY
locClustPully, optionsClustPully, indSequencedPully= clustersLocation('Pully', sequenced)

optionsClustPully
#cluster n°207 for the north of Pully: one of the first location with clusters in Vaud, start of the epidemic, duration is long, high number of already sequenced individuals
indPully1,selectedClustPully1=ind2sequenced(optionsClustPully,207,indSequencedPully)
selectedClustPully1

#Select last significant cluster of the first vague
lastClust=gpd.read_postgis("SELECT * FROM satscan_clusters_firstvague WHERE end_date=(SELECT MAX(end_date) FROM satscan_clusters_firstvague WHERE significant=True)",conn, geom_col='geometry')
lastClust

indLast,selectedClustLast=ind2sequenced(False,1683,False)

#MOUDON
locClustMoudon, optionsClustMoudon, indSequencedMoudon= clustersLocation('Moudon', sequenced)








#List & caracteristics of clusters located in Nyon and containing individuals already sequenced
optJoux=seqJoux.groupby(['id','duration','significant','observed'])['id_demande'].agg('count').reset_index().sort_values('id_demande',ascending=False)
print(optJoux.loc[optJoux.id==461,'id_demande'].values[0])


sqlJoux="select c.* FROM covid_tests_vd c, satscan_clusters_firstvague s WHERE s.id=461 AND (c.date_reception BETWEEN s.start_date AND s.end_date) AND st_within(c.geometry,s.geometry) AND res_cov=1"
#Extract nb of individuals if we sequenced this cluster (total=350)
indJoux=gpd.read_postgis(sqlJoux, conn, geom_col='geometry')

print('Number of new individuals to sequenced: ' + str(indJoux.shape[0]-2))
print('Number of sequenced individuals left: ' + str(5-3))


sql="SELECT s1.id FROM (SELECT s.* FROM satscan_clusters_firstvague s, municipalities m WHERE m.name='Nyon' AND st_intersects(s.geometry, m.geometry)) s1, satscan_clusters_firstvague s2 where s2.id=461 AND st_disjoint(s1.geometry, s2.geometry)"
others=pd.read_sql_query(sql, conn)['id']

others
seqJoux[seqJoux.id.isin(others)]
