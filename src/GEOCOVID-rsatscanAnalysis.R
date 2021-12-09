#GEOCOVID-rsatscanAnalysis.R


# LIBRARIES ---------------------------------------------------------------
# Basic
library(plyr)
library(tidyverse)
# Spatial
library(sf)
library(rsatscan)

# IMPORT DATA --------------------------------------------------------------------
geosan_db_dir <- '/mnt/data/GEOSAN/GEOSAN DB/data/'
output_dir <- '/mnt/data/GEOSAN/RESEARCH PROJECTS/GEOCOVID @ CHUV/GEOCOVID-genomics/'
tests <- st_read('../processed_data/s2_geocoded_data/s2_notfiltered_tests_geo_wreli.gpkg')
reli <- st_read(paste0(geosan_db_dir, "STATPOP 2020/processed_data/statpopVD_poly.gpkg"))
reli_date <- read_csv(paste0(output_dir, 'processed_data/s2_reli_date.csv'))


# SATSCAN INITIAL PARAMS ------------------------------------------------------
ssenv <<- new.env(parent=emptyenv())

#options for satscan can be found here https://rdrr.io/cran/rsatscan/src/R/zzz.R or by extract a .rpm
x=c("[Input]",
  ";case data filename",
  "CaseFile=",
  ";control data filename",
  "ControlFile=",
  ";time precision (0=None, 1=Year, 2=Month, 3=Day, 4=Generic)",
  "PrecisionCaseTimes=3",
  ";study period start date (YYYY/MM/DD)",
  "StartDate=2020/3/2",
  ";study period end date (YYYY/MM/DD)",
  "EndDate=2020/07/01",
  ";population data filename",
  "PopulationFile=",
  ";coordinate data filename",
  "CoordinatesFile=",
  ";use grid file? (y/n)",
  "UseGridFile=n",
  ";grid data filename",
  "GridFile=",
  ";coordinate type (0=Cartesian, 1=latitude/longitude)",
  "CoordinatesType=1",
  "",
  "[Analysis]",
  ";analysis type (1=Purely Spatial, 2=Purely Temporal, 3=Retrospective Space-Time, 4=Prospective Space-Time, 5=Spatial Variation in Temporal Trends, 6=Prospective Purely Temporal, 7=Seasonal Temporal)",
  "AnalysisType=4",
  ";model type (0=Discrete Poisson, 1=Bernoulli, 2=Space-Time Permutation, 3=Ordinal, 4=Exponential, 5=Normal, 6=Continuous Poisson, 7=Multinomial, 8=Rank, 9=UniformTime)",
  "ModelType=0",
  ";scan areas (1=High Rates(Poison,Bernoulli,STP); High Values(Ordinal,Normal); Short Survival(Exponential); Higher Trend(Poisson-SVTT), 2=Low Rates(Poison,Bernoulli,STP); Low Values(Ordinal,Normal); Long Survival(Exponential); Lower Trend(Poisson-SVTT), 3=Both Areas)",
  "ScanAreas=1",
  ";time aggregation units (0=None, 1=Year, 2=Month, 3=Day, 4=Generic)",
  "TimeAggregationUnits=1",
  ";time aggregation length (Positive Integer)",
  "TimeAggregationLength=1",
  "",
  "[Output]",
  ";analysis main results output filename",
  "ResultsFile=",
  ";output Google Earth KML file (y/n)",
  "OutputGoogleEarthKML=n",
  ";output shapefiles (y/n)",
  "OutputShapefiles=n",
  ";output cartesian graph file (y/n)",
  "OutputCartesianGraph=n",
  ";output cluster information in ASCII format? (y/n)",
  "MostLikelyClusterEachCentroidASCII=n",
  ";output cluster information in dBase format? (y/n)",
  "MostLikelyClusterEachCentroidDBase=n",
  ";output cluster case information in ASCII format? (y/n)",
  "MostLikelyClusterCaseInfoEachCentroidASCII=n",
  ";output cluster case information in dBase format? (y/n)",
  "MostLikelyClusterCaseInfoEachCentroidDBase=n",
  ";output location information in ASCII format? (y/n)",
  "CensusAreasReportedClustersASCII=n",
  ";output location information in dBase format? (y/n)",
  "CensusAreasReportedClustersDBase=n",
  ";output risk estimates in ASCII format? (y/n)",
  "IncludeRelativeRisksCensusAreasASCII=n",
  ";output risk estimates in dBase format? (y/n)",
  "IncludeRelativeRisksCensusAreasDBase=n",
  ";output simulated log likelihoods ratios in ASCII format? (y/n)",
  "SaveSimLLRsASCII=n",
  ";output simulated log likelihoods ratios in dBase format? (y/n)",
  "SaveSimLLRsDBase=n",
  ";generate Google Maps output (y/n)",
  "OutputGoogleMaps=n",
  "",
  "[Multiple Data Sets]",
  "; multiple data sets purpose type (0=Multivariate, 1=Adjustment)",
  "MultipleDataSetsPurposeType=0",
  "",
  "[Data Checking]",
  ";study period data check (0=Strict Bounds, 1=Relaxed Bounds)",
  "StudyPeriodCheckType=0",
  ";geographical coordinates data check (0=Strict Coordinates, 1=Relaxed Coordinates)",
  "GeographicalCoordinatesCheckType=0",
  "",
  "[Locations Network]",
  ";locations network filename",
  "LocationsNetworkFilename=",
  ";use locations network file",
  "UseLocationsNetworkFile=n",
  ";purpose of locations network file (0=Coordinates File Override, 1=Network Definition)",
  "PurposeLocationsNetworkFile=1",
  "",
  "[Spatial Neighbors]",
  ";use neighbors file (y/n)",
  "UseNeighborsFile=n",
  ";neighbors file",
  "NeighborsFilename=",
  ";use meta locations file (y/n)",
  "UseMetaLocationsFile=n",
  ";meta locations file",
  "MetaLocationsFilename=",
  ";multiple coordinates type (0=OnePerLocation, 1=AtLeastOneLocation, 2=AllLocations)",
  "MultipleCoordinatesType=0",
  "",
  "[Spatial Window]",
  ";maximum spatial size in population at risk (<=50%)",
  "MaxSpatialSizeInPopulationAtRisk=50",
  ";restrict maximum spatial size - max circle file? (y/n)",
  "UseMaxCirclePopulationFileOption=n",
  ";maximum spatial size in max circle population file (<=50%)",
  "MaxSpatialSizeInMaxCirclePopulationFile=50",
  ";maximum circle size filename",
  "MaxCirclePopulationFile=",
  ";restrict maximum spatial size - distance? (y/n)",
  "UseDistanceFromCenterOption=n",
  ";maximum spatial size in distance from center (positive integer)",
  "MaxSpatialSizeInDistanceFromCenter=1",
  ";include purely temporal clusters? (y/n)",
  "IncludePurelyTemporal=n",
  ";window shape (0=Circular, 1=Elliptic)",
  "SpatialWindowShapeType=0",
  ";elliptic non-compactness penalty (0=NoPenalty, 1=MediumPenalty, 2=StrongPenalty)",
  "NonCompactnessPenalty=1",
  ";isotonic scan (0=Standard, 1=Monotone)",
  "IsotonicScan=0",
  "",
  "[Temporal Window]",
  ";minimum temporal cluster size (in time aggregation units)",
  "MinimumTemporalClusterSize=1",
  ";how max temporal size should be interpretted (0=Percentage, 1=Time)",
  "MaxTemporalSizeInterpretation=0",
  ";maximum temporal cluster size (<=90%)",
  "MaxTemporalSize=50",
  ";include purely spatial clusters? (y/n)",
  "IncludePurelySpatial=n",
  ";temporal clusters evaluated (0=All, 1=Alive, 2=Flexible Window)",
  "IncludeClusters=0",
  ";flexible temporal window start range (YYYY/MM/DD,YYYY/MM/DD)",
  "IntervalStartRange=2000/1/1,2000/12/31",
  ";flexible temporal window end range (YYYY/MM/DD,YYYY/MM/DD)",
  "IntervalEndRange=2000/1/1,2000/12/31",
  "",
  "[Cluster Restrictions]",
  ";risk limit high clusters (y/n)",
  "RiskLimitHighClusters=n",
  ";risk threshold high clusters (1.0 or greater)",
  "RiskThresholdHighClusters=1",
  ";risk limit low clusters (y/n)",
  "RiskLimitLowClusters=n",
  ";risk threshold low clusters (0.000 - 1.000)",
  "RiskThresholdLowClusters=1",
  ";minimum cases in low rate clusters (positive integer)",
  "MinimumCasesInLowRateClusters=0",
  ";minimum cases in high clusters (positive integer)",
  "MinimumCasesInHighRateClusters=2",
  "",
  "[Space and Time Adjustments]",
  ";time trend adjustment type (0=None, 2=LogLinearPercentage, 3=CalculatedLogLinearPercentage, 4=TimeStratifiedRandomization, 5=CalculatedQuadratic, 1=TemporalNonparametric)",
  "TimeTrendAdjustmentType=0",
  ";time trend adjustment percentage (>-100)",
  "TimeTrendPercentage=0",
  ";time trend type - SVTT only (Linear=0, Quadratic=1)",
  "TimeTrendType=0",
  ";adjust for weekly trends, nonparametric",
  "AdjustForWeeklyTrends=n",
  ";spatial adjustments type (0=None, 1=SpatiallyStratifiedRandomization, 2=SpatialNonparametric)",
  "SpatialAdjustmentType=0",
  ";use adjustments by known relative risks file? (y/n)",
  "UseAdjustmentsByRRFile=n",
  ";adjustments by known relative risks file name (with HA Randomization=1)",
  "AdjustmentsByKnownRelativeRisksFilename=",
  "",  
  "[Inference]",
  ";p-value reporting type (Default p-value=0, Standard Monte Carlo=1, Early Termination=2, Gumbel p-value=3) ",
  "PValueReportType=0",
  ";early termination threshold",
  "EarlyTerminationThreshold=50",
  ";report Gumbel p-values (y/n)",
  "ReportGumbel=n",
  ";Monte Carlo replications (0, 9, 999, n999)",
  "MonteCarloReps=999",
  ";adjust for earlier analyses(prospective analyses only)? (y/n)",
  "AdjustForEarlierAnalyses=n",
  ";prospective surveillance start date (YYYY/MM/DD)",
  "ProspectiveStartDate=2000/12/31",
  ";perform iterative scans? (y/n)",
  "IterativeScan=n",
  ";maximum iterations for iterative scan (0-32000)",
  "IterativeScanMaxIterations=10",
  ";max p-value for iterative scan before cutoff (0.000-1.000)",
  "IterativeScanMaxPValue=0.05",
  "",
  "[Cluster Drilldown]",
  ";perform detected cluster standard drilldown (y/n)",
  "PerformStandardDrilldown=n",
  ";perform detected cluster Bernoulli drilldown (y/n)",
  "PerformBernoulliDrilldown=n",
  ";minimum number of locations in detected cluster to perform drilldown (positive integer)",
  "DrilldownMinimumClusterLocations=2",
  ";minimum number of cases in detected cluster to perform drilldown (positive integer)",
  "DrilldownMinimumClusterCases=10",
  ";p-value cutoff of detected cluster to perform drilldown (0.000-1.000)",
  "DrilldownClusterPvalueCutoff=0.05",
  ";adjust for weekly trends, purely spatial Bernoulli drilldown",
  "DrilldownAdjustForWeeklyTrends=n",
  "",
  "[Miscellaneous Analysis]",
  ";calculate Oliveira's F",
  "CalculateOliveira=n",
  ";number of bootstrap replications for Oliveira calculation (minimum=100, multiple of 100)",
  "NumBootstrapReplications=1000",
  ";p-value cutoff for cluster's in Oliveira calculation (0.000-1.000)",
  "OliveiraPvalueCutoff=0.05",
  ";frequency of prospective analyses type (0=Same Time Aggregation, 1=Daily, 2=Weekly, 3=Monthy, 4=Quarterly, 5=Yearly)",
  "ProspectiveFrequencyType=0",
  ";frequency of prospective analyses  (positive integer)",
  "ProspectiveFrequency=1",
  "",
  "[Power Evaluation]",
  ";perform power evaluation - Poisson only (y/n)",
  "PerformPowerEvaluation=n",
  ";power evaluation method (0=Analysis And Power Evaluation Together, 1=Only Power Evaluation With Case File, 2=Only Power Evaluation With Defined Total Cases)",
  "PowerEvaluationsMethod=0",
  ";total cases in power evaluation",
  "PowerEvaluationTotalCases=600",
  ";critical value type (0=Monte Carlo, 1=Gumbel, 2=User Specified Values)",
  "CriticalValueType=0",
  ";power evaluation critical value .05 (> 0)",
  "CriticalValue05=0",
  ";power evaluation critical value .001 (> 0)",
  "CriticalValue01=0",
  ";power evaluation critical value .001 (> 0)",
  "CriticalValue001=0",
  ";power estimation type (0=Monte Carlo, 1=Gumbel)",
  "PowerEstimationType=0",
  ";number of replications in power step",
  "NumberPowerReplications=1000",
  ";power evaluation alternative hypothesis filename",
  "AlternativeHypothesisFilename=",
  ";power evaluation simulation method for power step (0=Null Randomization, 1=N/A, 2=File Import)",
  "PowerEvaluationsSimulationMethod=0",
  ";power evaluation simulation data source filename",
  "PowerEvaluationsSimulationSourceFilename=",
  ";report power evaluation randomization data from power step (y/n)",
  "ReportPowerEvaluationSimulationData=n",
  ";power evaluation simulation data output filename",
  "PowerEvaluationsSimulationOutputFilename=",
  "",  
  "[Spatial Output]",
  ";automatically launch map viewer - gui only (y/n)",
  "LaunchMapViewer=y",
  ";create compressed KMZ file instead of KML file (y/n)",
  "CompressKMLtoKMZ=n",
  ";whether to include cluster locations kml output (y/n)",
  "IncludeClusterLocationsKML=y",
  ";threshold for generating separate kml files for cluster locations (positive integer)",
  "ThresholdLocationsSeparateKML=1000",
  ";report hierarchical clusters (y/n)",
  "ReportHierarchicalClusters=y",
  ";criteria for reporting secondary clusters(0=NoGeoOverlap, 1=NoCentersInOther, 2=NoCentersInMostLikely,  3=NoCentersInLessLikely, 4=NoPairsCentersEachOther, 5=NoRestrictions)",
  "CriteriaForReportingSecondaryClusters=0",
  ";report gini clusters (y/n)",
  "ReportGiniClusters=n",
  ";gini index cluster reporting type (0=optimal index only, 1=all values)",
  "GiniIndexClusterReportingType=0",
  ";spatial window maxima stops (comma separated decimal values[<=50%] )",
  "SpatialMaxima=1,2,3,4,5,6,8,10,12,15,20,25,30,40,50",
  ";max p-value for clusters used in calculation of index based coefficients (0.000-1.000)",
  "GiniIndexClustersPValueCutOff=0.05",
  ";report gini index coefficents to results file (y/n)",
  "ReportGiniIndexCoefficents=n",
  ";restrict reported clusters to maximum geographical cluster size? (y/n)",
  "UseReportOnlySmallerClusters=n",
  ";maximum reported spatial size in population at risk (<=50%)",
  "MaxSpatialSizeInPopulationAtRisk_Reported=50",
  ";restrict maximum reported spatial size - max circle file? (y/n)",
  "UseMaxCirclePopulationFileOption_Reported=n",
  ";maximum reported spatial size in max circle population file (<=50%)",
  "MaxSizeInMaxCirclePopulationFile_Reported=50",
  ";restrict maximum reported spatial size - distance? (y/n)",
  "UseDistanceFromCenterOption_Reported=n",
  ";maximum reported spatial size in distance from center (positive integer)",
  "MaxSpatialSizeInDistanceFromCenter_Reported=1",
  "",
  "[Temporal Output]",
  ";output temporal graph HTML file (y/n)",
  "OutputTemporalGraphHTML=n",
  ";temporal graph cluster reporting type (0=Only most likely cluster, 1=X most likely clusters, 2=Only significant clusters)",
  "TemporalGraphReportType=0",
  ";number of most likely clusters to report in temporal graph (positive integer)",
  "TemporalGraphMostMLC=1",
  ";significant clusters p-value cutoff to report in temporal graph (0.000-1.000)",
  "TemporalGraphSignificanceCutoff=0.05",
  "",
  "[Other Output]",
  ";report critical values for .01 and .05? (y/n)",
  "CriticalValue=n",
  ";report cluster rank (y/n)",
  "ReportClusterRank=n",
  ";print ascii headers in output files (y/n)",
  "PrintAsciiColumnHeaders=n",
  ";user-defined title for results file",
  "ResultsTitle=",
  "",
  "[Elliptic Scan]",
  ";elliptic shapes - one value for each ellipse (comma separated decimal values)",
  "EllipseShapes=1.5,2,3,4,5",
  ";elliptic angles - one value for each ellipse (comma separated integer values)",
  "EllipseAngles=4,6,9,12,15",
  "",
  "[Power Simulations]",
  ";simulation methods (0=Null Randomization, 1=N/A, 2=File Import)",
  "SimulatedDataMethodType=0",
  ";simulation data input file name (with File Import=2)",
  "SimulatedDataInputFilename=",
  ";print simulation data to file? (y/n)",
  "PrintSimulatedDataToFile=n",
  ";simulation data output filename",
  "SimulatedDataOutputFilename=",
  "",
  "[Run Options]",
  ";number of parallel processes to execute (0=All Processors, x=At Most X Processors)",
  "NumberParallelProcesses=0",
  ";suppressing warnings? (y/n)",
  "SuppressWarnings=n",
  ";log analysis run to history file? (y/n)",
  "LogRunToHistoryFile=n",
  ";analysis execution method  (0=Automatic, 1=Successively, 2=Centrically)",
  "ExecutionType=0",
  "",
  "[System]",
  ";system setting - do not modify", "Version=10.0.1"
)

ssenv$.ss.params.defaults = x
ssenv$.ss.params = x



# DATA WRANGLING ----------------------------------------------------------

# Convert attribute's type and remove geometry
tests <- tests %>% mutate(date_reception = as.Date(date_reception)) %>% st_drop_geometry()

# Filter for the first study period
start <- as.Date("2020/03/02",format="%Y/%m/%d")
end <- as.Date("2020/06/30",format="%Y/%m/%d")
tests <- tests %>% filter(between(date_reception, start, end))

# LOOP --------------------------------------------------------------------

# Date of propesctive analysis
theDate <- as.Date("2020/03/10",format="%Y/%m/%d")

# Filter tests
tests_analysis <- tests %>% filter(between(date_reception, start, theDate))
tests_analysis <- tests_analysis %>% filter(institution == 'NO')

# Filter geographic locations
reli_date_analysis <- reli_date %>% filter(between(date, start, theDate))

# Filter duplicates
duplicates <- tests_analysis %>%filter(duplicated(.[["id_patient_study2"]])) %>% distinct(id_patient_study2)

for(i in 1:nrow(duplicates)){
  id_patient = duplicates[i,]
  subs = tests_analysis %>% filter(id_patient_study2 == id_patient) %>% select(id_demande_study2,id_patient_study2,date_reception, res_cov)
  if(all(subs$res_cov==0)){ # if patient was never tested positive, keep most recent test
    id_demande_to_keep = subs %>% arrange(date_reception) %>% do(tail(., n=1)) %>% pull(id_demande_study2)
  }
  else{ # if patient was ever tested positive, keep first positive test
    id_demande_to_keep = subs %>% filter(res_cov == 1) %>% slice(which.min(date_reception)) %>% pull(id_demande_study2)
  }
  # Remove other rows corresponding to the same patient from the initial dataframe
  tests_analysis <- tests_analysis %>% filter((id_patient_study2 != id_patient) | ((id_patient_study2 == id_patient) & (id_demande_study2 == id_demande_to_keep)))
}

# Count daily number of tests and cases per RELI
nb_per_reli <- tests_analysis %>% group_by(reli, date_reception) %>% summarise(nb_cases=sum(res_cov), nb_tests=n())

# Merge
data_analysis <- reli_date_analysis %>% left_join(nb_per_reli, by=c('reli'='reli', 'date'='date_reception'))
# Fill NA with 0
data_analysis <- data_analysis %>% mutate(nb_cases=replace_na(nb_cases,0), nb_tests=replace_na(nb_tests,0))

# Create graph for perc_positivity
stats <- data_analysis %>% group_by(date) %>% summarise(nb_cases=sum(nb_cases), nb_tests=sum(nb_tests))
stats <- stats %>% mutate(perc_positivity=nb_cases*100/nb_tests)
ggplot(data=stats, aes(x=date, y=perc_positivity)) + geom_line() + scale_x_date(date_breaks = "2 days", date_labels = "%d/%m")

# Define parameters
satscan_usr_param = list(
  CaseFile="cases.cas",
  PopulationFile="population.pop",
  CoordinatesFile="geo.geo",
  PrecisionCaseTimes="3",#3=Day
  StartDate=format(as.Date(start),"%Y/%-m/%-d"),
  EndDate=format(as.Date(end),"%Y/%-m/%-d"),
  CoordinatesType=1, #1=Lat/Lon
  AnalysisType=3, #4=Prospective Space-Time (3=Retrospective Space-Time)
  ModelType=0, #0=Discrete Poisson
  ScanAreas=1, #1=High rates
  ResultsFile=paste0(output_dir,'res.txt'),
  MaxSpatialSizeInPopulationAtRisk=10, #10%
  SpatialWindowShapeType=0, #0=Circles,
  SpatialAdjustmentType=0,
  MaxTemporalSize=7, #50% (14 days)
  MaxTemporalSizeInterpretation=1, #0=Percentage (1=Time)
  TimeAggregationUnits=3, #day
  TimeAggregationLength=1, 
  MaxSpatialSizeInPopulationAtRisk_Reported=1, #1%
  OutputTemporalGraphHTML='y',
  MonteCarloReps=999,
  AdjustForEarlierAnalyses='n',
  MinimumCasesInHighRateClusters=5, #Min 3 cases
  MinimumTemporalClusterSize=3,
  CriteriaForReportingSecondaryClusters=1, #No cluster centers in other clusters
  Version='10.0.1')

#Filter cases and controls for the given time period
cases_analysis <- data_analysis %>% filter(nb_cases > 0) %>% mutate(date=format(as.Date(date),"%Y/%-m/%-d")) %>% select(reli, nb_cases, date)
pop_analysis <- data_analysis %>% mutate(date=format(as.Date(date),"%Y/%-m/%-d"))  %>% select(reli, date, nb_tests)
geo_analysis <- data_analysis %>% distinct(reli, geometry) %>% st_as_sf(wkt='geometry')
st_crs(geo_analysis) <- 2056 
geo_analysis <- geo_analysis %>% mutate(geometry=st_transform(geometry, 4326)) %>% mutate(x = sf::st_coordinates(.)[,1], y = sf::st_coordinates(.)[,2])


# library(tmap)
# current.mode <- tmap_mode("view")
# tm_basemap(leaflet::providers$Stamen.Watercolor) + tm_shape(geo_analysis) + tm_dots(col = "red")

#Create temporary directory and store csv for satscan
td = tempdir()
write.cas(as.data.frame(cases_analysis),td, "cases")
write.pop(as.data.frame(pop_analysis),td,"population")
write.geo(as.data.frame(geo_analysis %>% st_drop_geometry()), td, "geo")

#Reset parameters
invisible(ss.options(reset=TRUE))
#print(satscan_usr_param)
ss.options(satscan_usr_param)
write.ss.prm(td, "geocovid")

#Run analysis
res = satscan(td,"geocovid",sslocation="/home/aladoy/SaTScan")
#Extract clusters shape 
clusters=st_as_sf(res$shapeclust)
#Save text file
capture.output(summary(res), file=paste0(dir,theDate,'_results.txt'),append=FALSE)
#Save to shapefile
st_write(clusters,paste0(dir,theDate,'_results.geojson'),driver='GeoJSON')



# Exclude 


# - Download satscan version v10
# - Add a non-parametrical spatial adjustment to detect localized temporal increases

