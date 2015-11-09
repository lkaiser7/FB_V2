### forest birds SDMs V 2.0 source script ###
### scripts to build and run biomod2 sdms ###
### master code to use in sdm FB_analysis ###

# remove all past worksheet variables
rm(list = ls()) 
# keep code from running after errors 
# options(error = stop) 

# load necessary packages
library(stringr)

###############################
##### SET SOURCE LOCATION #####
###############################

# assign source file (directory_registry file) to appropriate hardware (by machine number)
rootDir = "Y:/"
# location of code repository
codeDir = "C:/Users/lkaiser/Dropbox/GitHub/FB_V2/"  
# project directory to save all outputs to
resultsDir = paste0(rootDir, "PICCC_analysis/FB_analysis/model_results/biomod2/") 
# location of necessary data to run scripts (maxent.jar, species csvs, crop rasters, etc.)
necessary_run_data = paste0(resultsDir,"necessary_run_data/") 
# location of bioclimativ variables data
bioclimDataDir = paste0(necessary_run_data,"clim_data/") 
# historical data at 125 m resolution
fitting_clim_data_dir = paste0(bioclimDataDir,"all_baseline/125m/") 
# current bioclim modeling variables at 500 m resolution
clim_data_2000=paste0(bioclimDataDir,"all_baseline/500m/")
# future bioclim modeling variables at 500 m resolution
clim_data_2100=paste0(bioclimDataDir,"all_future/500m/")

#######################################
##### GENERAL MODEL CONFIGURATION #####
#######################################

# assign project name to current run
project_name = "FB_ARCH_4runs" 
# choose species of interest
spp_nm = c('Akekee', 'Akikiki')
# NOTE: Due to data sensitivity, pseudo data is provided for Akekee and Akikiki as example

# pick BIOMOD2 models to run ('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT')
models_to_run = c("GBM","MAXENT")
# select evalutation methods ('KAPPA', 'TSS', 'ROC')
eval_stats = c("ROC")
# choose bioclimatic variables of interest 
env_var_files = c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
# to plot graph results (T) or not (F)
plot_graphs = FALSE
# to solve memory problems (T) for larger models or not (F)
apply_biomod2_fixes = TRUE
# to overwrite past results (T) or not (F)
overwrite = FALSE
# to turn on multi-instance excution for 1 species per thread (T) or not (F)
paralelize = FALSE
# set high number of cores to run threads on available computure
cpucores = 20 

### MAIN SCRIPTS ###
# to run the model fitting (T) or not (F)
EM_fit = TRUE     
# to run the ensemble modelling (T) or not (F)
EM_ensemble = TRUE  
# to project the model results (T) or not (F)
EM_project = TRUE   

### AUXILIARY SCRIPTS ###
# to get variable importance and evaluation score (T) or not (F) (script a)
merge_all_var_importance_and_model_eval = FALSE  
# to graph evaluation scores (T) or not(F) (script b)
model_eval_graph = FALSE
# to graphy variable importance (T) or not (F) (script c)
var_importance_graph = FALSE
# to create model response curves (T) or not (F) (script d)
create_response_curves = FALSE

### BASELINE vs FUTURE ###
# to map analog climates (T) or not (F)
create_analog_climate_map = TRUE
# to create raster files (T) or not (F)
raster_output_creation = TRUE 
# to calculate distribution shifts (T) or not (F)
distribution_shift_calculations = TRUE 
# to create ensemble maps (T) or not (F)
spp_ensemble_maps = TRUE

######################################
##### OPTIONS FOR SPECIFIC STEPS #####
######################################

### EM_fit (script #1)
# number of evaluation runs for ensemble modeling
NbRunEval = 4 #number of evaluation runs for ensemble modeling
# to include absences (T) or not (F)  **NOTE: in test phase
include_Abs = FALSE
# to consider PAs outside the climate envelope (T) or not (F) of collected points
PseudoAbs_outside_CE = FALSE
# If 1, then PA density will equal point density of surveyed areas
dens_PAs_outside_CE = 1 
# number of repetitions for PA selection
PA.nb.rep = 4
# number of PAs to select (overridden if if PseudoAbs_outside_CE = T)
PA.nb.absences = 1000
# number of candidate PAs to select (only if PAs_outside_CE = F, otherwise use PA.nb.absences)
candidatePAperPA = 200  
# strategy for selecting PAs ('random', 'sre', 'disk' or 'user.defined')
PA.strategy = "random" 
# distance equivalence ratio
equiv_100m = 0.0009430131
# minimum distance (500 m) of PAs from actual data points 
PA.dist.min = 5*equiv_100m
# to run the full models (T) or not (F)
do.full.models = TRUE

### EM_ensemble (script #2)
# set minimum evaluation score below which models will be excluding when building ensembles
eval.metric.threshold = rep(0.5,length(eval_stats)) 

### EM_project (script #3)
# project by island to avoid memory issues (T) or not (F)
proj_by_island = FALSE
# choose to run baseline (1) or future (4) projections
baseline_or_future = 4
# to save clamping mask (T) or not (F)
clampingMask = FALSE
# to store mask in memroy (T) or not (F) (if ClampingMask = T, mask saved to hard drive)
memory = TRUE

### create_analog_climate_map (script e)
# to compare with current (1) of future (4) climate
toCompareWithCurrentClimate = 4

### raster_output_creation (script f)
# to select statistic type of ensemble models
spp_ensemble_type = "wmean" 
# to select evaluation statistic of ensemble models
spp_ensemble_eval_stats = c("ROC") 
# projects to compare
comp_projects = c('baseline', 'future')
# to plot species ensembles (T) or not
plot_spp_ensemble_CV = TRUE
# to create mask map for species ensembles (T) or not (F)
masked_spp_ensemble_map = FALSE

### distribution_shift_calculations (script g)
# set model resolution (km)
model_resolution = 0.5
# to exclude areas outside habitat (T) or not (F)
exclude_areas_beyond_primary_habitat = FALSE

### spp_ensemble_maps (script h)
# to overlay habitat range (T) or not (F)
habitat_overlay = FALSE 
# to add prehistoric distribution from landfire BPS (T) or not (F)
BPS = FALSE 

###########################
##### RUNNING SCRIPTS #####
###########################

# assign working directory for project
working_dir = paste0(resultsDir, project_name, "/")
# create working directory if it does not already exist
dir.create(working_dir, showWarnings = FALSE)
# set working directory
setwd(working_dir)

# assign directory for cropped raster files
crop_raster_dir = paste0(working_dir, "map_crop/")
# assign directory fo single species csv files
csv_dir = paste0(necessary_run_data,"single_sp_CSVs/")
# create directory for output rasters if it does not already exist
dir.create(paste0(working_dir,"output_rasters/"), showWarnings = FALSE)

# set plot options depending on saving to server (T) or not (F)
useRasterDef = TRUE
interpolateDef = FALSE

# create directory for temporary data to avoid memory errors
dir_for_temp_files<-paste0(rootDir,'temp/', project_name, "/", baseline_or_future, "/")

# apply BIOMOD2 fixes to run if TRUE to solve memory problems (Line 53)
if (apply_biomod2_fixes){
  # assign directory for temporary maxent data for scenario
  maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
  # create directory for temporary maxent data
  dir.create(dir_for_temp_files, showWarnings = FALSE, recursive = TRUE)
}

# assign projected climate data set to use depending on scenario (Line 125)
if (baseline_or_future == 1){
  clim_data_dir = clim_data_2000 
  proj_nm = 'baseline'}
if (baseline_or_future == 2){
  clim_data_dir = clim_data_2000wettest
  proj_nm ='baseline_wettest'}
if (baseline_or_future == 3){
  clim_data_dir = clim_data_2000driest
  proj_nm = 'baseline_driest'}
if (baseline_or_future == 4){
  clim_data_dir = clim_data_2100 
  proj_nm = 'future'}
if (baseline_or_future == 5){
  clim_data_dir = clim_data_2100wettest
  proj_nm = 'future_wettest'}
if (baseline_or_future == 6){
  clim_data_dir = clim_data_2100driest
  proj_nm = 'future_driest'}
if (baseline_or_future == 7){
  clim_data_dir = clim_data_2100rev
  proj_nm = 'future_rev'}

# run scripts for model fitting, ensemble modeling, and model projections with settings above

# start the clock to time processing
ptmOverallStart<-proc.time()

if (EM_fit){ # 1 - run fitting code
  source(paste0(codeDir,"1_BM2_FB_SDM_fitting_LK.r")) 
}
if (EM_ensemble){ # 2 - run ensemble code
  source(paste0(codeDir,"2_BM2_FB_SDM_EM_creation_LK.r")) 
}
if (EM_project){ # 3 - run projection code
  source(paste0(codeDir,"3_BM2_FB_SDM_EM_projection_byIsland_LK.r"))   
}

#auxiliary scripts (a - h)
if (merge_all_var_importance_and_model_eval){ # a - variable importance/evaluation statistics
  source(paste0(codeDir,"1opt_merge_all_var_importance_and_model_eval.R"))}
if (model_eval_graph){ # b - evaluatation statiscs graphs
  source(paste0(codeDir,"1opt_model_eval_graph.R"))}
if (var_importance_graph){ # c - importance variable graphs
  source(paste0(codeDir,"1opt_var_importance_graph.R"))}
if (create_response_curves){ # d - response curves
  source(paste0(codeDir,"2opt_BM2_FB_SDM_response_curves.r"))}
if (create_analog_climate_map){ # e - map analog climates
  source(paste0(codeDir,"4opt_create_analog_climate_map_LK.R"))}
if (raster_output_creation){ # f - create output rasters
  source(paste0(codeDir,"4_SDM_raster_output_creation_newBiomVer.R"))}
if (distribution_shift_calculations){ # g - calcuate distribution shift
  source(paste0(codeDir,"6_distribution_shift_calculations_newbiomod2.r"))}
if (spp_ensemble_maps){ # h - map species ensembles
  source(paste0(codeDir,"7_spp_ensemble_maps.r"))}

# stop the clock and calculate processing time
ptmOverallElaps = proc.time() - ptmOverallStart
# create temporary variable of numeric value of the time elapsed per species
jnk = as.numeric(ptmOverallElaps[3])/length(spp_nm) 
# convert elapsed time in seconds to minutes
jnk = jnk/60 
# sign-posting of amount of time it took to run all codes per species
cat('\n It took ', jnk, 'minutes (on average) to model each species with',
    length(models_to_run), 'model types')

#######################
### END CONFIG FILE ###
#######################