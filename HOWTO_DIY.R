### HOW TO DO IT YOURSELF ### 

# a basic guide to using the FB_V2 GitHub repository #
# steps to build and run SDMs in the BIOMOD2 package #
# directions for using the 0_config_file and scripts #

# Welcome to the GitHub repository for native Hawaiian Forest Bird Species Distribution Models!
# This archive of code contains the scripts necessary to run comparisons of current 
# and future climate scenarios for mulitple species. Here are basic directions how 
# to use this repository of scripts and build your own species distribution models.

###################
# LIST OF SCRIPTS #
###################

### MAIN SCRIPTS ###

# 0_config_file_STARTHERE.R
# 1_FB_SDM_fitting.R
# 2_FB_SDM_ensembles.R
# 3_FB_SDM_projections.R

### AUXILIARY SCRIPTS ###

# 1a_merge_var_imp_and_model_eval.R
# 1b_model_eval_graph.R
# 1c_var_imp_graph.R
# 2a_response_curves.R
# 3a_projection_mods.R
# 4_raster_output_creation.R
# 5_analog_climate_map.R
# 6_distribution_shift_calcs.R
# 7_ensemble_maps.R

##############
# DIRECTIONS #
##############

# The main configuration file (0_config_file) holds all the inputs for the other scripts.
# Any of these settings and configurations can be changed in that main script to reflect
# a certain desired output.  There are a few main lines that should generally be changed.
# Setting source locations of all data and scripts will vary by each individual computer.
# These basic directions use the example forest bird data provided for Akekee and Akikiki
# to compare baseline and future projected ranges of the distributions for these species.

# 1. assign your project run a name (line 40)
# this will be unique to the individual analysis currently running

# 2. change any basic model configurations you see fit (lines 41-60)
# here you can select certain inputs based on your desired outcome

# 3. decide which scripts to run (lines 64-91)
# set script names to either run (TRUE) or not (FALSE)
# set lines 75-81 to TRUE for baseline and FALSE for future
# set lines 85-91 to FALSE for baseline and TRUE for future

# 4. select number of runs/iterations (line 99)
# this will be done per each model per species

# 5. select which climate scenario to run (line 129)
# choose baseline (1) or future (4)

# 6. run the entire script once with the baseline settings (steps 3 & 5)

# 7. change to future settings and rerun the entire script (steps 3 & 5)

# By running the main configuration file (0_config_file) twice, once for baseline and 
# once for future, the main and auxiliary scripts should create two climate scenarios 
# to compare the distribution and range changes of the selected species.  From here, 
# an additional niche overlap analysis can be applied once this modeling is complete.

### NICHE OVERLAP SCRIPTS ###

#
#

##############
# DIRECTIONS #
##############

# 8.
