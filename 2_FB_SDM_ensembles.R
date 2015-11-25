### BIOMOD2: ensemble models ###
### see source script inputs ###
### -step 2 = sdm ensembles- ### 

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)

# increase designated memory limit size
# memory.limit(size = 4095) 

# reset sp_nm counter to first species
sp_nm = spp_nm[2] 

# initialize parallel function
sp_parallel_run = function(sp_nm){  
  # load packages in snowfall
  library(biomod2)
  library(stringr)
  
  # definte species name as a character string
  sp_nm = as.character(sp_nm) 
  # replace "_" with "." in sp_nm
  sp_dir = str_replace_all(sp_nm,"_", ".") 
  # sink console outputs to text log file
  sink(file(paste0(working_dir,sp_dir, "/", sp_dir,Sys.Date(),
                   "_EM_creation_log.txt"), open = "wt"))
  
  # sign-posting noting the beginning of the ensemble creation
  cat('\n',sp_nm,'ensemble creation...')
  
  # set name of file to load workspace data from model fitting
  workspace_name = paste0(sp_nm,"_FB_modelfitting.RData") 
  # set name of file to save workspace data from ensemble models
  workspace_name_out = paste0(sp_nm,"_FB_EM_fit.RData") 
  
  # check to see if workspace output has already been done or not
  if (file.exists(workspace_name_out) == F | overwrite == T){ 
    # start the clock for processing time
    ptmModule2Start <- proc.time()
    
    # load model fitting results for species of interest
    load(workspace_name) 
    # replace "_" with "." in sp_nm
    sp_nm = str_replace_all(sp_nm,"_", ".") 
    # read data frame from species info csv file
    spp_info = read.csv(paste0(csv_dir, "FB_spp_data.csv")) 
    
    ############################
    ##### modeling_summary #####
    ############################
    # return summary of BIOMOD2 model fitting results for species
    myBiomodModelOut 
    # showClass("myBiomodModelOut")
    
    ############################
    ##### model_evaluation #####
    ############################
    # get model evaluation as array with results for each species model
    myBiomodModelEval<-get_evaluations(myBiomodModelOut) 
    
    ###############################
    ##### variable_importance #####
    ###############################
    # return an array with model variable importances (i.e bio1, bio7, etc)
    get_variables_importance(myBiomodModelOut) 
    
    #############################
    ##### remove_bad_models #####
    #############################
    # myBiomodModelEval_sum_cutoffs=apply(myBiomodModelEval, c(2,3,4,5), sum)
    # myBiomodModelEval_sum_cutoffs=c(myBiomodModelEval_sum_cutoffs[2, , , ])
    # good_cutoffs = which(!is.na(myBiomodModelEval_sum_cutoffs))
    # models.computed=myBiomodModelOut@models.computed
    # remaining_models=models.computed[good_cutoffs]
    
    # models computed to use from BIOMOD2 model object output
    remaining_models = myBiomodModelOut@models.computed
    
    #############################
    ##### ensemble_modeling #####
    #############################
    # combine models and make ensemble predictions built with BIOMOD_Modeling in module 1
    myBiomodEM<-BIOMOD_EnsembleModeling( 
      modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
      chosen.models = remaining_models,  #vector of model runs to use
      em.by = 'all',  #how models will be combined
      eval.metric = eval_stats, #evaluation metrics to build ensemble
      eval.metric.quality.threshold = eval.metric.threshold,  #threshold to exclude models
      prob.mean = TRUE,  #estimate mean probabilities 
      prob.cv = TRUE,  #estimate coefficient of variation
      prob.ci = TRUE,  #estimate confidence interval of prob.mean
      prob.ci.alpha = 0.05,  #signficance level for estimating confidence interval
      prob.median = TRUE,  #estimate median
      committee.averaging = TRUE,  #estimate committee averaging
      prob.mean.weight = TRUE,  #estimate weighted sums
      prob.mean.weight.decay = 'proportional' ) #define relative importance of weights
    
    # sign-posting indicating that ensemble modeling is done
    cat('\n',sp_nm,'ensemble done.')

    #####################################
    ##### ensemble_modeling_outputs #####
    #####################################
    # return summary from ensemble modeling
    myBiomodEM
    
    # return evaluation scores(testing.data, cutoff, sensitivity, & specificity)
    get_evaluations(myBiomodEM)
    # save workspace of R environment 
    save("myBiomodEM", "myBiomodModelOut","remaining_models", file = workspace_name_out)
    
    # stop the clock and calculate processing time to run the code
    ptmModule2Elaps = proc.time() - ptmModule2Start 
    # assign temporary numeric value to the time elapsed
    jnk = as.numeric(ptmModule2Elaps[3]) 
    # convert time elapsed from seconds to minutes
    jnk = jnk/60
    # sign-posting of time it took to build ensemble model per species
    cat('\n It took', jnk, 'minutes to ensemble model', sp_nm)
    
    # otherwise
    }else{
      # sign-posting if ensemble modeling has previously been done for this species run
      cat('\n',sp_nm,'ensemble previously done.')
  }
  # return output to console
  sink(NULL)
}
# END snowfall function

########################################
##### SNOWFALL FUNCTION SCRIPT RUN #####
########################################

if (is.null(cpucores)){
  # determine total number of processors available
  cpucores = as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))  
}else{
  # determine minimum number of processors available
  cpucores = min(cpucores, as.integer(Sys.getenv('NUMBER_OF_PROCESSORS')))
}
# initialize parallel computing on minimum number of cpu cores
sfInit(parallel = TRUE, cpus = cpucores)
# export all environmentals variables to cores for cluster programming
sfExportAll() 
# time parallel run calculation of function across cores
system.time((sfLapply(spp_nm, fun = sp_parallel_run)))
# remove all global environmental variables from cluster cores
sfRemoveAll()
# end parallel computing on other cpu cores
sfStop()

#################################
##### END ENSEMBLE MODELING #####
#################################