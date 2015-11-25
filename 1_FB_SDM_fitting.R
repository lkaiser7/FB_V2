### BIOMOD2: fit sdm to data ###
### see source script inputs ###
### -step 1 = model fitting- ###

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)
library(tools)

# copy the necessary data to run the models into the working directory
dirs = list.dirs(gsub(".$", "", necessary_run_data), full.names = TRUE, recursive = TRUE)
for (dir in dirs){
  layers <- list.files(dir, pattern = NULL, full.names = FALSE, include.dirs = FALSE)
  for (layer in layers){
    layer_full_nm = paste(dir,layer, sep = "/")
    if (file.info(layer_full_nm)$isdir == FALSE){
      out_dir_nm = str_replace(dir, necessary_run_data, working_dir)
      dir.create(out_dir_nm, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      out_lyr_nm = str_replace(layer_full_nm, necessary_run_data, working_dir)
      if (file.exists(out_lyr_nm) == FALSE){
        cat('\n','found ', layer, 'in ', dir)
        file.copy(layer_full_nm, out_lyr_nm, overwrite = TRUE, recursive = FALSE,
                  copy.mode = TRUE)  # if'recursive' = T will result in warnings
        cat('\n','saved as ', out_lyr_nm)
      }
    }
  }
}
# sign-posting noting all necessary data has been copied to working directory
cat('\n Copying of necessary files is complete.')

# assign location of species information csv (e.g. raster directory, raster set, endemic...)
spp_info = read.csv(paste0(csv_dir, "FB_spp_data.csv"))

# create vector of bioclimatic variable names without the file extension (.grd)
var_names <- unlist(file_path_sans_ext(env_var_files))

# set sp_nm = spp_nm[1] for debugging
n_abs_removed = c()

# initialize snowfall parallel computing function
sp_parallel_run = function(sp_nm){
  # load required libraries
  library(biomod2)
  library(raster)
  library(randomForest)
  library(dismo)
  library(mda)
  library(stringr)
  library(tools)
  require(snowfall)
  
  # convert sp_nm to character in case species are numbered
  sp_nm = as.character(sp_nm) 
  # replace "_" with "." in sp_nm
  sp_dir = str_replace_all(sp_nm, "_", ".") 
  # create new directory per species
  dir.create(sp_dir, showWarnings = FALSE) 
  # write console outputs to text file to log processing
  sink(file(paste0(working_dir, sp_dir, "/", sp_dir, Sys.Date(),
                   "_fitting_log.txt"), open = "wt"))
  
  # sign-posting of start of model fitting per species
  cat('\n', sp_nm, 'model fitting...') 
  # save file name of variable importance
  FileName00<-paste0(sp_nm, "_VariImp.csv") 
  # check if analysis for species is already done or not
  if (file.exists(FileName00) == FALSE | overwrite == TRUE){  #allow for overwrite capacity    
    # start the clock to calculate processing time
    ptmModule1Start<-proc.time()
    # save file name for R workspace data after mode fitting
    workspace_name = paste0(sp_nm, "_FB_modelfitting.RData") 
    
    #####################
    ##### LOAD DATA #####
    #####################
    
    # sign-posting to laod raster_based_env_grid
    cat('\n loading rasters...') 
    
    # match species index in csv file to species being processed
    sp_index = which(spp_info[,"Species"] == sp_nm) 
    # locate raster directory for species being processed in csv file
    raster_res = paste0("/", spp_info[sp_index, "rasterdir"]) 
    # assign raster to crop bioclimatic variables by per species
    crop_raster = raster(paste0(crop_raster_dir, raster_res, ".grd")) 
    # assign bioclimatic variable mask as predictors object
    predictors = raster(paste0(fitting_clim_data_dir, env_var_files[1])) 
    # crop predictor grid mask to species raster
    predictors = crop(predictors, crop_raster)  
    # add additional bioclimatic variables to predictors raster object
    for (jj in 2:length(env_var_files)){ 
      # create temporary raster of each additional bioclimatic variable
      temp = raster(paste0(fitting_clim_data_dir, "/", env_var_files[jj]))
      # crop each additional bioclimatic variable to species raster
      temp = crop(temp, crop_raster)
      # add each additional bioclimatic variable to predictors raster stack
      predictors = addLayer(predictors, temp)
    }
    # assign names of bioclimatic variables to raster stack
    names(predictors)<-var_names
    # remove temporary variables
    rm("jj", "crop_raster", "temp")
    
    # save name of jpeg file to be created for selected environmental variables
    jpeg_name = paste0(sp_nm, "_env_vars_used.jpg")
    # create blank jpeg file in working directory
    jpeg(jpeg_name, width = 10, height = 10, units = "in", 
         pointsize = 12, quality = 90, bg = "white", res = 300)
    # plot predictors raster stack of bioclimatic variables 
    plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, 
         useRaster = useRasterDef, axes = TRUE, addfun = NULL)  
    # save jpeg file
    dev.off()
    
    # sign-posting loading species point data
    cat('\n loading species data...')
    # load presence and absence data points csv file per species
    mySpeciesOcc = read.csv(paste0(csv_dir, "/", sp_nm, '_pres_abs.csv')) 
    
    # extract X (lon), Y (lat) and presence (and absence) data from species csv
    mySpeciesOcc = cbind(mySpeciesOcc[, 2:3], pa = mySpeciesOcc[, 1]) 
    # check true occurrence presence/absence data
    head(mySpeciesOcc)
    # check if true absences should be included (T) or not (F)
    if (!include_Abs){
      # select only true presence data points if absences are excluded
      mySpeciesOcc = mySpeciesOcc[mySpeciesOcc$pa == 1, ]
    }
    
    # sign posting for pseudo-absence handling to define points
    cat('\n defining candidate PA points...')
    
    # extract only lat (Y) and long (X) points from occurrence data
    PA_XY = mySpeciesOcc[,1:2] 
    # build raster layer of response variable from environmental mask
    mySREresp<-reclassify(subset(predictors, 1, drop = TRUE), c(-Inf, Inf, 0)) 
    # assign shared cells from bioclimatic variables and occurrences to '1'
    mySREresp[cellFromXY(mySREresp, PA_XY)]<-1 
    # calculate the number of true absences
    Act_abs = dim(mySpeciesOcc[mySpeciesOcc$pa == 0,])[1] 

    # loop to assign pseudo absences outside climate envelope (T) or randomly (F)
    if (PseudoAbs_outside_CE){
      # set tail threshold for surface range envelop values extremes
      sre_tail = 0.025
      sp_CE = sre(Response = mySREresp, Explanatory = predictors, 
      # calculate SRE by removing extreme (tail) values and correct point density in SRE
                  NewData = predictors, Quant = sre_tail) 
      # count the number of cells with species occurrence data
      n_PandA = sum(as.matrix(mySREresp), na.rm = TRUE)*(1-sre_tail*2) 
      # count the number of cells with presence or absences predicted in climate envelope
      CE_cells = sum(as.matrix(sp_CE), na.rm = T) 
      # determine the point density of true data within the climate envelope
      CE_point_density = round(n_PandA/CE_cells, digits = 4) 
      # create raster of all cells outside the climate envelope
      neg_sp_CE = sp_CE == 0 
      # calculate the number of cells outside the climate envelope 
      neg_CE_cells = sum(as.matrix(neg_sp_CE), na.rm = T) 
      # account for actual absences in true data
      n_PseudoAbs_pts = round(neg_CE_cells*dens_PAs_outside_CE*CE_point_density) + Act_abs 
      # create matrix of potential pseudo absence candidate lon/lat (X, Y) points      
      PseudoAbs_cand_pts = rasterToPoints(neg_sp_CE, fun = function(x){x == 1}) 
      
      # assign pseudo absences anywhere without limitation to climate envelop     
    }else{
      # create raster of areas outside the true occurrence data
      neg_mySREresp = mySREresp == 0 
      # plot raster of areas outside the true known data
      # plot(neg_mySREresp)
      # create matrix of potential pseudo absence candidate lon/lat (X, Y) points      
      PseudoAbs_cand_pts = rasterToPoints(neg_mySREresp, fun = function(x){x==1}) 
      
      # determine the number of pseudo absence points to select
      if (candidatePAperPA == 0){
        # assign number of pseudo absences to 0_config_file inputs
        n_PseudoAbs_pts = PA.nb.absences  
        
      # determine the number of pseudo absence points to select
      }else{
        # assign number of pseduo absences based onn normalized data dimensions 
        n_PseudoAbs_pts = round(dim(PseudoAbs_cand_pts)[1]/candidatePAperPA)        
      }
    }  
    
    # extract only X (lon) and Y (lat) data for candidate pseudo absence points
    PseudoAbs_cand_pts = as.data.frame(PseudoAbs_cand_pts[,1:2]) 
    # check pseudo absence point data    
    head(PseudoAbs_cand_pts) 
    # return dimensions (points, rows) of pseudo absence point data
    dim(PseudoAbs_cand_pts) 
    # create data frame excluding any missing geographic location information
    PseudoAbs_cand_pts_noNA = PseudoAbs_cand_pts[complete.cases(PseudoAbs_cand_pts),] 
    # add 'pa' column of "NA" to designate pseudo absence poitns
    PseudoAbs_cand_pts_noNA = cbind(PseudoAbs_cand_pts_noNA,
                                    pa = rep('NA', dim(PseudoAbs_cand_pts_noNA)[1],1)) 
    # rename column headers 
    names(PseudoAbs_cand_pts_noNA) = c('X', 'Y', 'pa') 
    # check pseudo absences with all columns and no missing data
    head(PseudoAbs_cand_pts_noNA) 
    
    # merge true occurrence data with pseudo absences in new data frame
    mySpeciesOcc_w_Pseudo<-data.frame(rbind(mySpeciesOcc, PseudoAbs_cand_pts_noNA)) 

    # sign posting of extraction of environmental variables to point data
    cat('\n extracting env vars to points...')
    
    # create matrix with select bioclimatic variables and cell numbers of all data points
    relBioclimData<-extract(predictors, mySpeciesOcc_w_Pseudo[,1:2], cellnumbers = TRUE) 
    # create data frame with bioclimatic data for all data points
    XY_PresAbsPA_Bioclim = data.frame(cbind(mySpeciesOcc_w_Pseudo, relBioclimData)) 
    # check bioclimatic point data
    head(XY_PresAbsPA_Bioclim) 
    # remove any rows with missing data
    XY_PresAbsPA_Bioclim_noNA = XY_PresAbsPA_Bioclim[complete.cases(
      XY_PresAbsPA_Bioclim[4:dim(XY_PresAbsPA_Bioclim)[2]]),] 
    # check bioclimatic data points with no missing data
    head(XY_PresAbsPA_Bioclim_noNA) 
    
    # create temporary vector of 'pa' column to select, count, and remove duplicate cells
    jnk = c(1, 0, NA)
    # sort data to remove pseudo absences before true absence and presence points
    XY_PresAbsPA_Bioclim_sort<-XY_PresAbsPA_Bioclim_noNA[order(match(
      XY_PresAbsPA_Bioclim_noNA$pa, jnk)),] 
    # remove temporary junk vector
    rm(jnk)    
    # identify duplicate cells in 'pa' column 
    dupEntries<-duplicated(XY_PresAbsPA_Bioclim_sort$cells) 
    # count number of duplicate entries
    n_dups = length(dupEntries[dupEntries == TRUE]) 
    # sign-posting detailing number of entries removed per species 
    cat('\n Out of ', length(dupEntries), 'points,', n_dups, 
        'were removed because they were within the same raster cell for', sp_nm)
    # create data frame with duplicates removed
    PresAbsPA_noDup <- XY_PresAbsPA_Bioclim_sort[!dupEntries, -4] 
    # check all data with duplicates removed
    head(PresAbsPA_noDup)
    # save output table (too large!)
    # write.table(PresAbsPA_noDup, file=paste0(sp_nm,"_loc_data_table.csv"), col.names=NA)
    
    # table number of values for presences, absences, and pseudo absences
    sp_loc_summary = table(PresAbsPA_noDup$pa, useNA = "ifany")
    # store table as a data frame
    sp_loc_summary = as.data.frame(sp_loc_summary)
    # rename levels of 'pa' attributes
    levels(sp_loc_summary$Var1) = c(0, 1, NA, "n to select")
    # create temporary object of pseudo absences to be selected per run
    jnk = c("n to select", n_PseudoAbs_pts)
    # combine number of pseudo absences to select with table values
    sp_loc_summary = rbind(sp_loc_summary, jnk)
    # rename column headers
    names(sp_loc_summary) = c("Data_type", "Freq")
    # save output table 
    write.table(sp_loc_summary, file=paste0(sp_nm, "_loc_data_summary.csv"), col.names=NA)

    ##########################
    ##### PLOT DATA MAPS #####
    ##########################
    
    # save name of jpeg file to be create per species
    jpeg_name2 = paste0(sp_nm, "_loc_data_used.jpg")
    # create blank jpeg file in working directory
    jpeg(jpeg_name2, width = 10, height = 10, units = "in", 
         pointsize = 12, quality = 90, bg = "white", res = 300)
    # plot map with set coordinate system per species location data 
    plot(seq((min(PresAbsPA_noDup[,1])-0.1), (max(PresAbsPA_noDup[,1])+0.1), 
             by = ((max(PresAbsPA_noDup[,1])+0.1) - (min(PresAbsPA_noDup[,1])-0.1))/5), 
         seq((min(PresAbsPA_noDup[,2])-0.1), (max(PresAbsPA_noDup[,2])+0.1), 
             by = ((max(PresAbsPA_noDup[,2])+0.1) - (min(PresAbsPA_noDup[,2])-0.1))/5), 
         type = "n", xlab = "Lon", ylab = "Lat")
    # add all data points
    points(x = PresAbsPA_noDup[is.na(PresAbsPA_noDup[,3]), 1], 
           y = PresAbsPA_noDup[is.na(PresAbsPA_noDup[,3]), 2], 
           type = "p", col = "grey", pch = 20, cex = 0.7)
    # add only true absence (0) points
    points(x = PresAbsPA_noDup[PresAbsPA_noDup[,3] == 0, 1], 
           y = PresAbsPA_noDup[PresAbsPA_noDup[,3] == 0, 2], 
           type = "p", col = "red", pch = 20,cex = 0.7)
    # add only true presence (1) points
    points(x = PresAbsPA_noDup[PresAbsPA_noDup[,3] == 1, 1], 
           y = PresAbsPA_noDup[PresAbsPA_noDup[,3] == 1, 2], 
           type = "p", col = "blue", pch = 20,cex = 0.7)
    # save jpeg file
    dev.off()
    
    #############################################
    ##### BIOMOD DATA FORMATING AND FITTING #####
    #############################################
    
    # sign-posting defining variables used by BIOMOD2 package
    cat('\n biomod model config...')
    
    # new data frame of XY coordinates only
    myRespXY = PresAbsPA_noDup[,1:2] 
    # new data frame for presences (1), absences (0), or PAs (NA)
    myResp <- data.frame(Sp_Bio = PresAbsPA_noDup[,3]) 
    # re-assign all 'NA' pseudo-absences to real <NA> for formating
    myResp[myResp =='NA'] = NA 
    # new data frame with bioclimatic point data only
    myExpl = PresAbsPA_noDup[,4:dim(PresAbsPA_noDup)[2]]

    # load BIOMOD2 data for formatting
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = myResp,  #pres/abs/pa points
      expl.var = myExpl,  #bioclim values
      resp.xy = myRespXY, #xy coordinates
      resp.name = sp_nm,  #species name
      PA.nb.rep = PA.nb.rep,  #number of PAs selections
      PA.nb.absences = n_PseudoAbs_pts,  #number of PAs to select
      PA.strategy = PA.strategy,  #how to select PAs
      PA.dist.min = PA.dist.min)  #minimum distance to presences

    # sign-posting of completed biomod data formating
    cat('\n biomod formatting data complete')
    
    # plot data with selected pseudo absences (takes too long!)
    if (plot_graphs == TRUE & PA.nb.rep < 9){
      # assign location for jpeg file per species
      jpeg_name3 = paste0(sp_nm, "_loc_data_used2.jpg") 
      # create blank jpeg file in assigned location
      jpeg(jpeg_name3, width = 10, height = 10, units = "in", 
           pointsize = 12, quality = 90, bg = "white", res = 300)
      # plot BIOMOD2 data used (true occurrences, true absences, and psuedo absences)
      plot(myBiomodData) 
      dev.off()
    }
    
    # increase the limit for total allocation of memory
    # memory.limit(size = 24000000)

    # set different options for selected modeling techniques from source script
    myBiomodOption <- BIOMOD_ModelingOptions(
      # GBM: Generalized Boosted Regression
      GBM = list( distribution = 'bernoulli', interaction.depth = 7,  shrinkage = 0.001, 
                  bag.fraction = 0.5, train.fraction = 1, n.trees = 100, cv.folds = 10),
      # MARS: Multivariate Adaptive Regression Splines
      MARS = list( degree = 2, penalty = 2, thresh = 0.001, prune = TRUE),
      # RF: Random Forest Classification and Regression
      RF = list(do.classif = TRUE, ntree = 100, mtry = 'default', 
                max.nodes = 10, corr.bias = T), 
      # MAXENT: Maximum Entropy
      MAXENT = list(maximumiterations = 100, visible = F, linear = TRUE, 
                    quadratic = TRUE, product = TRUE, threshold = TRUE, hinge = TRUE, 
                    lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, 
                    beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
                    beta_hinge = -1, defaultprevalence = 0.5))
    
    # sign-posting for ensemble model fitting
    cat('\n model fitting...')
    
    # create ensemble model from formatting data and modeling options
    myBiomodModelOut <- BIOMOD_Modeling(
      myBiomodData,  #formatted biomod data  
      models = models_to_run,  #select models to run          
      models.options = myBiomodOption,  #biomod options object                
      NbRunEval = NbRunEval,  #number of evaluation runs      
      DataSplit = 80,  #amount of data to use for training
      Yweights = NULL,  #response points weights                  
      VarImport = 10,  #permuations to estimate variable importance
      do.full.models = do.full.models,  #calibrate and evalute to all data
      models.eval.meth = eval_stats,  #evaluation metrics                            
      SaveObj = TRUE,  #save output object
      rescal.all.models = TRUE)  # scale to binomial GLM
    
    # return summary output of biomod model runs
    myBiomodModelOut
    
    # return output model evaluation metrics results
    myBiomodModelEval<-get_evaluations(myBiomodModelOut)
    # return names of model evaluations
    dimnames(myBiomodModelEval)
    
    # Validate selectec metrics for all tests (TSS, ROC, or KAPPA)
    if ("TSS" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["TSS", "Testing.data",,,] 
      # create data frame with variable importance values
      Spp_TSS<-data.frame(myBiomodModelEval["TSS", "Testing.data",,,]) 
      # assign file path name for results 
      FileName<-paste0(sp_nm, "_TSS.csv") 
      # create .csv file and save TSS outputs
      write.table(Spp_TSS, file = FileName, sep = ",", col.names = NA) 
    }  
    if ("ROC" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["ROC", "Testing.data",,,]
      # create data frame with variable importance values
      Spp_ROC<-data.frame(myBiomodModelEval["ROC", "Testing.data",,,])
      # assign file path name for results 
      FileName<-paste0(sp_nm, "_ROC.csv")
      # create .csv file and save ROC outputs
      write.table(Spp_ROC, file = FileName, sep = ",", col.names = NA)
    }
    if ("KAPPA" %in% eval_stats){
      # return variable importance for each model
      myBiomodModelEval["KAPPA", "Testing.data",,,]
      # create data frame with variable importance values
      Spp_KAP<-data.frame(myBiomodModelEval["KAPPA", "Testing.data",,,])
      # assign file path name for results 
      FileName<-paste0(sp_nm, "_KAP.csv")
      # create .csv file and save KAPPA outputs
      write.table(Spp_KAP, file = FileName, sep = ",", col.names = NA)
    }
    
    # get the variable importance for selected bioclimatic variables
    get_variables_importance(myBiomodModelOut) 
    # create data frame with model variable importances
    Spp_VariImp<-data.frame(get_variables_importance(myBiomodModelOut)) 
    # save variable importance as csv file
    write.table(Spp_VariImp, file = FileName00, sep = ",", col.names = NA) 
    
    # save model fitting R workspace
    save("myBiomodModelOut", file = workspace_name)
    
    # calculate total processing time to run code 
    ptmModule1Elaps = proc.time() - ptmModule1Start 
    # create temporary object to store numeric value of time elapsed
    jnk = as.numeric(ptmModule1Elaps[3]) 
    # convert time from seconds to minutes
    jnk = jnk/60
    # sign-posting of processing time per species 
    cat('\n','It took ', jnk, "minutes to model", sp_nm)
  }else{
    # sign-posting if file for variable importance has already been created 
    cat('\n fitting for', sp_nm, 'already done...') 
    # indicates that this species has already been run 
  }    
  # reset sink to console outpu
  sink(NULL)
} # END snowfall function

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

#############################
##### END MODEL FITTING #####
#############################