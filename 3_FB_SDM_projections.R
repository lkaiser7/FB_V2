### BIOMOD2: sdm projections ###
### see source script inputs ###
### -step 3 = project sdms - ###

##################
##### SET UP #####
##################

# load necessary packages
require(snowfall)

# create data frame from species csv info file
spp_info = read.csv(paste0(csv_dir,'FB_spp_data_HI.csv')) 

# reset sp_nm counter to first species
sp_nm = spp_nm[2] 

# initialize parallel function
sp_parallel_run = function(sp_nm){
  # load packaged in snowfall
  library(biomod2)
  library(stringr)
  library(colorRamps)
  library(rasterVis)
  library(tools)
  library(ncdf)
  library(gbm) 
  
  # set options for biomod2 fixes in code (if assigned TRUE in source script)
  if (apply_biomod2_fixes){
    rasterOptions(tmpdir = dir_for_temp_files, timer = T, progress = "text", todisk  = T) #set options for raster package
    source(paste0(codeDir,"3a_projection_mods.R")) #all of fixes to biomod2 code created by AV
  }
  
  # definte species name as a character string
  sp_nm = as.character(sp_nm) 
  # replace "_" with "." in sp_nm
  sp_dir = str_replace_all(sp_nm,"_", ".") 
  
  # sink console outputs to text log file
  sink(file(paste0(working_dir, sp_dir, "/", sp_dir,Sys.Date(),
                   "_proj_log.txt"), open="wt"))
  
  # sign-posting beginning model projections
  cat('\n',sp_nm,'model projection...') 
  # sign-posting loading previous workspace data from ensemble fitting
  cat('\n','loading EM_fit workspace file...')
  
  # set file name of workspace data to load from ensemble fitting
  workspace_name = paste0(sp_nm,"_FB_EM_fit.RData") 
  # load workspace from previous step
  load(workspace_name) 
  # assign location for all ensemble model plot outputs
  plots = paste0(working_dir, "AllEMplots_pmw") 

  # run if overwrite function in source script is off (F)  
  if (file.exists(plots) == FALSE | overwrite == TRUE) { 
    # create directory for plot outputs
    dir.create(plots, showWarnings = FALSE)
  }
  
  # set name of file to save workspace data from projections
  workspace_name_out = paste0(sp_nm,"_FB_EM_proj_", proj_nm, ".RData") 
  
  # check if analysis for species is already done or not
  if (file.exists(workspace_name_out) == F | overwrite == T){ #allow for overwrite capacity
    
    # sign-posting of loading bioclimatic variables data
    cat('\n','loading predictors...')
    
    # assign index for species from info file
    sp_index = which(spp_info[,"Species"] == sp_nm) 
    # assign raster resolution directory for species
    raster_res = spp_info[sp_index, "rasterdir"] 
    
    # sign-posting indicating which raster is used for projections
    cat('\n using these env files for projection raster:', env_var_files, 
        '\n from dir:', clim_data_dir)
    
    # create raster layer from resolution extent .grd file
    crop_raster = raster(paste0(crop_raster_dir,raster_res,".grd")) 
    # create raster layer from bioclimatic variables
    predictors = raster(paste0(clim_data_dir, env_var_files[1])) 
    # crop predictors to raster resolution extent
    predictors = crop(predictors,  crop_raster) 
    # assign temporary variable to count number of bioclimatic variables for runs
    jnk0 = length(env_var_files)
    
    # add all other bioclimatic variables raster layers to raster stack
    for (jj in 2:jnk0){
      # create temporary raster layer of each new bioclim variable layer to add
      temp = raster(paste0(clim_data_dir, env_var_files[jj]))
      # crop temporary raster layer to raster resolution extent
      temp = crop(temp,  crop_raster)
      # add new bioclimatic variable layer to raster stack
      predictors = addLayer(predictors, temp)
    }
    
    # create vector of bioclimatic variables names without file extensions
    var_names<-unlist(file_path_sans_ext(env_var_files)) 
    # assign names of raster stack layers to bioclimatic variable names
    names(predictors)<-var_names 
    # remove temporary variables
    rm("crop_raster" ,"temp", "jnk0") 
    # return summary of raster stack of predictors
    predictors 

    # sign-posting noting completion of loading predictors
    cat('\n done loading predictors.')   
    
    # to project island by island
    if (proj_by_island){
      # sign-posting to project by indivdual island
      cat('\n','creating predictors by island...')
      
      # Define extent of each different island
      Kauai = c(-159.82,-159.26, 21.84, 22.25)
      Oahu = c(-158.32, -157.62,  21.22, 21.73)
      Molokai = c(-157.34, -156.69, 21.03, 21.25)
      Lanai = c(-157.08, -156.78, 20.70, 20.92)
      Maui= c(-156.8, -155.53, 20.46, 21.05)
      Hawaii = c(-156.10,-154.74, 18.87, 20.30)
      Kahoolawe = c(-156.8, -156.51, 20.46, 20.62)
      
      # identify which islands the species is found on based on info file row number
      sp_row<-which(spp_info[, "Species"] == sp_nm) 
      # return data frame of islands species is located on if on more than one island
      spIslands<-spp_info[sp_row,(6:length(names(spp_info)))] 
      # list names of islands where species is found
      spIslandNames<-names(spIslands)[spIslands > 0] 
      
      # add Kahoolawe to species island list if Maui is included
      if ("Maui" %in% spIslandNames) {
        spIslandNames<-append(spIslandNames, "Kahoolawe")
      }
      
      # cut out each island separately
      for (i in 1:length(spIslandNames)){
        # select extent for specific island
        e = extent(get(spIslandNames[i]))
        # crop and stack predictors to island extent
        Isras = stack(crop(predictors, e, snap = 'in'))
        # rename bioclimatic variables in raster stack
        names(Isras)<-var_names
        # assign raster to list of island names where species is found
        assign(spIslandNames[i], Isras)
      }
      
      # remove Kahoolawe from Maui to avoid extent and reclassification issues
      if ("Maui" %in% spIslandNames) {
        # combine values into a vector for reclassification
        rcl<-c(0, 10000, -1)
        # create a matrix from set vector values
        rcl<-matrix(rcl, ncol = 3, byrow = TRUE)
        # define Kahoolawe raster values as -1 to separate from Maui once merged
        a = reclassify(Kahoolawe, rcl)  
        # merge Maui and reclassifed Kahoolawe data
        Maui = try(merge(a, Maui), TRUE)
        # combine values into a vector for reclassification with NAs
        rcl<-c(-1, NA)
        # create a matrix from set vector values with NAs
        rcl<-matrix(rcl, ncol = 2, byrow = TRUE)
        # merge Maui data with reclassified NA data
        Maui = try(reclassify(Maui, rcl), TRUE)   
        # create raster stack from new reclassified Maui data
        Maui = stack(Maui)
        # rename bioclimatic variables in raster stack
        names(Maui)<-var_names
        # remove Kahoolawe from list of islands where species is found
        spIslandNames <- spIslandNames[which(spIslandNames != "Kahoolawe")] 
      }
      
      # sign-posting of completion of predictors per island
      cat('\n done creating predictors by island.') 
      
    }else{
      # otherwise create blank list of island names for species
      spIslandNames = c("")
    }
    
    ######################################
    ##### RUN PROJECTIONS PER ISLAND #####
    ######################################
    
    # set islands per species counter to 1 for debugging
    spIsland = spIslandNames[1] 
    # loop through each island per species
    for (spIsland in spIslandNames){
      # if species exists on more than one island
      if (length(spIslandNames) > 1) {
        # project species per individual island
        projection_name = paste0(proj_nm, "_", spIsland)
        # rename species per island 
        sp_is = paste0(sp_nm, "_", spIsland)
      } else {
        # otherwise set projection name to scenario
        projection_name = proj_nm
        # keep species name per island
        sp_is = sp_nm
      }
      
      # set location to save R workspace data after model projections
      workspace_name_out0 = paste0(sp_is, "_FB_all_model_proj_", proj_nm, ".RData") 
      # check if file has already been created or not 
      if (file.exists(workspace_name_out0) == F | overwrite == T){ #allow to overwrite
        # list islands per species if projection per island is TRUE
        if (proj_by_island){predictors<-get(spIsland)}
                
        # change raster for the "wettest" and "driest" scenarios 
        alt_scen = c(2,3,5,6)
        # if cliate scenario uses an alternate wet or dry scenario
        if (baseline_or_future %in% alt_scen){
          # create new raster stack of predictors adjust to alternate scenario
          predictors<-stack((subset(predictors, 1)),
                            (subset(predictors, 2)),
                            (subset(predictors, 3))*4,
                            (subset(predictors, 4)))
          # rename preditors by bioclimatic variables
          names(predictors)<-var_name
        }
        
        # sign-posting on ongoing plotting of predictors
        cat('\n','plot predictors...')
        
        # for baseline (current) runs
        if (baseline_or_future == 1){ 
          # save name of jpeg file to be created for basline scenarios
          jpeg_name = paste0(sp_is, "_env_vars_used_for_projection_", projection_name, ".jpg") #assigns location for map jpeg file
          # create blank jpeg file in directory
          jpeg(jpeg_name, width = 10, height = 10, units = "in", 
               pointsize = 12, quality = 90, bg = "white", res = 300) 
          # plot predictors for baseline projections
          plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, axes = TRUE,
               useRaster = useRasterDef, addfun = NULL, interpolate = interpolateDef) 
          # save jpeg file
          dev.off()
        }
        
        # sign-posting of created raster stack projection
        cat('\n', sp_is,'projection raster stack created.')
        # reclaim memory no longer in use and return memory usage summary
        gc()
        
        # sign-posting to project BIOMOD2 outputs
        cat('\n','run biomod_projection...')
        
        # for baseline projections
        myBiomodProj_baseline<-BIOMOD_Projection(
          modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
          new.env = predictors,  #new explanatory variables to project onto for future
          proj.name = projection_name,  #projection directory
          selected.models = remaining_models,  #which models to use for projections
          binary.meth = eval_stats,  #evaluation method statistics 
          compress = 'xz',  #compression format of files on hard drive
          build.clamping.mask = clampingMask,  #if clamping mask should be saved or not
          keep.in.memory = memory) #if clamping mask should be saved to hard disk or not
        
        # reclaim memory no longer in use and return memory usage summary
        gc() #reclaims memory that is no longer used.
        # sign-posting of completed projections
        cat('\n',sp_nm,'projection complete.')
        # sign-posting of memory size and availability
        cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn') 
        
        # save projection R workspace environement 
        save("myBiomodProj_baseline", file = workspace_name_out0)
        # sign-posting of completed BIOMOD2 projections
        cat('\n done with running biomod_projection.')
        
        # otherwise if projection workspace already exists
      } else {
        # load existing R workspace
        load(workspace_name_out0)
        # sign-posting indicating models are loaded from past runs
        cat('\n', sp_nm, 'projection of individual models loaded from past run.')
      }
      
      # run BIOMOD2 fixes if TRUE in source script
      if (apply_biomod2_fixes){
        # load baseline projections manually from directory
        myBiomodProjection<-LoadProjectionManually(myBiomodProj_baseline)        
      } else {
        # use previously created baseline projections
        myBiomodProjection<-myBiomodProj_baseline
      }
      
      # sign-posting to run forecasting of ensemble models
      cat('\n run ensemble forecasting...')
      
      # set location to save R workspace data and ensemble projections
      workspace_name_out1 = paste0(sp_is, "_FB_all_model_proj_EF", proj_nm, ".RData") 
      # check if ensemble projections workspace already exists or not
      if (file.exists(workspace_name_out1) == F | overwrite == T){ #allow to overwrite
        
        
        ################################
        ##### ensemble_forecasting #####
        ################################
        # run ensemble projections for species
        myBiomodEF <- BIOMOD_EnsembleForecasting(
          projection.output = myBiomodProjection,  #BIOMOD.projection.out from projections
          total.consensus = TRUE,  #mean of all combined model projections
          EM.output = myBiomodEM,  #BIOMOD.EnsembleModeling.out from ensemble modeling
          binary.meth = eval_stats,  #evaluation method statistics 
          keep.in.memory = memory)  #if output should be saved to hard disk or not
        
        # sign-posting of completed ensemble projections
        cat('\n', sp_nm, 'ensemble projection done.') 
        
        # save ensemble projections R workspace environement 
        save("myBiomodProjection", "myBiomodEF", file = workspace_name_out1)
        
        # remove any temporary files
        removeTmpFiles(h=1)
        
        # sign-posting of projection done for species per island
        cat('\n' ,sp_nm, "_", spIsland,'done.')
      } else {
        # sign-posting if projection for species is already done
        cat('\n', sp_nm, "_", spIsland, 'previously done.')
      }   
    }  # -=-=- end projections per individual islands -=-=-#
    
    #####################################################
    ##### MERGE RASTERS PER ISLAND INTO ARCHIPELAGO #####
    #####################################################
    # if projections by island have been previously done
    if (proj_by_island){
      # sign-posting of merging individual island forecasts 
      cat('\n','merging island by island forecasts')

      # combine species data found on more than one island 
      if (length(spIslandNames) > 1) { 
        # sign-posting of merging single island projections 
        cat('\n',"merging individual island projections for ", sp_nm)
        
        # set directory to save combined output rasters
        combinedDir = paste0(working_dir, sp_dir, "/proj_", proj_nm) 
        # create directory if does not already exist
        dir.create(combinedDir, showWarnings = FALSE) 
        # set first island directory
        island1Dir = paste0(working_dir, sp_dir,"/proj_", proj_nm,"_", spIslandNames[1],"/") 
        # list all files within single island directory
        fileList<-list.files(island1Dir, pattern = "*.grd") 
        # loop through each file in island directory
        for (file in fileList) {
          # create raster stack from .gri files in directory list
          rasterStack<-stack(paste0(island1Dir, file)) 
          # assign names to combined raster stack
          rasterStackNames = names(rasterStack) 
          # for each additional island 
          for (islandNum in 2:length(spIslandNames)) { 
            # sign-posting listing combining of island rasters per species
            cat('\n doing', file, spIslandNames[islandNum])
            # name new island being processed
            newIslandName = spIslandNames[islandNum] 
            # create temporary island directory 
            tempIslandDir = str_replace_all(island1Dir, spIslandNames[1], newIslandName) 
            # list files in temporary island directory
            tempFile = str_replace_all(file, spIslandNames[1], newIslandName) 
            # create raster stack from .gri files in temporary directory list
            tempStack<-stack(paste0(tempIslandDir, tempFile)) 
            # merge each island with first  raster stack created above for species 
            rasterStack = stack(merge(rasterStack, tempStack)) 
          }
          # assign names to merged raster stack
          names(rasterStack)<-rasterStackNames 
          # create file name for merged raster stack
          combinedFileName = str_replace_all(file, paste0("_", spIslandNames[1]), "") 
          # set file location for mergerd raster stack
          combinedFileLoc = paste0(combinedDir, "/", combinedFileName) 
          # save merged raster stack file to location
          writeRaster(rasterStack, combinedFileLoc, overwrite = TRUE) 
        }
      }
    }  # -=-=- end merge raster per island into one -=-=- #
    
    # otherwise combined raster was created above
  } else {
    # sign-posting that raster was previously calculated
    cat('\n',sp_nm,'previously calculated.')
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
##### END MODEL PROJECTIONS #####
#################################