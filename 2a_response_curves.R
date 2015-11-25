### ensemble response curves ###
### after 2_FB_SDM_ensembles ### 

##################
##### SET UP #####
##################

# load necessary packages
library(biomod2)
library(stringr)

# select first models run by last name
last_model = models_to_run[length(models_to_run)]

# create directories for outputs
dir.create("output_rasters/", showWarnings = FALSE)
dir.create("output_rasters/response_curves/", showWarnings = FALSE)
dir.create("output_rasters/response_curves/combo/", showWarnings = FALSE)

# loop through all species
for (sp_nm in spp_nm){
  # convert species name to character
  sp_nm = as.character(sp_nm)
  # set initial species name
  sp_nm0 = sp_nm
  
  # sign-posting of ongoing modeling
  cat('\n',sp_nm,'modeling...')
  
  # set name of file to save workspace output data after model run
  workspace_name = paste0(sp_nm,"_FB_modelfitting.RData") 
  # check to see if file already exists
  if (file.exists(workspace_name)){
    # if so, load workspace
    load(workspace_name)    
  }else{
    # otherwise set name of file to load workspace data from model run    
    workspace_name = paste0(sp_nm0,"_FB_run.RData")
    # load workspace from previous ensemble modeling
    load(workspace_name)    
  }

  # replace any species names with "_" to "."
  sp_nm = str_replace_all(sp_nm,"_", ".") 
  
  # loop through all models
  model = models_to_run[1]
  # for each model
  for (model in models_to_run){
    # set temporary name for jpeg image file
    temp_jpeg_name = paste0('output_rasters/response_curves/combo/', sp_nm, "_", 
                            "response_curves", "_", model, "_all_vars.jpg")
    # check to see if jpeg file already exists
    if (file.exists(temp_jpeg_name) == FALSE | overwrite == TRUE){
      # if so, load BIOMOD2 model output object per model
      loaded_models<-BIOMOD_LoadModels(myBiomodModelOut, models = model)
      # remove last column
      loaded_models = loaded_models[1:length(loaded_models) -1]
      
      # plot 2D response plots
      myRespPlot2D<-response.plot2(models = loaded_models,
                                   Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                   show.variables= get_formal_data(myBiomodModelOut,
                                                                   'expl.var.names'),
                                   do.bivariate = FALSE,
                                   fixed.var.metric = 'mean',
                                   save.file="no", 
                                   name="response_curve", 
                                   ImageSize=480, 
                                   plot=FALSE)
      # NOTE: "myRespPlot2D" is now a list object rather than an array ###
      
      # set to 1 for debugging
      bioclim_cnt = 1 
      # loop through each bioclimatic variable
      for (bioclim_cnt in 1:length(myRespPlot2D)) {
        # set y limits
        ymax_lim = 0
        ymin_lim = 10000
        # set y limits per bioclimatic variable
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
          ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
        }
        # set x limits per bioclimatic variable
        xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
        xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
        
        # rename variable names for plot
        var_name = names(myRespPlot2D)[bioclim_cnt]
        
        # save name of jpeg file for response curves
        jpeg_name = paste0('output_rasters/response_curves/', sp_nm, "_", 
                           "response_curve", "_", model, "_", var_name, ".jpg")
        # create blank jpeg file in directory 
        jpeg(jpeg_name, width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        
        # take mean of for first bioclimatic variable
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        # repeat for all other bioclimatic variables
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          # take mean of each additional bioclimatic variable
          pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
          
          # for second bioclimatic variable
          if (rep == 2){
            # plot bioclimatic variable response curve
            plot(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), ylab = "Response",
                 ylim = c(ymin_lim, ymax_lim), xlab = var_name, col = "grey")
          } else {
            # otherwise add additional variable lines to response curve
            lines(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), 
                  ylim = c(ymin_lim, ymax_lim), col = "grey")
          }
        }
        
        # average response curve
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        # average all other variable response curves
        pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
        # add average response curve line to plot
        lines(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), 
              ylim = c(ymin_lim, ymax_lim), lwd = 3)    
        # save jpeg file
        dev.off()
      }
      
      ########################
      ##### combo figure #####
      ########################
      # use same formatted setting above to combine response variable curve plots
      jpeg_name_combined = paste0('output_rasters/response_curves/combo/', sp_nm, "_", 
                                  "response_curves","_", model, "_all_vars.jpg")
      jpeg(jpeg_name_combined, width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300) 
      ymax_lim = 0
      ymin_lim = 10000      
      for (bioclim_cnt in 1:length(myRespPlot2D)) {=
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
          ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
        }
      }
      
      # set graphical parameters fro combined plot
      par(mfrow = c(2,2), oma = c(0,0,3,0))
      
      bioclim_cnt = 1
      for (bioclim_cnt in 1:length(myRespPlot2D)) {
        xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
        xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
        var_name = names(myRespPlot2D)[bioclim_cnt]
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
          pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
          if (rep == 2){
            plot(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), ylab = "Response", 
                 ylim = c(ymin_lim, ymax_lim), xlab = var_name, col = "grey")
          } else {
            lines(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), 
                  ylim = c(ymin_lim, ymax_lim), col = "grey")
          }
        }
        
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
        lines(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), 
              ylim = c(ymin_lim, ymax_lim), lwd = 3)          
      } 
      
      # create temporary title for model response curves per species 
      temp_title = paste0(model, " modeled response for ", sp_nm)
      # add temporary text title to plot
      mtext(temp_title, adj = 0.5, side = 3, outer = TRUE) 
      dev.off()
      
      # sign-posting of completed response curve plots
      cat('\n response curves for ',sp_nm, model, 'finished.')  
    }else{
      # sign-posting if response curve plots have already been done
      cat('\n response curves for ',sp_nm, model, 'already done.')  
    }      
  }
}

####################################
##### END RESPONSE CURVE PLOTS #####
####################################