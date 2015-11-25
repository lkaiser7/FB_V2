### variable importance graph ###
### use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# load necessary packages
library(ggplot2)
library(reshape2)
library(plyr)

# open variable importance csv file
varImp0 = read.csv('all_VariImp.csv')
# store column names
colnames = names(varImp0)
# create empty list of models
model_cols = list()

# loop through all models
for (model in models_to_run){
  # store each model as a temporary object
  jnk = grep(paste0(model, "+"), colnames, perl = TRUE, value = FALSE)  
  # add each model to empty list created
  model_cols[[length(model_cols)+1]]=jnk
}

# select first species
sp_nm = spp_nm[1]

# create directory in output rasters if it does not already exist
dir.create("output_rasters/varImp/", showWarnings = FALSE)

# loop through all species
for (sp_nm in spp_nm){
  # loop through all models
  for (model in models_to_run){
    # select which models to analyze from runs
    model_n = which(models_to_run == model)
    # select variable importance per model
    varImp = varImp0[varImp0[,1] == sp_nm,model_cols[[model_n]]]
    # transform the data in a data frame
    varImp = as.data.frame(t(varImp))
    # rename columns by models
    names(varImp) = varImp0[varImp0[,1] == sp_nm, "rownames.Spp_VariImp."]
    # reshape and stack the data
    long = melt(varImp)
    # rename columns
    names(long) = c("Predictor", "Value")
    
    jpeg_name=paste0("output_rasters/varImp/", sp_nm, "_",
                     model, "_variable_importance_box_plot.jpg")
    a = qplot(Predictor, Value, data = long, geom = c("boxplot", "jitter"), 
            fill = Predictor, main = "", xlab = "", ylab = "Variable importance" )
    ggsave(filename = jpeg_name, plot = a)
    
    cat('done with', sp_nm, model, 'variable importance box plot \n')
  }
}

# load and format data in loop as done previously done above for violin plots
sp_nm=spp_nm[1] 
dir.create("output_rasters/varImp/", showWarnings=F)
for (sp_nm in spp_nm){
  for (model in models_to_run){
    model_n=which(models_to_run==model)
    varImp=varImp0[varImp0[,1]==sp_nm,model_cols[[model_n]]]
    varImp=as.data.frame(t(varImp))
    names(varImp)=varImp0[varImp0[,1]==sp_nm,"rownames.Spp_VariImp."]
    long=melt(varImp)
    names(long)=c("Predictor", "Value")
    
    jpeg_name=paste0("output_rasters/varImp/", sp_nm, "_",
                     model, "_variable_importance_violin_plot.jpg")
    a = qplot(Predictor, Value, data = long, geom = c("violin"),
            fill = Predictor, main = "", xlab = "", ylab = "Variable importance" )
    a = a + geom_violin(scale = "width")
    a = a + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    ggsave(filename = jpeg_name, plot = a)  
    
    cat("done with ", sp_nm, " ",model, " variable importance box plot", '\n')
  }
}

####################################
##### END MODEL FITTING GRAPHS #####
####################################