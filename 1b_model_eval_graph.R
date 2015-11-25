### evaluation statistics graphs ###
### to use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# set working directory 
setwd(working_dir)
# create directory in output rasters if it does not already exist
dir.create("output_rasters/evalMat/", showWarnings = FALSE)

# load necessary packages
library(ggplot2)
library(reshape2)
library(plyr)

# select first evaluation statistic
eval_stat = eval_stats[1]
# select first model
model = models_to_run[1]

# loop through all evaluation statistics
for (eval_stat in eval_stats){
  # open stored evaluation statistics csv file
  evalMat0 = read.csv(paste0('all_eval_mat_',eval_stat,'.csv'))
  # loop through all models
  for (model in models_to_run){
    # select evaluation statistics per model
    evalMat = evalMat0[evalMat0[,2]==model,3:ncol(evalMat0)]
    # transform the data in a data frame
    evalMat = as.data.frame(t(evalMat))
    # rename columns by models
    names(evalMat) = evalMat0[evalMat0[,2]==model,1]
    # reshape and stack the data
    long = melt(evalMat)
    # rename columns
    names(long) = c("Species", "Value")
    # select data per species
    long=long[long[,1] %in% spp_nm,]
    # reorder data frame
    long=long[order(long[,1]),]
    # sort data frame per species as factor levels
    long$Species<-factor(long$Species, levels = sort(levels(long$Species)))
    
    # store output jpeg file
    jpeg_name = paste0("output_rasters/evalMat/",eval_stat, 
                     "_",model, "_variable_importance_box_plot.jpg")
    # store basic qplot boxplot
    a = qplot(Species, Value, data = long, geom = c("boxplot"), fill = Species, main = "",
            xlab = "", ylab = paste(eval_stat, "model evaluation")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    # remove legend
    a = a + theme(legend.position = "none")
    # save jpeg file
    ggsave(filename = jpeg_name, plot = a)   
    
    # sign-posting of clompleted box plot
    cat('done with', eval_stat, model, 'evaluation statistic box plot \n')
  }
}

# load and format data in loop as done previously done above for violin plots
for (eval_stat in eval_stats){
  evalMat0 = read.csv(paste0('all_eval_mat_', eval_stat, '.csv'))
  for (model in models_to_run){
    evalMat = evalMat0[evalMat0[,2] == model,3:ncol(evalMat0)]
    evalMat = as.data.frame(t(evalMat))
    names(evalMat) = evalMat0[evalMat0[,2] == model,1]
    lon = melt(evalMat)
    names(long) = c("Species", "Value")
    long = long[long[,1] %in% spp_nm,]
    long = long[order(long[,1]),]
    long$Species<-factor(long$Species, levels = sort(levels(long$Species)))
    
    jpeg_name = paste0("output_rasters/evalMat/", eval_stat, 
                     "_", model, "_variable_importance_violin_plot.jpg")
    a = qplot(Species, Value, data=long, geom=c("violin"), fill = Species, main = "",
            xlab="", ylab=paste(eval_stat, "model evaluation"))
    #+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    a = a + geom_violin(scale = "width")
    a = a + guides(fill = guide_legend(keywidth = 1, keyheight = 1.2))
    a = a + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    a = a + coord_flip()
    a = a + guides(fill = guide_legend(reverse=TRUE), guide = guide_legend(title = NULL))
    a
    ggsave(filename = jpeg_name, plot = a)    
    
    cat('done with', eval_stat, model, 'evaluation statistic box plot \n')
  }
}

####################################
##### END MODEL FITTING GRAPHS #####
####################################