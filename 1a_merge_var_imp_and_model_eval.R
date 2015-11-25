### tables from model fitting ###
### use after 1_model_fitting ### 

##################
##### SET UP #####
##################

# load necessary packages
library(biomod2)
library(reshape2)
library(plyr)

# set sp_nm to the first value
sp_nm = spp_nm[1]
# set i counter to 1
i = 1

# loop through all species to get variable importance and model evalutation data
for (sp_nm in spp_nm){
  # set file name of saved workspace data from model fitting
  workspace_name = paste(sp_nm,"_FB_modelfitting.RData", sep = "") 
  # load model fitting workspace data
  load(workspace_name)
  # get model evaluations
  myBiomodModelEval<-getModelsEvaluations(myBiomodModelOut)  
  # return dimension names of model evaluations 
  dimnames(myBiomodModelEval)
  
  # loop through all evaluation statistics 
  for (eval_stat in eval_stats){
    # select evaluation statistic values
    myBiomodModelEval[eval_stat,"Testing.data",,,]
    # store statistic values as a data frame
    Spp_eval<-data.frame(myBiomodModelEval[eval_stat,"Testing.data",,,])
    # combine with column and row names for metrics and model runs
    Spp_eval = cbind(matrix(sp_nm,dim(Spp_eval)[1],1),rownames(Spp_eval),Spp_eval)
    
    # for the first species in the loop
    if (i == 1){
      # create a new data frame to store all values 
      assign(paste0("all_eval_mat_", eval_stat), Spp_eval)
      # otherwise for all other species in the loop 
    }else{
      # create a temporary object to store evaluation statistics
      jnk = rbind(get(paste0("all_eval_mat_", eval_stat)), Spp_eval)
      # add additional statistic values to created data frame
      assign(paste0("all_eval_mat_", eval_stat), jnk)      
    }
  }
  
  # select variable importance values
  getModelsVarImport(myBiomodModelOut)
  # store variable importance values as a data frame
  Spp_VariImp<-data.frame(getModelsVarImport(myBiomodModelOut))
  # combine column and row names for variables and model runs
  Spp_VariImp=cbind(matrix(sp_nm,dim(Spp_VariImp)[1],1),rownames(Spp_VariImp),Spp_VariImp)
  
  # for the first species in the loop
  if (i == 1){
    # create a new data frame to store all values
    all_var_imp_mat = Spp_VariImp
    # otherwise for all other species in the loop
  }else{
    # add additional variable importance values to created data frame
    all_var_imp_mat = rbind(all_var_imp_mat, Spp_VariImp)
  }
  
  # add 1 to the i counter for next loop
  i = i + 1
}

# store file name for all variable importance values
FileName<-paste("all_VariImp.csv")
# write variable importance file
write.table(all_var_imp_mat, file = FileName, sep=",", row.names = FALSE)

# take the mean of all the importance values per variable
all_var_imp_mean = all_var_imp_mat[,1:2]
# combine mean values
all_var_imp_mean = cbind(all_var_imp_mean, meanVarImp = 
                         apply(all_var_imp_mat[,3:dim(all_var_imp_mat)[2]],1,mean, na.rm=T))
# rename columns 
names(all_var_imp_mean)=c("species","var","meanVarImp")
# remove row names
row.names(all_var_imp_mean)<-NULL 
# reshape and aggregate the data per variable
all_var_imp_mean = dcast(all_var_imp_mean, species ~ var, value.var = "meanVarImp")

# store file name for mean variable importance values
FileName<-paste("all_VariImp_mean.csv")
# write variable importance means file
write.table(all_var_imp_mean, file = FileName, sep=",", row.names = FALSE)

# loop through all evaluation metrics
for (eval_stat in eval_stats){
  # store file name for all evaluation statistcs
  FileName<-paste0("all_eval_mat_",eval_stat,".csv")
  # write evaluation statistics file
  write.table(get(paste0("all_eval_mat_",eval_stat)), 
              file = FileName, sep=",", row.names = FALSE)

  # retrun evaluation statistic object
  tmp_eval_map = get(paste0("all_eval_mat_",eval_stat))
  # rename columns
  names(tmp_eval_map)[c(1:2)]=c("species","model")
  # reformat the statistics
  tmp_eval_map<-reshape(tmp_eval_map, timevar=c("model"), idvar=c("species"), dir="wide")
  # select the first column
  tmp_eval_map2 = tmp_eval_map[,1]
  # create data frame with mean values
  tmp_eval_map2 = data.frame(species = tmp_eval_map2, meanEval = 
                             apply(tmp_eval_map[,2:dim(tmp_eval_map)[2]],1,mean, na.rm=T))
  # remove row names
  row.names(tmp_eval_map2)<-NULL 
  # save temporary object with column names
  # jnk = names(tmp_eval_map)[c(3:dim(tmp_eval_map)[2])]
  # reshape and aggregate data
  # all_var_imp_mean = dcast(all_var_imp_mean, species ~ model, value.var=jnk, na.rm=T)
  # store file name for mean evaluation statistics
  FileName<-paste0("all_eval_mean_mat_",eval_stat,".csv")
  # write mean evaluation statistics file
  write.table(tmp_eval_map2, file = FileName, sep=",", row.names = FALSE)
}

####################################
##### END MODEL FITTING TABLES #####
####################################