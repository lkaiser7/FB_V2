
if (toCompareWithCurrentClimate==1){
  clim_surface_to_use=clim_data_2000 
  proj_nm0='baseline'}
if (toCompareWithCurrentClimate==2){
  clim_surface_to_use=clim_data_2000wettest
  proj_nm0='baseline_wettest'}
if (toCompareWithCurrentClimate==3){
  clim_surface_to_use=clim_data_2000driest
  proj_nm0='baseline_driest'}
if (toCompareWithCurrentClimate==4){
  clim_surface_to_use=clim_data_2100 
  proj_nm0='future'}
if (toCompareWithCurrentClimate==5){
  clim_surface_to_use=clim_data_2100wettest
  proj_nm0='future_wettest'}
if (toCompareWithCurrentClimate==6){
  clim_surface_to_use=clim_data_2100driest
  proj_nm0='future_driest'}
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

setwd(working_dir)
library(biomod2)
library(stringr)

#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
#memory.limit
(size=240000)
#sp_nm=spp_nm[1]
spp_info=read.csv(paste(csv_dir,'FB_spp_data_HI.csv', sep = ""))

dir.create("output_rasters/",showWarnings=F)
dir.create("output_rasters/analog_clim/",showWarnings=F)


####create analog climate map
sp_nm=spp_nm[1]
n_abs_removed=c()
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,' analog climate mapping...')
  raster_name=paste0("output_rasters/analog_clim/", sp_nm, "_analog_clim_", comp_projects[2])
  tif_name=paste(raster_name,".tif", sep="")
  jpeg_name=paste(raster_name,".jpg", sep="")
  if (file.exists(tif_name)==F | overwrite==1){ #check to see if the analysis for this species was already done    
    all_clim_data=c("clim_data_2000", "clim_surface_to_use")
    clim_stacks=c("biovars2000", "biovars2100")
    
    sp_index=which(spp_info[,"Species"]==sp_nm)
    raster_res= spp_info[sp_index,"rasterdir"]
    clim_data_dir=fitting_clim_data_dir 
    jnk0=length(env_var_files)
    crop_raster=raster(paste(crop_raster_dir,raster_res,".grd",sep=""))
    for (dfd in 1:length(all_clim_data)){
      clim_data=all_clim_data[dfd]
      clim_data_dir=get(clim_data)
      predictors_temp = raster( paste(clim_data_dir, env_var_files[1], sep=""))
      predictors_temp=crop(predictors_temp,  crop_raster)
      for (jj in 2:jnk0){
        temp=raster(paste(clim_data_dir, env_var_files[jj], sep=""))
        temp=crop(temp,  crop_raster)
        predictors_temp = addLayer(predictors_temp, temp)
      }
      names(predictors_temp)<- var_name
      assign(clim_stacks[dfd], predictors_temp)
    }
    
    jnk=raster(paste(clim_data_dir, "bio_12.tif", sep=""))
    jnk=crop(jnk,  crop_raster)
    Response_var  =  reclassify(jnk,  c(0.1,Inf,1))
    
    pred <- sre(Response_var, biovars2000, biovars2100, Quant = 0.0001)
    analog_climates2100=subset(pred, 1)
    rm("crop_raster" ,"temp")         
    
    jpeg(jpeg_name,
         width = 10, height = 10, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(analog_climates2100)
    dev.off()

    writeRaster(analog_climates2100, tif_name, format="GTiff", overwrite=TRUE)
    
    cat('\n',sp_nm,'analog climate raster created...')
  }
}



### extra coding for map figure output
## open state coastlines for mapping
#mapDir<-"C:/Users/lkaiser/Dropbox/FB_SDM/necessary_run_data/coastlines"
## projection coordinates for shapefiles 
#coords<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
## 100 m equivalency factor
#equiv_100m = 0.0009430131  
## open shpapefile for island coastlines 
#state<-readOGR(mapDir, "state_coast", p4s = coords)
## create tiff file of analog climate map output
#tiff(file = "output_rasters/analog_clim/Kauai_climate_envelope.tif")
#plot(pred, legend = FALSE, 
     #main = "Climate Envelope of Kaua`i", 
     #xlab = "Longitude", ylab = "Latitude")
#plot(state, add = T)
#legend("bottomleft", legend = "Current Analog Climates",
       #pch = 15, col = "darkgreen")
#dev.off()
### end extra coding