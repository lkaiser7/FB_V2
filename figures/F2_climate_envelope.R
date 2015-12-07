# FB Manuscript V2 - Figure 2
# archipelago-wide analog climate envelop mapping
# baseline (current) Kauai climate mapped on other islands
# from script 4 from ensembles

# load necessary packages
library("raster")
library("rgdal")
library("biomod2")
library("stringr")

# set species names to be analyzed
sp_nm<-c("Akekee", "Akikiki")

# climate mapping root directory
root_Dir<-"Y:/PICCC_analysis/FB_analysis/model_results/biomod2/"
# set working directory
setwd(root_Dir)
# set output folder location
out_Dir<-paste0(root_Dir, "clim_map/output_rasters/analog_clim/")

# raster directory
crop_raster_dir<-paste0(root_Dir, "clim_map/map_crop/")
# raster grid for Kauai extent
kauai_ras = raster(paste0(crop_raster_dir, "Kauai Endemics.grd"))
# raster grid for archipelago extent
arch_ras = raster(paste0(crop_raster_dir, "Apapane_Iiwi.grd"))

# list environmental variable files
env_var_files<-c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif")

# set fitting climate data at 125 m
clim_data_dir<-paste0(root_Dir, "necessary_run_data/clim_data/all_baseline/125m/")
# stack fitting climate variables at 125 m
fit_stack<-stack(paste0(clim_data_dir, env_var_files))
# crop fit variables to Kauai
fit_crop<-crop(fit_stack, kauai_ras)

# set baseline climate data at 500 m 
clim_data<-paste0(root_Dir, "necessary_run_data/clim_data/all_baseline/500m/")
# stack 4 bioclimatic factors for the archipelago
predictors<-stack(paste0(clim_data, env_var_files))

# create raster mask from bio12 variable
jnk_ras<-raster(paste0(clim_data_dir, "bio12.tif"))
# crop mask to Kauai extent
jnk_ras<-crop(jnk_ras,  kauai_ras)
# reclassify values of raster mask
resp_var<-reclassify(jnk_ras, c(0.1, Inf, 1))

# ?sre
# create surface range envelope using raster maks, fitted Kauai climate variables, and all predictors
sre_pred<-sre(resp_var, fit_crop, predictors, Quant = 0.0001)
# plot climate envelope output
plot(sre_pred)

# open state coastlines for mapping (*do not add / to end of map path!)
mapDir<-paste0(root_Dir, "necessary_run_data/Coastlines")
# projection coordinates for shapefiles 
coords<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
# 100 m equivalency factor
equiv_100m = 0.0009430131  
# open shpapefile for island coastlines 
state<-readOGR(mapDir, "state_coast", p4s = coords)

# plot and save jpeg image
jpeg(file = paste0(out_Dir, "Kauai_CE.jpeg"), res = 300, units = "in",
     width = 10, height = 10, pointsize = 12, quality = 90)
plot(sre_pred, legend = FALSE, main = "Climate Envelope of Kaua`i", 
     cex.main = 2.5, xlab = "Longitude", ylab = "Latitude", cex.lab = 1.5)
plot(state, add = TRUE, border = "slategray")
legend("bottomleft", legend = "Current Analog Climates", 
       pch = 15, col = "darkgreen")
dev.off()

# plot and save tiff file
tiff(file = paste0(out_Dir, "Kauai_CE.tif"), res = 300, units = "in",
     width = 10, height = 10, pointsize = 12, compression = "lzw")
plot(sre_pred, legend = FALSE, main = "Climate Envelope of Kaua`i", 
     cex.main = 2.5, xlab = "Longitude", ylab = "Latitude", cex.lab = 1.5)
plot(state, add = TRUE, border = "slategray")
legend("bottomleft", legend = "Current Analog Climates",
       pch = 15, col = "darkgreen", cex = 1.5)
dev.off()
