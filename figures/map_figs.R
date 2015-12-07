# map figures for species distribution results 
# based on ensemble models in 'Y:/' directory

library(raster)

rootDir<-"Y:/PICCC_analysis/FB_analysis/model_results/biomod2/"



### FIGURE 1 ###
# archipelago-wide present and future distribtuions of Kauai endemic species
isDir<-paste0(rootDir, "FB_ALL_100runs/")
setwd(isDir)

# subset of Hawaii (11) and Maui (8) native forest bird species 
sp_nm = c('Akekee', 'Akikiki')

past<-"/proj_baseline"
pred<-"/proj_future"

# main hawaiian islands outline shapefile
mask_layer = shapefile("Main_Hawaiian_Islands_simple3.shp")

par(mar = c(2.5, 3, 2, 1), mfrow=c(2, 1))

# AKEKEE
plot(raster((paste0(sp_nm[1], past, past, "_", sp_nm[1], "_ensemble_ROCbin.grd"))),
     col = c("gray",  heat.colors(1)))
par(new = T)
plot(raster((paste0(sp_nm[1], pred, pred, "_", sp_nm[1], "_ensemble_ROCbin.grd"))), 
     col = c(NA, "darkgreen"))
title("(a) Akekee Species Distribution", line = -1, cex.main = 1)
# legend 
legend("bottomleft", legend = c("present", "projected"), col = c("red", "darkgreen"), pch = 15)

# main plot title
title("Archipelago-Wide Distribution of Akekee and Akikiki", line = 0.5, cex.main = 1.5)

# AKIKIKI
plot(raster((paste0(sp_nm[2], past, past, "_", sp_nm[2], "_ensemble_ROCbin.grd"))),
     col = c("gray",  heat.colors(1)))
par(new = T)
plot(raster((paste0(sp_nm[2], pred, pred, "_", sp_nm[2], "_ensemble_ROCbin.grd"))), 
     col = c(NA, "darkgreen"))
title("(b) Akikiki Species Distribution", line = -1, cex.main = 1)
# legend 
legend("bottomleft", legend = c("present", "projected"), col = c("red", "darkgreen"), pch = 15)


### FIGURE 2 ###
# future projections for Akekee and Akikiki on Maui and Hawaii Island
maDir<-paste0(rootDir, "FB_MA_projection_100_runs/output_rasters/main/")  # Maui
biDir<-paste0(rootDir, "FB_BI_projection_100_runs/output_rasters/main/")  # Big Island
setwd(maDir)

plot(raster((paste0(sp_nm[1], "_response_zones_ROC_wmean_future.tif"))), 
     col = c(terrain.colors(1), "lightgreen"))

