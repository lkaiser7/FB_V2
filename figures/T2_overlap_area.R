# FB Manuscript V2 - Table 2
# overlapping habitat area of endemic Kauai Forest Birds
# calculating overlap area and percent with native birds
### CHANGE LINES 25-31 FOR DIFFERENT ISLANDS & SCENARIOS


# load necessary packages
library("raster")

# set root and output directories
rootDir<-"Y:/PICCC_analysis/FB_analysis/model_results/biomod2/"
outDir<-paste0(rootDir, "TabFig_Outputs/")
setwd(rootDir)

# set island and archipelago run paths
maDir<-paste0(rootDir, "FB_MA_projection_100_runs/")
biDir<-paste0(rootDir, "FB_BI_projection_100_runs/")

# set island species (ma = 8, bi = 11)
ma_sp = c('Akekee', 'Akikiki', 'Akohekohe', 'Apapane', 'Hawaii_Amakihi',  
          'Iiwi', 'Maui_Alauahio', 'Maui_Parrotbill')
bi_sp = c('Akekee', 'Akikiki', 'Akiapolauu', 'Apapane', 'Hawaii_Akepa', 'Hawaii_Amakihi',
          'Hawaii_Creeper', 'Hawaii_Elepaio', 'Iiwi', 'Omao', 'Palila')

### CHOOSE: island run [Maui or Hawaii]
isDir<-maDir
is_sp<-ma_sp
is_run<-"ma"

### CHOOSE: baseline or future run
bof<-"baseline"

# set raster location 
ras<-paste0(isDir, "output_rasters/")

# set vegetation biome limits for species ranges
vegDir<-"Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/current_veg_mask/"

# create directory to save rasters
dir.create(paste0(outDir, is_run, "_", bof, "_cropped_binary_rasters"))
rasDir<-paste0(outDir, is_run, "_", bof, "_cropped_binary_rasters/")

# loop through species and create binary rasters
for (i in 1:length(is_sp)) {  # set i = 1 for debugging single run
  # open suitability tiff file
  sp_rast<-raster(paste0(ras, is_sp[i],"_clipped_suitability_", bof, "_ROC_wmean.tif"))

  # open vegetation overlay mask archipelago-wide for species 
  veg_mask<-raster(paste0(vegDir, is_sp[i],"_current_veg_mask.tif"))
  
  # calculate ratio to crop raster
  veg_ratio = round(res(sp_rast)[1]/res(veg_mask)[1])
  if (veg_ratio > 1) {
    veg_mask<-aggregate(veg_mask, fact = veg_ratio, fun = max)
  }
  
  # crop vegetation mask to raster extent for island run
  veg_mask<-crop(veg_mask, sp_rast)
  # resample values between rasters using nearest neighbor method
  veg_samp = resample(veg_mask, sp_rast,  method="ngb")
  # crop to only include avaiable vegetation 
  resp_sp = sp_rast*veg_samp
  
  # convert all species range to binary values 
  resp_sp[resp_sp[] == 0]<-NA
  resp_sp[resp_sp[] > 0]<-1
  
  # save binary raster file
  writeRaster(resp_sp, overwrite = T, filename = paste0(rasDir, is_sp[i], "_raster_area.grd"))
}

# select rasters (.grd) files
f_names<-list.files(rasDir)
grd_files<-grep(".grd", f_names, value = TRUE)

# isolate Kauai endemic rasters
ake<-raster(paste0(rasDir,grep("Akekee", grd_files, value = TRUE)))
aki<-raster(paste0(rasDir,grep("Akikiki", grd_files, value = TRUE)))

# loop through rasters(.grd files) and calculate habitat overlap area
for (j in 1:length(grd_files)) {  # set j = 1 for debugging single run
  # open first raster file
  sp_comp<-raster(paste0(rasDir, grd_files[j]))
  
  # sum species raster with Kauai endemics to determine overlap area
  ake_sum<-sum(ake, sp_comp)
  aki_sum<-sum(aki, sp_comp)
  
  # plot and save overlap area for visual confirmation
  jpeg(file = paste0(outDir, "overlap_maps/", is_run, "_", bof, "_", is_sp[j], "_overlap.jpg"))
  par(mfrow=c(2,1))
  # Akekee
  plot(ake, col = "yellow", main = paste(is_sp[j], "vs. Akekee", sep = " "))
  plot(sp_comp, add = T, col = "blue")
  # overlap
  plot(ake_sum, add = T, col = "green")
  legend("bottomleft", legend = c("Akekee", is_sp[j], "Overlap"),
         col = c("yellow", "blue", "green"), pch = 15)
  # Akikiki
  plot(ake, col = "yellow", main = paste(is_sp[j], "vs. Akikiki", sep = " "))
  plot(sp_comp, add = T, col = "blue")
  # overlap
  plot(aki_sum, add = T, col = "green")
  legend("bottomleft", legend = c("Akikiki", is_sp[j], "Overlap"),
         col = c("yellow", "blue", "green"), pch = 15)
  dev.off()
  
  # create object with frequency/4 to calculate overlap area (sp_area, ake_ol, aki_ol)
  sp_ol<-data.frame(grd_files[j], 
                    freq(sp_comp, useNA = 'no')[2]/4, 
                    freq(ake_sum, useNA = 'no')[2]/4, freq(aki_sum, useNA = 'no')[2]/4)
  
  # calculate percent overlap of habitat area
  sp_ol[1 , 5]<-cbind(sp_ol[3]/sp_ol[2]*100)
  sp_ol[1 , 6]<-sp_ol[4]/sp_ol[2]*100
  names(sp_ol)<-c("sp_file", "sp_area", "ake_ol", "aki_ol", "ake_pct", "aki_pct")
  
  # build output table per island and scenario
  if (j == 1) {
    ol_dat<-sp_ol
    names(ol_dat)<-c("sp_file", "sp_area", "ake_ol", "aki_ol", "ake_pct", "aki_pct")
  } else {
    if (j > 1) {
      ol_dat<-rbind(ol_dat, sp_ol)
    }
  }  
}
# save final table with all species per scenario
write.csv(ol_dat, file = paste0(outDir, is_run, "_", bof, "_habitat_overlap_area.csv"))


### END SCRIPT: REPEAT LOOP 4x FOR OTHER SCENARIOS (baseline or future) AND ISLANDS (maui or big island)