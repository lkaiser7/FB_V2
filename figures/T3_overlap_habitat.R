# FB Manuscript V2 - Table 3
# overlapping habitat area for endemic Kauai Forest Birds
# calculating km^2 of overlap area with native birds
### CHANGE LINES 24 - 27 FOR DIFFERENT ISLANDS

# load necessary packages
library("raster")

# set root and output directories
rootDir<-"Y:/PICCC_analysis/FB_analysis/model_results/biomod2/"
outDir<-paste0(rootDir, "TabFig_Outputs/")
setwd(rootDir)

# set island and archipelago run paths
maDir<-paste0(rootDir, "FB_MA_projection_100_runs/")
biDir<-paste0(rootDir, "FB_BI_projection_100_runs/")

# set island species 
ma_sp = c('Akohekohe', 'Apapane', 'Hawaii_Amakihi',  
          'Iiwi', 'Maui_Alauahio', 'Maui_Parrotbill')
bi_sp = c('Akiapolauu', 'Apapane', 'Hawaii_Akepa', 'Hawaii_Amakihi',
          'Hawaii_Creeper', 'Hawaii_Elepaio', 'Iiwi', 'Omao', 'Palila')

### CHOOSE: island run [Maui or Hawaii]
isDir<-maDir
is_sp<-ma_sp
is_run<-"ma"

# set raster location 
ras<-paste0(isDir, "output_rasters/")

# set species names (Kauai endemics)
sp_nm<-c("Akekee", "Akikiki")

# set scenario names (current or projected)
bof<-c("baseline", "future")

# set vegetation biome limits for species ranges
vegDir<-"Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/current_veg_mask/"

# open raster and crop for Kauai endemics and save
for (k in 1:length(sp_nm)) {
  for (f in 1:length(bof)) {
    # open suitability tiff file
    sp_rast<-paste0(ras, sp_nm[k],"_clipped_suitability_", bof[f], "_ROC_wmean.tif")
    sp_rast<-raster(sp_rast)
    
    # open vegetation overlay mask archipelago-wide for species 
    veg_mask<-paste0(vegDir, sp_nm[k],"_current_veg_mask.tif")
    veg_mask<-raster(veg_mask)
    
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
    resp_veg = sp_rast*veg_samp
    
    # convert all species range to binary values 
    resp_veg[resp_veg[] == 0]<-NA
    resp_veg[resp_veg[] > 0]<-1
    
    # save each raster for Kauai species and scenarios
    if (k == 1 && f == 1) {
      ake_bas<-resp_veg  # akekee baseline
    } else {
      if (k == 1 && f == 2) {
        ake_fut<-resp_veg  # akekee future
      } else {
        if (k == 2 && f == 1) {
          aki_bas<-resp_veg  # akikiki baseline
        } else {
          if (k == 2 && f == 2) {
            aki_fut<-resp_veg  # akikiki future
          }
        }
      }
    }
    # end of else statements to save endemic objects
  }
}

# stack Kauai endemic raster outputs
end_sp<-c("ake_bas", "ake_fut", "aki_bas", "aki_fut")
end_stack<-stack(ake_bas, ake_fut, aki_bas, aki_fut)

# open rasters for species and scenarios [set i = 1 and n = 1 for debugging]
for(i in 1:length(is_sp)) {
  for (n in 1:length(bof)){
    # open suitability ROC wmean
    sp_rast<-paste0(ras, is_sp[i],"_clipped_suitability_", bof[n], "_ROC_wmean.tif")
    sp_rast<-raster(sp_rast)
    
    # open vegetation overlay mask archipelago-wide for species 
    veg_mask<-paste0(vegDir, is_sp[i],"_current_veg_mask.tif")
    veg_mask<-raster(veg_mask)
    
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
    resp_veg = sp_rast*veg_samp
    
    # convert all species range to binary values 
    resp_veg[resp_veg[] == 0]<-NA
    resp_veg[resp_veg[] > 0]<-1
    
    for (j in 1:4){
      # sum rasters to get overlap area
      sp_sum<-sum(end_stack[[j]], resp_veg)
      
      # count frequency of presence pixels
      ol_pix<-freq(sp_sum)[1, 2]
      
      # divide by 4 (500 m) to get km^2
      ol_area<-sp_pix/4
      
      # create object with output 
      ol_info<-c(end_sp[[j]], is_sp[i], bof[n], ol_area)
      
      # save outputs of endemic comparisons 
      if (j == 1) {
        ol_dat<-ol_info
      } else {
        if (j > 1) {
          all_ol<-rbind(ol_dat , ol_info)
        }
      }
    }
        
    # save output for scenarios
    if (n == 1) {
      s_out<-all_ol
    } else {
      if (n > 1) {
        time_out<-rbind(s_out, all_ol)
      }
    }
  }
  # add output per species and save 
  if (i == 1) {
    ol_out<-time_out
  } else {
    if (i > 1) {
      all_out<-rbind(ol_out, time_out)
    }
  }
}

# calculate percent change from baseline to future
area<-as.numeric(all_out[ , 3])
change_1<-(area[1]-area[2])/area[1]*100
change_2<-(area[3]-area[4])/area[3]*100

# add percent change to matrix output
all_out<-cbind(all_out, c(NA, change_1, NA, change_2))

# add column headers to all_out object
colnames(all_out)<-c("species", "scenario", "area_km2", "%change")

# save output table
write.csv(all_out, file = paste0(outDir, is_run, "_", "endemics_habitat_area.csv"))
# repeat for other island runs [lines 19-21]
