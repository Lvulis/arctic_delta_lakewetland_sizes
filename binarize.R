## This script processes an individual GSW monthly water mask
## by filling in missing data and 
## Contact lvulis [at] uci.edu
## Last edited 2021-09-09

#### IMPORT LIBRARIES ####

library(raster)
library(EBImage)
library(reticulate)
library(stars)
library(sf)
ndimage <- import("scipy.ndimage")
measure <- import("skimage.measure")

source("custom_fcns.R")

#### Begin work ####
# Mask to binarize
tobinarize <- raster(paste0("deltas/",dname,"/tiffs_",Z,"/",Z, "_",mo,".tif"))
# Read in the subaerial delta (labelled as SHORELINE_FULL)
shoreline <- st_read(paste0("deltas/",dname,"/outlines/",dname,"_shoreline_FULL.shp"))

# Project to relevant UTM zone using stars
tobinarize_S <- st_as_stars(tobinarize)
tobinarize_S <- st_warp(tobinarize_S, crs = st_crs(shoreline)$epsg, cellsize = 30)
tobinarize <- as(tobinarize_S, "Raster")

# Use raster for iterative NA-filling (only really relevant within extent of delta)
# first make sure to get out only the subaerial delta
r2 <- tobinarize
r2[1:prod(dim(tobinarize))] <- 1
r2 <- clip(r2, shoreline)
studyzone <- which(r2[]>0)
if(length(which(clip(tobinarize, shoreline)[]==0))>0){
  stop = FALSE
} else {
  stop = T
}
NA_ind <- which(is.na(tobinarize[]))
tobinarize[NA_ind] <- -1
tobinarize[tobinarize==0] <- NA

# Then run the iterative filling

while(stop == FALSE) {
  tobinarize <- focal(tobinarize, w = matrix(1, 3, 3), fun = mean, pad = T, na.rm = T, NAonly= T)
  tobinarize2 <- clip(tobinarize, shoreline)
  if(sum(is.na(tobinarize2[studyzone])) >= 1) {
    stop <- FALSE
  } else {
    stop <- TRUE
  }
}


tobinarize <- tobinarize-1
tobinarize[NA_ind] <- NA
tobinarizer <- t(as.matrix(tobinarize))
tobinarize[] <- t(tobinarizer)

# "channelmask" includes lakes, channelmask_clean *does not*
fn <- paste0("deltas/",dname,"/",Z,"_",mo,"_channelmask.tif")

writeRaster(tobinarize, filename = fn,
            format = "GTiff", options=c("COMPRESS=LZW", "TFW=NO"),
            overwrite = TRUE, datatype = "INT1U")

#### Extract channel network ####

tobinarizer[NA_ind] <- 0

### keep largest basically
cncmp <- measure$label(tobinarizer, connectivity = 2)
cncmp[is.na(tobinarizer)] <- 0
mode(cncmp) <- 'integer'
props <- measureProps(cncmp)
tokp <- names(which.max(props[, 's.area']))
clean_mask <- rmObjects(cncmp, setdiff(rownames(props), tokp), reenumerate = F)
clean_mask[clean_mask>= 1] <- 1
clean_mask[NA_ind] <- NA
## Save the largest component for later
clean_maskr <- raster(t(clean_mask), crs = crs(tobinarize))
extent(clean_maskr) <- extent(tobinarize)

fn2 <- paste0("deltas/",dname,"/",Z,"_07_channel_clean.tif")

writeRaster(clean_maskr, filename = fn2,
            format = "GTiff", options=c("COMPRESS=LZW", "TFW=NO"),
            overwrite = TRUE, datatype = "INT1U")
clean_maskr[is.na(clean_maskr)] <- 0

writeRaster(clean_maskr, filename = paste0("deltas/",dname,"/",Z,"_07_channel_cleanNONA.tif"),
            format = "GTiff", options=c("COMPRESS=LZW", "TFW=NO"),
            overwrite = TRUE, datatype = "INT1U")

###########
rm(fn,fn2, clean_mask, clean_maskr, tokp, props, tobinarize,tobinarizer,
   r2, studyzone, cncmp, DCN)
