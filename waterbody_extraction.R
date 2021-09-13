### Runs waterbody extraction on an individual monthly mask

#### load relevant libraries ####
library(raster)
library(EBImage)
library(sf) 
library(stars)
library(reticulate)
source("custom_fcns.R")
measure <- import("skimage.measure")

#### Declare constants ####
## Constants or images that will be used to process, 
## we run each one at a time. 
## Need to have *done dname*, one *Z* and one *mo* uncommented.
###
# dname <- 'Yukon'
# Z     <- 2017
# Z     <- 2016
# mo    <- '07'

###
# dname <- 'Kobuk'
# Z     <- 2002
# Z     <- 2012
# mo    <- '07'

###
# dname <- 'Nadym'
# Z     <- 1999
# Z     <- 2018
# mo    <- '07'

###
# dname <- 'Ob'
# Z     <- 2000
# Z     <- 2013
# mo    <- '07'

###
# dname <- 'Pur'
# Z     <- 2010
# Z     <- 2016
# mo    <- '07'

###
# dname <- "Mackenzie"
# Z     <- 2004
# Z     <- 2014
# mo    <- '07'

###
# dname <- "Yenisei"
# Z     <- 2013
# Z     <- 2016
# mo    <- "07"

###
dname <- "Colville"
# Z     <- 2008
Z     <- 2011
mo    <- '07'
###
###
# dname <- 'Kolyma'
# Z     <- 2009
# Z     <- 2013
# mo    <- '07'

###
# dname <- "Lena"
# Z     <- 2013
# Z     <- 2007
# mo    <- '07'

###
# dname <- 'Yana'
# Z     <- 2011
# Z     <- 2014
# mo    <- '07'
###

# dname <- 'Indigirka'
# Z     <- 2007
# Z     <- 2016
# mo    <- '07'

## Size bins if want to do plots of log (area) from this script
bins <- seq(3.7, 8.1, by = 0.2)

#### START WITH binarize.R. to interpolate regions of the mask where necessary ####

# source("~/deltas/r_work/scripts/binarize.R")

## Load subaerial delta outline
shoreline <- st_read(paste0("G:/My Drive/Arctic/scripts_from_laptop/scripts/size_pdf_paper_scripts/",dname,"/outlines/",dname,"_shoreline_FULL.shp"))

#### load water mask (including channels)
water_mask <- raster(paste0("G:/My Drive/Arctic/scripts_from_laptop/scripts/size_pdf_paper_scripts/",dname, "/", Z, "_", mo, "_channelmask.tif"))

#### Load DCN 
DCN <- raster(paste0("G:/My Drive/Arctic/scripts_from_laptop/scripts/size_pdf_paper_scripts/", dname, "/", Z, "_", mo, "_channel_clean.tif"))
DCN <- clip(DCN, shoreline) # this is the equivalent of QGIS/ArcGIS clip, not base R clip
## Matrix form.
DCNt <- t(as.matrix(DCN))
# Pixel resolution is 30 m (Landsat) but we take it from the image
pixres = res(DCN)[1]
on_CN <- which(DCNt>0)

## Clip up and set up watermask for analysis
water_mask <- clip(water_mask, shoreline)
water_mask <- t(as.matrix(water_mask))
water_mask[water_mask<0] <- NA


# add a buffer around delta area so objects on edge don't get identical label

mask_hat_pre <- water_mask
mask_hat_pre[!is.na(mask_hat_pre)] <- 1
mask_hat_pre[is.na(mask_hat_pre)]  <- 0
kern <- makeBrush(3, shape = 'diamond')
mask_hat <- dilate(mask_hat_pre, kern)
# display(mask_hat - mask_hat_pre)
buffer_indices <- which((mask_hat - mask_hat_pre) == 1)
rm(mask_hat_pre, mask_hat)
gc()


## Removing pixels next to the DCN
water_mask_temp <- water_mask
water_mask_temp[buffer_indices] <- 0
water_mask_temp[on_CN] <- 0

# # CONNCOMPONETS USING RETICULATE

conncomp_july <- measure$label(water_mask_temp, connectivity = 2)
conncomp_july[is.na(water_mask_temp)] <- NA
conncomp_july[water_mask_temp == 0 ] <- NA
mode(conncomp_july) <- 'integer'

# Pick up object properties
lake_props <- measureProps(conncomp_july, pixres)

minpix <- 5
too_small <- rownames(lake_props[lake_props[, "s.area"] <= minpix*pixres^2, ])
conncomp_july <- rmObjects(conncomp_july, too_small, reenumerate = T)
mode(conncomp_july) <- 'integer'

rm(buffer_indices, lake_props)
gc()

#### Set NA to zeor for labelling the image 
conncomp_july[is.na(conncomp_july)] <- 0L
mode(conncomp_july) <- 'integer'
lake_props <- measureProps(conncomp_july, pixres)
obj_moments <- computeFeatures.moment(conncomp_july)

lake_props <- data.frame(lake_props)

### add in occurrence image (made separately using occurrence_and_watermaskplots.R)
occurrence <- raster(paste0("deltas/",dname,"/", month.name[as.integer(mo)], "_occurrence.tif"))
occurrence_m <- t(as.matrix(occurrence))

# Get xy indices of object on the collocated raster for occurrence computations
obj_inds <- splitObjects(conncomp_july)
obj_inds <- obj_inds[rownames(lake_props)]

# There is an na.rm because some of the occurrence map has 
# literally NO DATA (i.e. zero classifications) on the entire record
# but we interpolated the actual mask. Applies in every case
o_seasonality <- sapply(obj_inds, function(x) {
  mean(occurrence_m[x], na.rm = T)
})
names(o_seasonality) <- rownames(lake_props)

### Compute ephemerality around the lake for climate section
lake_map2t <- conncomp_july # [2996:3067, 1661:1719]

outer_peri <- diladd(lake_map2t, makeBrush(3, 'box'))
outer_peri[is.na(occurrence_m)] <- 0

outer_peri_inds <- splitObjects(outer_peri)
outer_peri_inds <- outer_peri_inds[rownames(lake_props)]
outer_peri_occ <- sapply(outer_peri_inds, function(x) {
  mean(occurrence_m[x], na.rm = T)
})

thresh <- 0.85
lake_props$perm <- o_seasonality
lake_props$perm[o_seasonality >= thresh] <- 1
lake_props$perm[o_seasonality < thresh]  <- 0
lake_props$o_seasonality <- o_seasonality
lake_props$s.area <- round(lake_props$s.area)
lake_props$EP1_occur <- outer_peri_occ
### EXPORT DATA
lake_props2 <- cbind(lake_props, obj_moments[rownames(lake_props), ])

fn = paste0("deltas/",dname,"/",Z,"_", month.name[as.integer(mo)], "_lake_properties.RData")
#### SAVE: 
# save(lake_props2, file=fn)

