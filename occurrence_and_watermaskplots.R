library(raster)
library(sf)
library(stars)
library(EBImage)

source("custom_fcns.R")

# dname <- 'Yukon'
dname <- 'Kobuk'
# dname = 'Nadym'
# dname = 'Ob'
# dname = 'Pur'
# dname <- "Mackenzie"
# dname <- "Yenisei"
# dname <- "Colville"
# dname <- 'Kolyma'
# dname <- "Lena"
# dname <- 'Yana'
# dname <- 'Indigirka'
#### Load in shoreline
shoreline <- st_read(paste0("deltas/",dname,"/outlines/",dname,"_shoreline_FULL.shp"))


### Load up unprojected tiffs from locally downloaded GSW data
dirname <- dir("deltas/", paste0(dname, "_GSW*"))

mo <- 7
pixres <- 30
fl <- list.files(paste0("deltas/",dirname,"/"),paste0("*0",mo,".tif$"))

# setwd("..")
# setwd("..")


setwd(paste0("deltas/",dname,"/",month.name[mo],"_clips"))

## run this FIRST to save clipped data to generate clipped occurrence in local UTM zone
lapply(fl, function(x) {
  # dirname <- dir("deltas/", paste0(dname, "_GSW*"))
  newm <- raster(paste0("deltas/",dirname,"/", x))

  mv <- cellStats(newm, max)
  if(mv == 2) {
    newm <- as(st_warp(st_as_stars(newm),  crs = st_crs(shoreline)$epsg, cellsize = pixres), "Raster")
    newm2 <- clip(newm, shoreline)
    writeRaster(newm2, x, format = "GTiff", options=c("COMPRESS=LZW", "TFW=YES"),
                datatype = "INT1U", overwrite = TRUE)
  }
})

masks <- list.files(paste0("deltas/",dname,"/",month.name[mo],"_clips"),paste0("*0",mo,".tif$"))

#### raster

ptm <- proc.time()
layers <- (stack(masks))

for(i in 1:length(masks)) {
  print(i)
  layers[[i]][layers[[i]]==0] <- NA
  layers[[i]] <- layers[[i]]-1
  
}

todrop <- which(as.numeric(substr(masks, 1, 4)) < 1999)
if(length(todrop)>0) {
  layers <- dropLayer(layers, todrop)
}
pixres <- res(layers)[1]
# layers <- round(layers)

######### temporary re-adding in the water 5 km farther
# Draw a 5 km buffer and fill in any water that was clipped out
# mostly only for Ob :|
if(dname == "Ob") {
  id <- match(c(2016:2018), substr(names(layers), 2, 5))
  # id <- id[!is.na(id)]
  for(i in id) {
    print(i)
    # plot(layers[[i]])
    ibuff2 <- (!t(as.matrix(layers[[i]])))*1
    # ibuff2 <- (!ibuff)*1
    ibuff2[is.na(ibuff2)] <- 0
    dmap <- EBImage::distmap((!ibuff2)*1)*pixres
    ibuff2[dmap>=5e3] <- 0
    nL <- raster(x = t(!ibuff2)*1, crs = crs(layers))
    extent(nL) <- extent(layers) #, crs = crs(layers), extent = extent(layers))
    nL2 <- mask(nL, shoreline)
    layers[[i]] <- nL2
    rm(nL2, ibuff2,dmap)
  }
  names(layers) <- paste0("x", substr(masks[-todrop],1,7))
}

## Computing occurrence for that month:
# Drop layers w/ poor collocation
badcoll <- list(Yukon = 0,
     Kobuk = 0,
     Nadym = 0,
     Ob = 0,
     Pur = 0,
     Mackenzie = 0,
     Yenisei = 0,
     Colville = 0,
     Kolyma = 0,
     Lena = 2016:2018,
     Yana = 2016:2018,
     Indigirka = 0)

setwd(paste0("deltas/",dname))
layers <- round(layers)
names(layers) <- paste0("x", substr(masks[-todrop],1,7))
if(length(todrop)==0) {
  names(layers) <- paste0("x", substr(masks, 1, 4))
}

proc.time() - ptm


nobs <- sapply(1:nlayers(layers), function(i) sum(layers[[i]][]>=0, na.rm=T))
names(nobs) <- as.numeric(substr(names(layers), 2, 5))
if(length(todrop)==0) {
  names(nobs) <- substr(masks, 1, 4)
}

Rp <- raster(matrix(1, nrow = nrow(layers), ncol = ncol(layers)), crs = crs(layers))
extent(Rp) <- extent(layers)
Npixtotal <- sum(!is.na(mask(Rp, shoreline)[]))

obs_frac <- nobs/max(Npixtotal)

toplot <- obs_frac >= 0.1

obs_frac_withverylow <- obs_frac
nobs <- nobs[toplot]
obs_frac <- obs_frac[toplot]
if(length(which(!toplot)) >= 1){
  layers <- round(dropLayer(layers, which(!toplot)))
  names(layers) <- names(obs_frac)
}


#### To get map of how many observations 
gc()
npixobs <- apply(as.array(dropLayer(layers, match(badcoll[[dname]], substr(names(layers), 2, 5)))),
                 c(1, 2), function(x) {sum(!is.na(x))})
npixobs_r <- raster((npixobs), crs = crs(layers))
extent(npixobs_r) <- extent(layers)

writeRaster(npixobs_r, paste0(month.name[mo], "_npix.tif"), format = "GTiff", options=c("COMPRESS=LZW", "TFW=NO"),
            overwrite = TRUE)


occurrence <- mean(dropLayer(layers, match(badcoll[[dname]], substr(names(layers), 2, 5))),
                   na.rm = T)
occurrence_m <- t(as.matrix(occurrence))
writeRaster(occurrence, paste0(month.name[mo], "_occurrence.tif"), format = "GTiff", options=c("COMPRESS=LZW", "TFW=NO"),
            overwrite = TRUE)


ind_bad <- obs_frac < .99
ind_bad_withverylow <- obs_frac_withverylow < .99

# w_frac: fraction of resolved pixels which are water
w_frac <- sapply(1:nlayers(layers), function(i) {
  mean(layers[[i]][], na.rm=T)
})

names(w_frac) <- substr(names(layers), 2, 5)

if(dname == "Lena") {
  ind_bad["2007"] <- F
} else if (dname == "Indigirka") {
  ind_bad["2007"] <- F
} else if (dname == "Yukon") {
  ind_bad[c("2016", "2017")] <- F
}
# fn <- list.files(paste0("deltas/",dname), "*channelmask.tif")
# id_yr <- substr(fn, 1, 4)


### means:
weight_mean <- sum(w_frac[]*obs_frac[]/sum(obs_frac[]), na.rm=T)
if(dname == 'Indigirka') {
  verywet <- match(2017, names(obs_frac))
  weight_mean <- sum(w_frac[-verywet]*obs_frac[-verywet]/sum(obs_frac[-verywet]), na.rm=T)
}

sort(abs(weight_mean - w_frac[!ind_bad]))

ystar <- list(Yukon = 2017,
              Pechora = 2013,
              Kobuk = 2012,
              Nadym = 1999,
              Ob = 2000,
              Pur = 2010,
              Taz = 2016,
              Mackenzie = 2004,
              Yenisei = 2016,
              Colville = 2011,
              Kolyma = 2009,
              Lena = 2013,
              Yana = 2011,
              Indigirka = 2016)

ystar2 <- list(Yukon = 2016,
               Pechora = 2009,
               Kobuk = 2002,
               Nadym = 2018,
               Ob = 2013,
               Pur = 2016,
               Taz = 2013,
               Mackenzie = 2014,
               Yenisei = 2013,
               Colville = 2008,
               Kolyma = 2013,
               Lena = 2007,
               Yana = 2014,
               Indigirka = 2007)


filename = paste0("deltas/",dname,"/plots/",month.name[mo],"_Wetness_Curve.png")
png(file = filename, width = 800, height = 400, pointsize = 14)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 2))

plot(names(w_frac[]), w_frac[],  xlim = c(1999, 2018),
     xlab = "Year", ylab = "Water Fraction", main = bquote(.(dname)~.(month.name[mo])~"Water Fraction"),
     cex = .5, cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.6, pch = 16)
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))

abline(h = weight_mean, lwd = 3, lty = 3)

## all years
points(names(w_frac)[ind_bad], w_frac[ind_bad],
       pch = 16, cex = 1.6, col = 'red')

points(names(w_frac)[!ind_bad], w_frac[!ind_bad],
       pch = 16, cex = 1.6, col = 'black')

points(ystar[[dname]], w_frac[match(ystar[[dname]], names(w_frac))],
       pch = 24, cex = 1.6, col = 'blue', bg = 'blue')

points(ystar2[[dname]], w_frac[match(ystar2[[dname]], names(w_frac))],
       pch = 25, cex = 1.6, col = 'blue', bg = 'blue')


points(badcoll[[dname]], w_frac[match(badcoll[[dname]], names(w_frac))],
       pch = 15, cex = 1.7, col = 'black')

badcoll_and_bad <- match(badcoll[[dname]], as.integer(names(which(ind_bad))))
if(length(badcoll_and_bad) > 0) {
  points(names(which(ind_bad)[badcoll_and_bad]), w_frac[which(ind_bad)[badcoll_and_bad]],
         pch = 15, cex = 1.7, col = 'red')
}


plot(names(obs_frac_withverylow), obs_frac_withverylow*100, xlab = "Year", ylab = "% Resolved",
     ylim = c(0, 100), pch = 16, cex = 1.6, cex.lab = 1.6, cex.main = 1.6,
     cex.axis = 1.6)
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
abline(h=99, lty = 2, lwd = 3)

points(names(obs_frac_withverylow)[!ind_bad_withverylow], obs_frac_withverylow[!ind_bad_withverylow]*100,
       pch = 16, cex = 1.6, col = 'black')

points(names(obs_frac_withverylow)[ind_bad_withverylow], obs_frac_withverylow[ind_bad_withverylow]*100,
       pch = 16, cex = 1.6, col = 'red')

points(ystar[[dname]], obs_frac_withverylow[match(ystar[[dname]], names(obs_frac_withverylow))]*100,
       pch = 24, cex = 1.7, col = 'blue', bg = 'blue')

points(ystar2[[dname]], obs_frac_withverylow[match(ystar2[[dname]], names(obs_frac_withverylow))]*100,
       pch = 25, cex = 1.7, col = 'blue', bg = 'blue')


points(badcoll[[dname]], obs_frac_withverylow[match(badcoll[[dname]], names(obs_frac_withverylow))]*100,
       pch = 15, cex = 1.7, col = 'black')

if(length(badcoll_and_bad) > 0) {
  points(names(which(ind_bad)[badcoll_and_bad]), obs_frac[which(ind_bad)[badcoll_and_bad]]*100,
         pch = 15, cex = 1.7, col = 'red')
}


dev.off()

plot(0, 0)

legend('center', legend = expression(y^"*",
  {y^"*"}["alt"],
  "High data quality",
  "Poor data quality",
  "Mis-collocated",
  "Mis-collocated and poor data quality"), pch = c(24, 25, 16, 16, 15, 15),
  col = c("blue", "blue", "black", "red", "black", "red"),
  pt.bg = c("blue", "blue", "black", "red", "black", "red"))






