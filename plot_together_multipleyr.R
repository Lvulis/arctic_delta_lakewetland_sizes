### Implement goodness of fit testing for competing distributions
library(EBImage)
library(sf)
library(poweRlaw)
library(truncdist)
library(parallel)
library(raster)
library(stars)

source("custom_fcns.R")

dname <- c("Yukon", "Kobuk", "Nadym", "Ob", "Pur",  "Mackenzie", "Yenisei", "Colville", "Kolyma", "Lena", "Yana", "Indigirka")

Z_mat <- matrix(c(2017, 2016, 2012, 2002, 1999, 2018, 2000, 2013, 2010, 2016, 2004, 2014, 2016, 2013, 2011, 2008, 2009, 2013, 2013, 2007, 2011, 2014, 2016, 2007), byrow = T, ncol = 2)

rownames(Z_mat) <- dname
pixres<-30

delta_props <- data.frame(read.csv("deltas/delta_properties.csv",header = T,
                                   stringsAsFactors = F, row.names = 1))

# delta_props <- delta_props[1:6, ]
delta_props <- delta_props[, dname]

MAAT <- as.numeric(delta_props["MAAT_degC", ])
names(MAAT) <- dname
# color_scale <- colorRampPalette(c("blue", "red"))
# color_scale_vals <- color_scale(length(MAAT))
# color_ranks <- round(rank(MAAT))
# names(color_ranks) <- names(MAAT)

LAT <- as.numeric(delta_props["Apex Latitude_degN", ])
names(LAT) <- dname
color_scale <- colorRampPalette(c("red", "blue"))
color_scale_vals <- color_scale(length(LAT))
color_ranks <- round(rank(LAT))
names(color_ranks) <- names(LAT)


arealimall <- c(3.7,8.1)
bins <- seq(arealimall[1], arealimall[2], by = 0.2)


P_ET <- as.numeric(delta_props["P_ET_JUNJUL", ])
names(P_ET) <- dname


shoreline <- lapply(dname, function(x) {
  shoreline <- st_read(paste0("deltas/",x,"/outlines/",x,"_shoreline_FULL.shp"))
})
names(shoreline) <- dname

##### Load up ice cover data from Brown et al. 1998
mp <- stars::read_stars("G:/My Drive/Arctic/imagery/ARPA/nhipa.byte", crs = 3408)
st_crs(mp) <- 3408

rcl_mat <- rbind(do.call(rbind, lapply(0:2, function(i) {
  cbind((1:4)+i*4, i+1)
})), cbind((1:4)+(3*4), 4),
cbind((1:4)+4*4, 5),
cbind(21, 6), # glacier
cbind(22, 7), # relict
cbind(c(23,24), 8), # Water
cbind(25, 9)) # land

mp_2 <- st_reclassify(mp, rcl_mat) # create low/med/high ice

shoreline_polar <- lapply(shoreline, st_transform, st_crs(3408)) # project shorelines into brown projection
# Extract map values for each delta
mpval <- lapply(shoreline_polar, function(x) {
  # extract(mp_2, x)
  mp[x]
})

# Summary table (counts)
mp_cnt <- lapply(mpval, function(x) table(x))
dmp <- as.integer(names(sapply(mpval, function(x) which.max(table(x)))))
# 
names(dmp) <- dname
iceval <- lapply(shoreline_polar, function(x) {
  # extract(mp_2, x)
  mp_2[x]
})


names(iceval) <- NULL

dice_cnt <- sapply(iceval, function(x) table(x))

dice <- as.integer(names(sapply(iceval, function(x) which.max(table(x)))))
names(dice) <- dname
icetable <- c("high", "medium", "low")
icetable[dice]

icemap_delta <- lapply(shoreline_polar, function(x) {
  # clip(mp_2, x)
  mp_2[x]
  # st_crop(mp_2, x)
})

names(icemap_delta) <- dname
#######

###### SIZE DATA
setwd("G:/My Drive/Arctic/scripts_from_laptop/scripts/size_pdf_paper_scripts")
obj_data <- lapply(dname, function(x) {
  setwd(x)
  fn = list.files(path=".", pattern = "*July_lake_properties.RData$", ignore.case = T)
  
  obj_dat1 <- lapply(fn, function(z) {
    load(z)
    return(lake_props2)
  })
  names(obj_dat1) <- substr(fn, 1, 4)
  
  setwd("..")
  return(obj_dat1)
 
})
names(obj_data) <- dname


### 
yr_list <- (sapply(obj_data, names))
mode(yr_list) <- 'integer'

cxv <- t(sapply(1:ncol(yr_list), function(i) {
  match(yr_list[, i], Z_mat[i, ])
})) 

rownames(cxv) <- dname
# if cxv = 1, y*. If cxv = 2, y* replicate.


######
obj_summaries <- lapply(obj_data, function(x) {

  sapply(x, function(y) {
    kk <- y[, 's.area']
    kk <- kk[kk>=6*pixres^2]
    return(c(mean = mean(kk),
             median = median(kk),
             var = var(kk),
             N = length(kk),
             skew = skewness(kk),
             kurt = kurtosis(kk),
             max = max(kk),
             A_total = sum(kk)))
  })
})



lapply(obj_summaries, function(x) x[, 1] > x[, 2])

## Just a wrapper to get the empirical ccdf from poweRlaw 
alldat <- lapply(obj_data, function(y) {
  lapply(y, function(x) {
    z <- displ$new(round(x[, 's.area']))
    plot(z, draw = F)
  })
})

waterbody_area_hists <- lapply(obj_data, function(y) {
  lapply(y, function(x) {
    
    histo <- hist(log10(round(x[, 's.area'])), breaks = bins, prob = T)
    histo$relcounts <- histo$counts/sum(histo$counts)
    histo
  })
})

### subset perennial lakes

theta <- 0.85
mperm_l <- lapply(obj_data, function(y) {
  mperm <- lapply(y, function(x) {
    
    o_seasonality3 <- x$o_seasonality
    o_seasonality3[o_seasonality3 >= theta] <- 1
    o_seasonality3[o_seasonality3  < theta] <- 0
    z <- x[(which(o_seasonality3 == 1)), 's.area']
    pp <- dislnorm$new(z/(pixres^2))
    pp$setXmin(6)
    pp$setPars(estimate_pars(pp))
    pp$internal$gof <- get_distance_statistic(pp, xmax = max(z/pixres^2)+1)
    pp$internal$AIC <- 2*length(pp$pars) - 2 * dist_ll(pp)
    pp
  })
  names(mperm) <- names(y)
  return(mperm)
})

mperm_l2 <- lapply(mperm_l, function(x) {
  lapply(x, function(y) {
    y2 <- dislnorm$new(y$dat*pixres^2)
    y2$setXmin(y$getXmin()*pixres^2)
    y2$setPars(y$getPars())
    y2$pars[1] <- y2$pars[1] + log(pixres^2)
    return(y2)
  })
})

names(mperm_l2) <- dname

permdat <- lapply(mperm_l2, function(y) {
  lapply(y, function(x) {

    plot(x, draw = F)
  })
})

lake_area_hists <- lapply(mperm_l, function(y) {

  lapply(y, function(xnew) {
    hist(log10(xnew$dat * (pixres^2)), breaks = bins, prob = T)
    
  })

})

perm_KS_val <- sapply(mperm_l, function(x) {
  sapply(x, function(y) {
    get_distance_statistic(y)
  })
})

colnames(perm_KS_val) <- dname

### KS_pval:
# 1: Calculate point estimates of LN parameters and their goodness of fit, the KS value (KSd) from N data points
# 2: For i in 1:M, 
#   a: simulate N values from LN distribution w/ fitted par
#   b: fit new parameters and KS statistic (KSsim)
#   c: If KSd > KSsim, P = P + 1
# 3: p = P/M


M <- 5000

no_cores <- detectCores() - 1
ptm <- proc.time()
cl <- makeCluster(no_cores)

clusterExport(cl, c("M", "pixres"))
clusterEvalQ(cl, library(poweRlaw))


perm_KS_list_lilliefors <- parLapply(cl, mperm_l, function(x) {

  KS_1 <- lapply(x, function(y) {
    # parS
    sapply(1:M, function(i) {
      print(i)
      newy <- dislnorm$new(dist_rand(y, y$internal$n))
      newy$setXmin(y$getXmin())
      newy$setPars(estimate_pars(newy))
      ptm <- proc.time()
      KS_est <- get_distance_statistic(newy, xmax = max(newy$dat)+1)
      proc.time() - ptm
      return(c(newy$pars, KS_est))
    })

  })
})

stopCluster(cl)

proc.time() - ptm


sapply(1:2, function(k) {
  mean(KS_1[[k]] >= x[[k]]$internal$gof)
})

perm_KS_pval_bootstrap <- sapply(seq_along(dname), function(i) {
  sapply(1:nrow(perm_KS_val), function(k) {
    mean(perm_KS_list_lilliefors[[i]][[k]][3, ] >= perm_KS_val[k, i])
    })
})

colnames(perm_KS_pval_bootstrap) <- dname


perm_KS_pval_bootstrap_reorg <- t(sapply(seq_along(dname), function(i) {
  t(perm_KS_pval_bootstrap)[i, cxv[i, ]]
}))

write.table(perm_KS_pval_bootstrap_reorg, "clipboard", sep="\t")
mperm_pl_l <- lapply(mperm_l, function(x) {
  lapply(x, function(y) {
    mp <- displ$new(y$dat)
    mp$setXmin(y$getXmin())
    mp$setPars(estimate_pars(mp))
    mp$internal$gof <- get_distance_statistic(mp, xmax = max(mp$dat)+1)
    mp$internal$AIC <- 2*length(mp$pars) - 2 * dist_ll(mp)
    
    return(mp)
  })
})

#### comparison of AICs for lakes

mperm_ln_AIC <- sapply(mperm_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
mperm_pl_AIC <- sapply(mperm_pl_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
# which is better
sign(mperm_pl_AIC - mperm_ln_AIC)

dAIC_lakes <- sapply(seq_along(dname), function(i) {
  mperm_pl_AIC[cxv[i, yrpick], i] - mperm_ln_AIC[cxv[i, yrpick], i]
})

write.table(dAIC_lakes, "clipboard", sep="\t", row.names =F, col.names = F)

###### split ephemeral data

meph_l <- lapply(seq_along(obj_data), function(i) {
  meph <- lapply(seq_len(ncol(cxv)), function(j) {
    x <- obj_data[[i]][[j]]
    thresh3 <- test_zone[min_id[j]]
    o_seasonality3 <- x$o_seasonality
    o_seasonality3[o_seasonality3 >= thresh3] <- 1
    o_seasonality3[o_seasonality3  < thresh3] <- 0
    x <- x[(which(o_seasonality3 == 0)), 's.area']
    x <- x[x>= (6 * pixres^2)]
    
    mp <- displ$new(x / (pixres^2))

    mp$setXmin(estimate_xmin(mp, xmax = max(mp$dat)+1))

    mp$internal$gof <- get_distance_statistic(mp, xmax = max(mp$dat)+1)
    mp$internal$AIC <- 2*length(mp$pars) - 2 * dist_ll(mp)
    return(mp)
  })
  names(meph) <- names(y)
  return(meph)
})
names(meph_l) <- dname

# # copy sample sizes
write.table(sapply(meph_l, function(x) sapply(x, function(y) y$internal$n)), "clipboard", sep="\t")

plot(MAAT, sapply(meph_l, function(x) sapply(x, function(y) y$internal$n))[1, ],
     ylab = 'n')

ephdat <- lapply(meph_l, function(y) {
  lapply(y, function(x) {
    x2 <- x$copy()
    x2$dat <- x2$dat * (pixres ^ 2)
    plot(x2, draw = F)
  })
})


wetland_area_hists <- lapply(meph_l, function(y) {
  lapply(y, function(xnew) {
    hist(log10(xnew$dat * (pixres ^ 2)), breaks = bins, prob = T)
    
  })

})


eph_KS_val <- sapply(meph_l, function(x) {
  sapply(x, function(y) {
    y$internal$gof
  })
})
print(eph_KS_val)


################
M <- 5000
# 
no_cores <- detectCores() - 1
ptm <- proc.time()
cl <- makeCluster(no_cores)
pars <- NULL
xmins <- NULL
clusterExport(cl, c("M", "pixres", "pars", "xmins", 'bootstrap_p_helper',
                    'sample_p_helper'))
clusterEvalQ(cl, library(poweRlaw))



eph_KS_bootstrap <- parLapply(cl, meph_l, function(x) {
  lapply(x, function(y) {
    m_cpy = y$copy()
    x = m_cpy$dat
    x_lower = x[x < m_cpy$xmin]
    xmax <- max(x) + 1
    
    ptm <- proc.time()
    bootstraps <- as.data.frame(t(sapply(1:M, bootstrap_p_helper, m_cpy, x_lower, xmins, pars, xmax, "ks")))
    proc.time() - ptm
    return(bootstraps)
    
  })
})

stopCluster(cl)
proc.time() - ptm


bootstraps <- bootstrap_p(y, no_of_sims = M, threads = 3)

eph_KS_pval_bootstrap <- sapply(seq_along(dname), function(i) {
  sapply(seq_along(meph_l[[i]]), function(k) {
    mean(eph_KS_bootstrap[[i]][[k]][, 1] >= eph_KS_val[k, i])
  })
})

###### procedure to simulate PL values w/ given slope/xmin
###### fit lnorm and compare

M <- 1000
# 
no_cores <- detectCores() - 2
ptm <- proc.time()
cl <- makeCluster(no_cores)

pars <- NULL
xmins <- NULL
clusterExport(cl, c("M", "pixres", "pars", "xmins", 'bootstrap_p_helper',
                    'sample_p_helper'))
clusterEvalQ(cl, library(poweRlaw))


eph_bootstrap_distributionfit_estxmin <- parLapply(cl, meph_l, function(x) {
  lapply(x, function(y) {
    n <- length(y$dat[y$dat>=y$xmin])
    xmax <- max(y$dat) + 1
    ptm <- proc.time()
    
    x_lower = y$dat[y$dat < y$xmin]
    
    
    bootstraps <- sapply(1:M, function(i) {
      ## only sample ABOVE xmin 
      ## sample and reestimate xmin as is done in the clauset paper
      newx <- sample_p_helper(i, y, x_lower)
      new_pl <- displ$new(newx)

      new_pl$setXmin(estimate_xmin(new_pl))
      new_pl$internal$gof <- get_distance_statistic(new_pl, xmax = xmax)
      new_pl$internal$AIC <- 2*length(new_pl$pars) - 2 * dist_ll(new_pl)
      l1  = c(new_pl$pars, dist_ll(new_pl), new_pl$internal$gof, new_pl$internal$AIC)
      gc()
      new_ln <- dislnorm$new(newx)
      new_ln$setXmin(new_pl$xmin)
      new_ln$setPars(estimate_pars(new_ln))
      new_ln$internal$gof <- get_distance_statistic(new_ln, xmax = xmax)
      new_ln$internal$AIC <- 2*length(new_ln$pars) - 2 * dist_ll(new_ln)
      l2 = c(new_ln$pars, dist_ll(new_ln), new_ln$internal$gof, new_ln$internal$AIC)
      vuongp <- compare_distributions(new_pl, new_ln)$p_one_sided
      rm(newx, new_pl, new_ln)
      gc()
      return(c(l1, l2, vuongp))
    })
    
    proc.time() - ptm
    return(bootstraps)
    
  })
})

stopCluster(cl)
proc.time() - ptm

# Pick one of the years to highlight in the supplementary
y <- eph_bootstrap_distributionfit_estxmin[[11]][[1]]

filename = paste0("deltas/plots/heavytailedcomparisons/plvsln_GOF.png")
png(file = filename, width = 800, height = 400, pointsize = 14)

par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 2)) # ORIGINAL FOR SQUARE

dAIC <- y[4, ] - y[9, ]
plot(density(dAIC), xlab = bquote(Delta*"AIC"), lwd = 2,
     cex.lab = 1.6,cex.axis=1.6,cex.main=1.6,
     main=bquote(bold(f[x]*"("*Delta*"AIC)")))

AIC_likelihood <- exp(-abs(dAIC)/2)
plot(density(AIC_likelihood), xlab = bquote("AIC p."), lwd = 2,
     cex.lab = 1.6,cex.axis=1.6,cex.main=1.6,
     main=bquote(bold(f[x]*"(AIC p.)")))

dev.off()

write.table(sumStat(dLL), "clipboard", sep="\t")
write.table(sumStat(dAIC), "clipboard", sep="\t")
write.table(sumStat(AIC_likelihood ), "clipboard", sep="\t")
write.table(sumStat(y[10, ]), "clipboard", sep="\t")

### fit lognorm to ephemeral waterbodies
meph_lnorm_l <- lapply(meph_l, function(x) {
  lapply(x, function(y) {
    zz <- dislnorm$new(y$dat)

    zz$setXmin(y$getXmin())
    zz$setPars(estimate_pars(zz))
    zz$internal$gof <- get_distance_statistic(zz, xmax = max(zz$dat)+1)
    zz$internal$AIC <- 2*length(zz$pars) - 2 * dist_ll(zz)
    
    return(zz)
  })
})


#### comparison of AICs:

meph_ln_AIC <- sapply(meph_lnorm_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
meph_pl_AIC <- sapply(meph_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
# which is better
sign(meph_pl_AIC - meph_ln_AIC)
# p val:

#### goodness of fit of the lnorm distributions fit to wetland sizes
M <- 5e3

no_cores <- detectCores() - 1
ptm <- proc.time()
cl <- makeCluster(no_cores)

clusterExport(cl, c("M", "pixres", "meph_lnorm_l"))
clusterEvalQ(cl, library(poweRlaw))


### If the lognorm mean in ln scale is less than -10 this doesn't converge
### be weary of anything with meanlog < -10 !
eph_lnorm_KS_list_lilliefors <- parLapply(cl, meph_lnorm_l, function(x) {

  KS_1 <- lapply(x, function(y) {
    if(y$pars[1] > - 10){
    sapply(1:M, function(i) {
      print(i)
      newy <- dislnorm$new(dist_rand(y, y$internal$n))
      newy$setXmin(y$getXmin())
      newy$setPars(estimate_pars(newy))
      ptm <- proc.time()
      KS_est <- get_distance_statistic(newy, xmax = max(newy$dat)+1)
      proc.time() - ptm
      return(c(newy$pars, KS_est))
    })
    } else {
      return(NULL)
    }

  })
})

stopCluster(cl)

proc.time() - ptm


ephlnorm_KS_val <- sapply(meph_lnorm_l, function(x) {
  sapply(x, function(y) {
    y$internal$gof
  })
})

ephlnorm_KS_pval_bootstrap <- sapply(seq_along(dname), function(i) {
  sapply(seq_along(meph_lnorm_l[[i]]), function(k) {
    if(!is.null(eph_lnorm_KS_list_lilliefors[[i]][[k]])) {
      mean(eph_lnorm_KS_list_lilliefors[[i]][[k]][3, ] >= ephlnorm_KS_val[k, i], na.rm = T)
    } else {
      NULL
    }

  })
})

#### PICK YEAR:

## pars to take down: 
# xmin
# N
# xmax (already done)
# PL exponent
# PL KS
# LN mean
# LN sd
# LN KS
yrpick<-1
wetlandpars <- t(sapply(seq_along(dname), function(i) {
  il <- cxv[i, yrpick]
  yr <- yr_list[il, i]
  pl_xmin_pix <- meph_l[[i]][[il]]$getXmin()
  pl_xmin <- meph_l[[i]][[il]]$getXmin() * (pixres ^ 2) / 1e5
  pl_xmax_pix <- max(meph_l[[i]][[il]]$dat)
  pl_xmax <- max(meph_l[[i]][[il]]$dat) * (pixres ^ 2) / 1e5
  N_wetland <- meph_l[[i]][[il]]$internal$n
  plpars <- meph_l[[i]][[il]]$pars
  pl_KS <- eph_KS_val[il, i]
  pl_AIC <- meph_l[[i]][[il]]$internal$AIC
  ln_pars <- meph_lnorm_l[[i]][[il]]$pars
  ln_pars[1] <- ln_pars[1] + log(pixres ^ 2)
  ln_pars <- log10(exp(ln_pars))
  ln_KS <- meph_lnorm_l[[i]][[il]]$internal$gof
  ln_AIC <- meph_lnorm_l[[i]][[il]]$internal$AIC
  AIC_comp <- exp(-abs((pl_AIC - ln_AIC))/2)

  frac <- fractal_D[il, i]
  
  return(c(yr, pl_xmin_pix, pl_xmin, pl_xmax_pix, pl_xmax,
           N_wetland, plpars, pl_KS, pl_AIC, ln_pars,
           ln_KS, ln_AIC, AIC_comp, frac))
}))


rownames(wetlandpars) <- dname
colnames(wetlandpars) <- c("Year", "xmin_pix", 
                           "xmin_10^5m2", "xmax_pix",
                           "xmax_10^5m2",
                           "Ntail", "powlaw_exp", 
                           "powlaw_KS",
                           "powlaw_AIC",
                           "lognorm_mu", "lognorm_sd", 
                           "lognorm_KS", "lognorm_AIC",
                           "AIC_pval")

write.csv(wetlandpars, paste0("deltas/ephvper/wetland_pars.csv"))


###### lognormal and power law for all waterbodies ##############

mall_l <- lapply(obj_data, function(y) lapply(y, function(x) {

  mall <- dislnorm$new(round(x[, 's.area'])/(pixres^2))
  mall$setPars(estimate_pars(mall))
  mall$internal$gof <- get_distance_statistic(mall, xmax = max(mall$dat)+1)
  mall$internal$AIC <-  2*length(mall$pars) - 2 * dist_ll(mall)
  return(mall)
}))

mall_pl_l <- lapply(obj_data, function(y) lapply(y, function(x) {
  mall <- displ$new(round(x[, 's.area'])/(pixres^2))
  mall$setXmin(6)
  mall$setPars(estimate_pars(mall))
  mall$internal$gof <- get_distance_statistic(mall, xmax = max(mall$dat)+1)
  mall$internal$AIC <-  2*length(mall$pars) - 2 * dist_ll(mall)
  return(mall)
}))

all_KS_val <- sapply(mall_l, function(x) {
  sapply(x, function(y) {
    y$internal$gof
  })
})

colnames(all_KS_val) <- dname


### KS_pval:
# 1: Calculate point estimates of LN parameters and their goodness of fit, the KS value (KSd) from N data points
# 2: For i in 1:M, 
#   a: simulate N values from LN distribution w/ fitted par
#   b: fit new parameters and KS statistic (KSsim)
#   c: If KSd > KSsim, P = P + 1
# 3: p = P/M

M <- 5000

# no_cores <- detectCores() - 1 
no_cores <- 3
ptm <- proc.time()
cl <- makeCluster(no_cores)
clusterExport(cl, c("M", "pixres"))
clusterEvalQ(cl, library(poweRlaw))


all_KS_list <- parLapply(cl, mall_l, function(x) {

  KS_1 <- lapply(x, function(y) {
    # parS
    sapply(1:M, function(i) {
      newy <- dislnorm$new(dist_rand(y, y$internal$n))
      newy$setXmin(y$getXmin())
      # newy$setPars(y$getPars())
      newy$setPars(estimate_pars(newy))
      ptm <- proc.time()
      KS_est <- get_distance_statistic(newy, xmax = max(newy$dat)+1)
      proc.time() - ptm
      return(KS_est)
    })
    
    
  })
})

stopCluster(cl)

proc.time() - ptm

all_KS_pval_bootstrap <- sapply(seq_along(dname), function(i) {
  sapply(1:nrow(all_KS_val), function(k) {
    mean(all_KS_list[[i]][[k]] >= all_KS_val[k, i])
  })
})


colnames(all_KS_pval_bootstrap) <- dname

all_KS_pval_bootstrap_reorg <- t(sapply(seq_along(dname), function(i) {
  t(all_KS_pval_bootstrap)[i, cxv[i, ]]
}))

write.table(all_KS_pval_bootstrap_reorg, "clipboard", sep="\t")

#### comparison of AICs for lakes

mall_ln_AIC <- sapply(mall_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
mall_pl_AIC <- sapply(mall_pl_l, function(x) sapply(x, function(y) {
  y$internal$AIC
}))
# which is better
sign(mall_pl_AIC - mall_ln_AIC)

dAIC_waterbodies <- sapply(seq_along(dname), function(i) {
  mall_pl_AIC[cxv[i, yrpick], i] - mall_ln_AIC[cxv[i, yrpick], i]
})
write.table(dAIC_waterbodies, "clipboard", sep="\t", row.names =F, col.names = F)


#########

mallpars <- sapply(mall_l, function(y) {
  pars = sapply(y, function(x) {
    z = x$pars
    z[1] = z[1]+log(pixres^2)
    z
  })
})


mpermpars <- sapply(mperm_l, function(y) {
  pars = sapply(y, function(x) {
    z = x$pars
    z[1] = z[1]+log(pixres^2)
    z
  })
})

mephpars <- sapply(meph_l, function(x) sapply(x, function(y) {
  y$pars
}))

yrpick <- 1
MAATr <- range(MAAT)
MAATr <- MAATr * c(1.1, 0.9)  

mpermxmean <- sapply(mperm_l, function(y) {
  pars = sapply(y, function(x) mean(x$internal$dat[x$internal$dat>=x$getXmin()] * (pixres^2)))
})


#### PICK YEAR:
yrpick <- 1
###
exportpars <- t(sapply(seq_along(dname), function(i) {
  il <- cxv[i, yrpick]
  yr <- yr_list[il, i]
  N_lake <- mperm_l[[i]][[il]]$internal$n
  lnpars <- mperm_l[[i]][[il]]$pars
  lnpars[1] <- lnpars[1] + log(pixres ^ 2)
  lnpars <- log10(exp(lnpars))
  ln_KS <- perm_KS_val[il, i]
  lnp <- perm_KS_pval_bootstrap[il, i]
  N_wetland <- meph_l[[i]][[il]]$internal$n
  pl_xmin <- meph_l[[i]][[il]]$getXmin() * (pixres ^ 2) / 1e5
  pl_xmax <- max(meph_l[[i]][[il]]$dat) * (pixres ^ 2) / 1e5
  
  plpars <- meph_l[[i]][[il]]$pars
  
  pl_KS <- eph_KS_val[il, i]
  plp <- eph_KS_pval_bootstrap[il, i]
  return(c(yr, N_lake, lnpars, lnp, N_wetland, pl_xmin, pl_xmax, plpars,plp))
}))


rownames(exportpars) <- dname
colnames(exportpars) <- c("Year", "N_lake", "lake_log_mu", "lake_log_sd", "lake_lnorm_p",
                          "Ntail_wetland", "wetland_powlaw_xmin_10^5m2","wetland_powlaw_xmax_10^5m2",
                          "wetland_powlaw_exp", "wetland_powlaw_p")

write.csv(exportpars, paste0("deltas/ephvper/seperate_distr_pars.csv"))


exportpars_all <- t(sapply(seq_along(dname), function(i) {
  il <- cxv[i, yrpick]
  yr <- yr_list[il, i]
  
  N_lake <- mall_l[[i]][[il]]$internal$n
  lnpars <- mall_l[[i]][[il]]$pars
  lnpars[1] <- lnpars[1] + log(pixres ^ 2)
  lnpars <- log10(exp(lnpars))
  ln_KS <- all_KS_val[il, i]
  lnp <- all_KS_pval_bootstrap[il, i]
 
  return(c(yr, N_lake, lnpars, lnp))
}))


rownames(exportpars_all) <- dname
colnames(exportpars_all) <- c("Year", "N_waterbody", "waterbody_log_mu", "waterbody_log_sd", "waterbody_lnorm_p")
write.csv(exportpars_all, paste0("deltas/ephvper/all_distr_pars.csv"))


yrpick <- 1
# perm <- 3
# use width = 300 height = 700 for the wetland CCDF
for(perm in c(1,3,2)) {
    filename = paste0("deltas/plots/",plotname[perm],"_area_ccdf.png")
    png(file = filename, width = 300, height = 700, pointsize = 14)
    # png(file = filename, width = 500, height = 500, pointsize = 14)
    
    # par(mar = c(5.1, 5.1, 4.1, 2.1)) # ORIGINAL FOR SQUARE
    if(perm == 2){
      par(mar = c(5.1, 5.1, 2.1, 2.1)) # Tall
      xlim <- c(arealimall[1], 7) #7 for wetlands, 8 for lakes
      
    } else if (perm == 3) {
      par(mar = c(5.1, 2.1, 2.1, 5.1)) # Tall
      xlim <- c(arealimall[1], 8) #7 for wetlands, 8 for lakes
      
    }
    
    # xlim = arealimall
    # xlims_1 <- pretty(xlim)
    # ylims_1 <- 10^pretty(c(-5, 0), by = 1) # ylims -5 wihout the scaling thing!
    
    
    xlim <- c(arealimall[1], 8) #7 for wetlands, 8 for lakes
    xlims_1 <- pretty(xlim, n = 4)
    ylims_1 <- 10^seq(-12, 0, by = 3)

    # xlab = bquote("log"[10]*"(Area [m"^2*"])"),
    #main =bquote(bold("July"~"1-"~F(A[L])))
    plot(0, xlim = 10^xlim,  ylab = "",
         main = "", xlab = "", type = "l", 
         lwd = 2, ylim = range(ylims_1), log = "xy",
         xaxt = "n", yaxt = "n", cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6,
         )
    abline(h=ylims_1, lty = 3, lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
    abline(v=10^xlims_1, lty = 3, lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
    
    if(perm <=2) {
      axis(2, at=ylims_1, lab=sapply(ylims_1, function(x) {
        as.expression(bquote(10^.(log10(as.numeric(x)))))
      }), cex.axis = 1.8, cex.lab = 1.6, las = 1)
      title(ylab = bquote({P^"*"}(A > a)), line = 3.5, cex.lab = 1.6, cex.main = 1.6)
      
    } else if(perm == 3) {
      axis(4, at=ylims_1, lab=sapply(ylims_1, function(x) {
        as.expression(bquote(10^.(log10(as.numeric(x)))))
      }), cex.axis = 1.8, cex.lab = 1.6, las = 1)
      mtext(bquote({P}(A > a)), side = 4, line = 4.0, cex = 1.6)
      
    }
    
    axis(1, at=10^xlims_1, lab=sapply(xlims_1, function(x) {
      as.expression(bquote(.(as.numeric(x))))
    }), cex.axis = 1.8, cex.lab = 1.6)
    title(xlab = bquote("log"[10]*"(Area [m"^2*"])"), line = 3.3, cex.main = 1.6, cex.lab = 1.6)
    
    ccdfscaler <- 10^-seq(8, 0, length.out = 12)
    plrange <- sapply(seq_along(meph_l), function(i) {
      meph_l[[dname[i]]][[cxv[i, yrpick]]]$pars
    })
    plsort <- rank(-plrange)
    ccdfscaler <- ccdfscaler[plsort]
    
  if(perm == 1) {
    lapply(seq_along(dname), function(i) {
      lines(alldat[[i]][[cxv[i, yrpick]]][, 1],
            alldat[[i]][[cxv[i, yrpick]]][, 2], lwd = 2, type = "l",
            col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
    })
  } else if (perm == 2) {
    boundind <- c(which.min(plrange), which.max(plrange))
    iterator <- seq_along(dname)
    lapply(iterator, function(i) {
      xmin_fit <- meph_l[[i]][[cxv[i, yrpick]]]$getXmin() * (pixres^2)
      xd <- ephdat[[i]][[cxv[i, yrpick]]]
      above_xmin <- which(xd[, 1] >= xmin_fit)
      points(xd[-above_xmin, 1],
             xd[-above_xmin, 2]*ccdfscaler[i], lwd = 2, type = "p",
             col = 'grey', pch = i, cex = 0.25, lty = 1)
      
      points(xd[above_xmin, 1],
             xd[above_xmin, 2]*ccdfscaler[i], lwd = 2, type = "p",
             col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
     
      
    })
   
  } else if (perm == 3) {
    
    boundind <- c(which.min(plrange), which.max(plrange))
    iterator <- seq_along(dname)
    lapply(iterator, function(i) {
      xmin_fit <- mperm_l[[i]][[cxv[i, yrpick]]]$getXmin() * (pixres^2)
      # i <- plsort[j]
      xd <- permdat[[i]][[cxv[i, yrpick]]]
      above_xmin <- which(xd[, 1] >= xmin_fit)
      # if(i %in% boundind){
      points(xd[-above_xmin, 1],
             xd[-above_xmin, 2]*ccdfscaler[i], lwd = 2, type = "p",
             col = 'grey', pch = i, cex = 0.25, lty = 1)
      
      # col = color_scale_vals[color_ranks[dname[i]]]
      # col = rgb(.25, .25, .25)
      points(xd[above_xmin, 1],
             xd[above_xmin, 2]*ccdfscaler[i], lwd = 2, type = "p",
             col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
    })
    # ORIGINAL on squarish, all on top of one another
    # lapply(seq_along(dname), function(i) {
    #   lines(permdat[[i]][[cxv[i, yrpick]]][, 1],
    #         permdat[[i]][[cxv[i, yrpick]]][, 2], lwd = 2, type = "p",
    #         col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
    
    # })
  } else {
    stop()
  }
   
  
  # lapply(seq_along(dname), function(i) {
    # lines(alldat[[i]][[cxv[i, yrpick]]], lwd = 2, type = "l",
    #       col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
    # lines(mperm_l2[[i]][[cxv_max[[i]]]], lwd = 4, type = "l",
    #       col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 2)
    # lines(permdat[[i]][[cxv[i, yrpick]]], lwd = 2, type = "l",
    #       col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
    # lines(ephdat[[i]][[cxv[i, yrpick]]], lwd = 2, type = "l",
    #       col = color_scale_vals[color_ranks[dname[i]]], pch = i, cex = 0.5, lty = 1)
  # })
  dev.off()
}
# 
# 

plot_xseq <- seq(log10(6*pixres^2),arealimall[2],by=0.01)

### KEEP IN MIND 500 x 300 for FIGURE S3/S5 in paper and 500 x 500 for S^
for(perm in c(1, 2, 3)) {
  filename = paste0("deltas/plots/",plotname[perm],"_area_pdf.png")
  png(file = filename, width = 500, height = 500, pointsize = 14)
    
  xlim = arealimall
  xlims_1 <- pretty(xlim)
  
  # Can mess with this to get better ylim range if actually plotting wetland/waterbody cdf
  yr <- unlist(lapply(seq_along(dname), function(i) lake_area_hists[[i]][[cxv[i, yrpick]]]$density))

  ylims <- range(yr[yr>(.Machine$double.eps*10)])
  ylims[2] <- ylims[2]*1.1

  par(mar = c(5.1, 5.1, 4.1, 2.1))
  plot(999,999, xlim = xlim, ylab = "",
        type = "l", xlab = "",
       lwd = 2, ylim = ylims, log = "",
       xaxt = "n",  cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6)
  grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
  
  axis(1, at=xlims_1, lab=sapply(xlims_1, function(x) {
    as.expression(bquote(.(as.numeric(x))))
  }), cex.axis = 1.8, cex.lab = 1.6)
  title(xlab = bquote("log"[10]*"(Area [m"^2*"])"), line = 3.3, cex.main = 1.6, cex.lab = 1.6)
  title(ylab = bquote(bold("Density")), line = 3.5, cex.lab = 1.6)
  
  if(perm == 1) {
    lapply(seq_along(dname), function(i) {
      lines(waterbody_area_hists[[i]][[cxv[i, yrpick]]]$mids, waterbody_area_hists[[i]][[cxv[i, yrpick]]]$density,
            lwd = 2, type = "l", lty = 1,
            col = color_scale_vals[color_ranks[dname[i]]], pch = i )
    })
  } else if (perm == 2) {
    lapply(seq_along(dname), function(i) {
      lines(wetland_area_hists[[i]][[cxv[i, yrpick]]]$mids, wetland_area_hists[[i]][[cxv[i, yrpick]]]$density,
            lwd = 2, type = "p", lty = 1,
            col = color_scale_vals[color_ranks[dname[i]]], pch = i )
    })
  } else if (perm == 3) {
    lapply(seq_along(dname), function(i) {
      lines(lake_area_hists[[i]][[cxv[i, yrpick]]]$mids, lake_area_hists[[i]][[cxv[i, yrpick]]]$density,
            lwd = 2, type = "l", lty = 1,
            col = color_scale_vals[color_ranks[dname[i]]], pch = i )
    })
  } else {
    stop()
  }
  
  dev.off()
}



### For making FIG 2 panel B/C polygon


dn <- 'Kolyma'
i <- match(dn, dname)

lp_t <- subset(obj_data[[dn]][[cxv[dn, yrpick]]], perm >= 0)

perm_frac <- sapply(1:(length(bins)-1), function(k) {
  o_id <- which((lp_t$s.area > 10^bins[k]) & (lp_t$s.area <= 10^bins[k+1]))
  # counts <- length(o_id)
  
  sum(lp_t$perm[o_id])/length(o_id)
})

col1 <- rgb(0 , 0.6, 0, 0.55)
col2 <- rgb(0, 0, 1, 0.55)


filename = paste0("deltas/plots/lake_wetland_pdf_partition.png")
png(file = filename, width = 500, height = 500, pointsize = 14)
par(mar = c(5.1, 5.1, 4.1, 2.1))

xlim = arealimall
xlims_1 <- pretty(xlim)

ylims <- c(0, 0.2)
ylims_1 <- pretty((ylims))


plot(999,999, xlim = xlim, xlab = bquote("log"[10]*"(Area [m"^2*"])"), ylab = "",
     main = bquote(bold(PDF~.(dn)~.(Z_mat[dn, yrpick]))), type = "l",
     lwd = 2, ylim = ylims, log = "",
     cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6)
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
title(ylab = bquote(bold("Frequency")), line = 3.5, cex.lab = 1.6)

lines(waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$mids, waterbody_area_hists[[dn]][[cxv[i, yrpick]]]$relcounts,
      lwd = 2, type = "l", lty = 1,
      col = 'black', pch = i )
relcounts_permfrac <- waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$relcounts*perm_frac
relcounts_permfrac[is.nan(relcounts_permfrac)] <- ylims[1]
relcounts_toplot <- waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$relcounts
relcounts_toplot[relcounts_toplot==0] <- ylims[1]
polygon(c(waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$mids, rev(waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$mids)),
        c(relcounts_permfrac, rev(rep(ylims[1], length(relcounts_permfrac)))),
        col = col1, border = NA)
polygon(c(waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$mids, rev(waterbody_area_hists[[dn]][[cxv[dn, yrpick]]]$mids)),
        c(relcounts_permfrac, rev(relcounts_toplot)),
        col = col2, border = NA)



dev.off()



#######
yrpick <- 2
znew_l <- lapply(0:2, function(perm) {
  znew <- matrix(0, nrow = length(dname), ncol = 10)
  
  
  for(i in 1:nrow(cxv)) {
    print(i)
    ### all
    if(perm == 0) {
      q <- obj_data[[i]][[cxv[i, yrpick]]][, 's.area']
      q <- q[q>=mperm_l[[i]][[cxv[i, yrpick]]]$getXmin() * (pixres^2)]

    } else if (perm == 1) {
      q <- meph_l[[i]][[cxv[i, yrpick]]]$dat * (pixres^2)
     
      ### ephemeral
    } else if (perm == 2) {
      q <- mperm_l[[i]][[cxv[i, yrpick]]]$dat[mperm_l[[i]][[cxv[i, yrpick]]]$dat>=mperm_l[[i]][[cxv[i, yrpick]]]$getXmin()] * (pixres^2)
      ### perennial

    }
    znew[i, 1] <- quantile((q), 0.25)
    znew[i, 2] <- quantile((q), 0.9)
    znew[i, 3] <- log10(mean((q)))
    znew[i, 4] <- log10(sd((q)))
    znew[i, 5] <- length(q)
  }
  znew <- as.data.frame(znew)
  rownames(znew) <- dname
  colnames(znew) <- c("25", "90", "mean", "sd", "N")
  return(znew)
})

#### mean area versus MAAT


ylim <- list(ystar = c(1.5, 25),
             ystar_alt = c(1.5, 25))

yr_set <- 'ystar'
plotname <- c("Waterbody", "Wetland", "Lake")


# perm <- 3
for(perm in c(1,3,2)) {
  
  filename = paste0("deltas/plots/",plotname[perm],"mean_MAAT.png")
  png(file = filename, width = 500, height = 500, pointsize = 14)
  par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,1))
  
  realdat <- data.frame(MAAT = as.numeric(MAAT[order(MAAT)])[], mu = (10^znew_l[[perm]]$`mean`[order(MAAT)])[])
  
  corspearman <- cor(realdat, method = 'spearman')[2]
  print(paste0("Spearman tau is: ", corspearman))              

    plot(as.numeric(MAAT), 10^znew_l[[perm]]$`mean`/1e4, 
       ylab = bquote("Mean"~.(plotname[perm])~"Area [10"^4~"m"^2*"]"),
       xlab = bquote("Mean Annual Air Temperature ["*degree*"C]"),
       cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.6, 
       pch = dice[]+14, cex = 3,
       main = bquote(bar(A)), xlim = rev(MAATr),
       ylim = ylim[[yr_set]],
       las = 1,
       col = color_scale_vals[color_ranks[dname]],
       panel.first = grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5)))
  
  
  ## trend line setup
  mod <- lm(mu ~ MAAT, data = realdat)
  lines(realdat$MAAT, ((mod$fitted.values))/1e4,  lwd = 3, lty = 2, col = 'black')
  
  ### Bootstrap the GOF
  slope_est <- sapply(1:1e4, function(i) {
    ## Randomly shuffle the Y data against X,
    rd <- data.frame(mu = sample(realdat$mu), MAAT = sample(realdat$MAAT))
    c(lm(mu ~ MAAT, data = rd)$coefficients, cor(rd, method = 'spearman')[2])
  })
  
  ## Shuffling y:
  # get pval of spearman 
  cp = mean((slope_est[3, ] > abs(corspearman)) | (slope_est[3, ] < -abs(corspearman)))
  print(paste0("bootstrap spearman p-val is: ", cp))
  
  # pval of regression
  qp = mean((slope_est[2, ] > abs(mod$coefficients[2])) | (slope_est[2, ] < -abs(mod$coefficients[2])))
  print(paste0("bootstrap p-val is: ", qp))
  
  ### plot regression
  coef1 <- round(mod$coefficients[1])
  coef2 <- round(mod$coefficients[2])
  r2 <- round(summary(mod)$r.squared*100)/100
  pval <- round(summary(mod)$coefficients[2, 4]*1e4)/1e4

  print(paste0("t-test p-val is: ", pval))
  if(pval>0.05) {
    legtrac <- legend("topleft", legend = bquote("R"^2 == .(r2)*"," ~ p == .(pval) ),
                      box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
                      cex = 1.6,trace = T, plot = T)
  } else {
    legtrac <- legend("topleft",
                      legend = bquote(atop(bar(Area) %prop% .(coef2)*"MAAT",
                                           "R"^2 == .(r2)*"," ~ p == .(pval) )),
                      box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
                      cex = 1.6,trace = F, plot = T)
  }
  ###
  
  
  legend("topright",legend = c("", "Low", "Medium", "High"), cex = 1.6, pch = c(NA, 14+(3:1)),
         col = 'black')
  dev.off()
}


### PLOT OF WETLAND EXPONENTS SORTED BY EXPONENT
yrpick<-2

filename = paste0("deltas/plots/wetlandPL_ranked.png")
png(file = filename, width = 200, height = 700, pointsize = 14)
par(mar = c(5.1, 0.5, 2.1, 0.55))

plot(-1, -1, ylab = '', yaxt = 'n',
     xlab = bquote(alpha), xlim = c(1.5, 3.5),
     ylim = c(1, length(dname)), cex.axis = 1.6,
     cex.main = 1.6, cex.lab = 1.6,
    font.lab = 2)
# panel.first = grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))

plrange <- sapply(seq_along(meph_l), function(i) {
  meph_l[[dname[i]]][[cxv[i, yrpick]]]$pars
})
names(plrange) <- dname
plsort <- rank(-plrange)

dn_topl <- names(plrange)
longname <- nchar(dn_topl) > 4

dn_topl <- substr(dn_topl, 1, 4)
dn_topl[longname] <- paste0(dn_topl[longname], ".")

text(plrange+0.35, plsort, dn_topl, 
     cex = 1, font = 2)
axis(4, at=seq_along(dname), lab = rep("", length(dname)),
     # lab = names(sort(plrange, decreasing = T)),
     cex.axis = 1, cex.lab = 1, las = 1, font.axis = 2, tck = 0.05)

points( plrange, plsort,
        pch = dice+14, cex = 2, col = color_scale_vals[color_ranks])

dev.off()

############ BOXPLOT OR OTHER SUMMARY STATISTICS ##########

yrpick <- 1
MAATlab <- pretty(rev(MAATr))
yr_set <- 'ystar'

# perm <- 3
for(perm in 1:3) {
  filename = paste0("deltas/plots/",plotname[perm],"_area_boxplot_MAAT.png")
  png(file = filename, width = 575, height = 575, pointsize = 14)
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  
  
  #perm = 0 is all, perm = 1 wetlands, perm = 2 is lakes
  if(perm == 3) {
    LL <- lapply(seq_along(dname), function(i) log10(mperm_l[[i]][[cxv[i, yrpick]]]$dat * pixres^2))
  } else if(perm == 2) {
    LL <- lapply(seq_along(dname), function(i) log10(meph_l[[i]][[cxv[i, yrpick]]]$dat * pixres^2))
  } else if(perm == 1) {
    LL <- lapply(seq_along(dname), function(i) log10(round(obj_data[[i]][[cxv[i, yrpick]]][, "s.area"])))
  }
  
  
  boxplot(LL, at = as.numeric(MAAT), names = as.numeric(MAAT),
          ylab = bquote("log"[10]~"(Area)"), xaxt = "n",
          xlab = bquote("Mean Annual Air Temperature ["~degree*"C]"), 
          main = "Lake Area Distribution (Drier)",
          cex.lab = 1.85, cex.axis = 2, cex.main = 1.6, pars = list(boxwex = 0.3),
          ylim = arealimall, xlim = rev(MAATr))
  grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))
  
  axis(1, at=MAATlab, lab = MAATlab, cex.axis = 1.85, cex.lab = 2)
  
  points(MAAT, znew_l[[perm]][, "mean"],  
         cex = 2, pch = 16)
  
  dev.off()
}


##### QQ plots of lognormal fit


filename = paste0("deltas/plots/Lake_area_lnormQQ.png")
png(file = filename, width = 500, height = 500, pointsize = 14)

xlim <- ylim <- 10^arealimall
xlims_1 <-  ylims_1 <- pretty(log10(xlim))
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 1))

plot(1,1, xlab = bquote("Theoretical Quantiles log"[10]*"(Area [m"^2*"])"),
     ylab = bquote("Empirical Quantiles log"[10]*"(Area [m"^2*"])"),
     cex.lab = 1.6, cex.axis = 1.6,cex.main = 1.6, xlim = xlim, ylim = ylim, 
     log = "xy",
     main = "Lake Lognormal QQ Plot", xaxt = "n", yaxt = "n")
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))


axis(1, at=10^xlims_1, lab=sapply(xlims_1, function(x) {
  as.expression(bquote(.(as.numeric(x))))
}), cex.axis = 1.8, cex.lab = 1.6)
axis(2, at=10^xlims_1, lab=sapply(xlims_1, function(x) {
  as.expression(bquote(.(as.numeric(x))))
}), cex.axis = 1.8, cex.lab = 1.6, las = 2)


abline(0, 1, lwd = 1.5)

lapply(seq_along(dname), function(i) {
  fake_cdf <- emp_cdf(mperm_l[[i]][[cxv[i, yrpick]]]$dat*pixres^2)
  
  Qtheo <- EnvStats::qlnormTrunc(fake_cdf[, 2], 
                                 meanlog = mperm_l[[i]][[cxv[i, yrpick]]]$getPars()[1] + log(pixres^2), 
                                 sdlog = mperm_l[[i]][[cxv[i, yrpick]]]$getPars()[2], 
                                 min = 6*pixres^2, max = Inf)
  QQmat <- cbind(Qtheo, fake_cdf[, 1])
  
  QQmat <- QQmat[!is.infinite(QQmat[, 1]), ]
  
  alpha <- 0.05
  if(perm_KS_pval_bootstrap[cxv[i, yrpick], i] >= alpha){
    lines(QQmat, col = color_scale_vals[color_ranks[dname[i]]], pch = 1, cex = 1, lwd= 2)
  } 
  if(perm_KS_pval_bootstrap[cxv[i, yrpick], i] < alpha) {
    lines(QQmat, col = rgb(71/255, 71/255, 71/255, 0.5), pch = 1, cex = 1, lwd= 2)
  }

})

dev.off()


######## qq plots



filename = paste0("deltas/plots/Waterbody_area_lnormQQ.png")
# png(file = filename, width = 500, height = 1000, pointsize = 14)
png(file = filename, width = 500, height = 500, pointsize = 14)

xlim <- ylim <- 10^arealimall
xlims_1 <-  ylims_1 <- pretty(log10(xlim))
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(1, 1))

plot(1,1, xlab = bquote("Theoretical Quantiles log"[10]*"(Area [m"^2*"])"),
     ylab = bquote("Empirical Quantiles log"[10]*"(Area [m"^2*"])"),
     cex.lab = 1.6, cex.axis = 1.6,cex.main = 1.6, xlim = xlim, ylim = ylim, log = "xy",
     main = "Waterbody Lognormal QQ Plot", xaxt = "n", yaxt = "n")
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))


axis(1, at=10^xlims_1, lab=sapply(xlims_1, function(x) {
  as.expression(bquote(.(as.numeric(x))))
}), cex.axis = 1.8, cex.lab = 1.6)
axis(2, at=10^xlims_1, lab=sapply(xlims_1, function(x) {
  as.expression(bquote(.(as.numeric(x))))
}), cex.axis = 1.8, cex.lab = 1.6, las = 2)

abline(0, 1, lwd = 1.5)

lapply(seq_along(dname), function(i) {
  fake_cdf <- emp_cdf(mall_l[[i]][[cxv[i, yrpick]]]$dat*pixres^2)
  
  Qtheo <- EnvStats::qlnormTrunc(fake_cdf[, 2], 
                                 meanlog = mall_l[[i]][[cxv[i, yrpick]]]$getPars()[1] + log(pixres^2), 
                                 sdlog = mall_l[[i]][[cxv[i, yrpick]]]$getPars()[2], 
                                 min = 6*pixres^2, max = Inf)
  QQmat <- cbind(Qtheo, fake_cdf[, 1])
  
  QQmat <- QQmat[!is.infinite(QQmat[, 1]), ]
  alpha <- 0.05
  if(all_KS_pval_bootstrap[cxv[i, yrpick], i] >= alpha){
    lines(QQmat, col = color_scale_vals[color_ranks[dname[i]]], pch = 1, cex = 1, lwd= 2)
  } 
  if(all_KS_pval_bootstrap[cxv[i, yrpick], i] < alpha) {
    lines(QQmat, col = rgb(71/255, 71/255, 71/255, 0.5), pch = 1, cex = 1, lwd= 2)
  }
  
})

dev.off()

### lake perimeter occurrence

lake_occur <- sapply(obj_data, function(x) sapply(x, function(y) {
    mean(subset(y, (perm == 1) & (s.area >= 1e5) & (s.area <= 1e6))[, "EP1_occur"], na.rm = F)
  }))


ylim_l <- list(c(0, 0.35))
tag = c("EP1")




filename = paste0("deltas/plots/EP1_occurrence.png")
png(file = filename, width = 500, height = 500, pointsize = 14)
par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow=c(1,1))
dapr <- sapply(seq_along(dname), function(i) lake_occur[cxv[i, yrpick], i])

plot(MAAT, dapr,
     pch = dice+14, col = color_scale_vals[color_ranks], 
     cex = 3,
     panel.first = grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5)),
     ylab = "Mean lake border occurrence",
     xlab = bquote("Mean Annual Air Temperature ["*degree*"C]"),
     xlim = rev(MAATr), 
     ylim = ylim_l,
     cex.lab = 1.6, cex.main = 1.6, cex.axis = 1.6)



df <- data.frame(TEMP = MAAT, 
                 EP = dapr)
df <- df[order(df$TEMP), ]
mod <- lm(EP ~ TEMP, data = df)

lines(df$TEMP, ((mod$fitted.values)),  lwd = 3, lty = 2, col = 'black')

corspearman <- cor(df, method = 'spearman')[2]

slope_est <- sapply(1:1e4, function(i) {
  
  ## Randomly shuffle the Y data against X,
  rd <- data.frame(EP = sample(df$EP), TEMP = sample(df$TEMP))
  c(lm(EP ~ TEMP, data = rd)$coefficients, cor(rd, method = 'spearman')[2])
})

# mean(slope_est[2, ] < 0)
# get pval of spearman 
cp = mean((slope_est[3, ] > abs(corspearman)) | (slope_est[3, ] < -abs(corspearman)))
print(paste0("bootstrap spearman p-val is: ", cp))

# pval of regression
qp = mean((slope_est[2, ] > abs(mod$coefficients[2])) | (slope_est[2, ] < -abs(mod$coefficients[2])))
print(paste0("bootstrap p-val is: ", qp))


coef1 <- round(mod$coefficients[1]*1e3)/1e3
coef2 <- round(mod$coefficients[2]*1e2)/1e2
r2 <- round(summary(mod)$r.squared*100)/100
pval <- round(summary(mod)$coefficients[2, 4]*1e4)/1e4

if(pval>0.05) {
  legtrac <- legend("topleft", legend = bquote("R"^2 == .(r2)*"," ~ p == .(pval) ),
                    box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
                    cex = 1.6,trace = T, plot = T)
} else {
  legtrac <- legend("topleft",
                    legend = bquote(atop(bar(Area) %prop% .(coef2)*"MAAT",
                                         "R"^2 == .(r2)*"," ~ p == .(pval) )),
                    box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
                    cex = 1.6,trace = F, plot = T)
}


dev.off()