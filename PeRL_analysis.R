## Analysis of the PeRL database provided in Muster et al., 2019, Frontiers.
## Author: Lawrence Vulis; lvulis [a] uci.edu
## Goal: Compute statistical moments of observed waterbody sizes in the PeRL and GSW-derived datasets

##### LOAD LIBRARIES AND CUSTOM FUNCTIONS ######
library(sf)
library(poweRlaw)
library(EnvStats)
library(plotly)

summary_stats <- function(kk) {
  return(c(mean = mean(kk),
           median = median(kk),
           var = var(kk),
           N = length(kk),
           skew = skewness(kk),
           kurt = kurtosis(kk),
           l.skew = EnvStats::skewness(kk, method = 'l.moment'),
           l.kurt = EnvStats::kurtosis(kk, method = 'l.moment'),
           max = max(kk),
           A_total = sum(kk)))
}

source("custom_fcns.R")

### To get GSW comparison, go to "plot_together_multipleyr" and run the sections
### up to fitting the lognormal. 

a <- 6*pixres^2
b <- 1e6
# b <- +Inf

lake_summaries <- lapply(obj_data, function(x) {
  # mean(x[[1]][, 's.area']) > mean(x[[2]][, 's.area'])
  
  sapply(x, function(y) {
    y <- y[y[, "perm"] == 1, ]
    mm <- y[, 's.area']
    mm <- mm[(mm>=6*pixres^2) & (mm<=b)]/1e6 #truncate above 6 pix & below 1e6 m2 and convert to km2
    summary_stats(mm)
    
  })
})

#lnorm theoretical mean and variance:

theo_moment_ln_perm <- t(sapply(seq_along(mperm_l), function(i) {
  mux <- mperm_l[[i]][[cxv[i, yrpick]]]$getPars()[1] + log(pixres^2)
  sdx <- mperm_l[[i]][[cxv[i, yrpick]]]$getPars()[2]
  
  muhat <- lnorm_meanx(mux, sdx, xmin = a, xmax = b)
  varhat <- lnorm_varx(mux, sdx, xmin = a, xmax = b)
  skewhat <- lnorm_skewx(mux, sdx, xmin = a, xmax = b)
  c(mean = muhat, var = varhat, skew = skewhat) # * c(p_obs, p_obs^2)
  
}))



cor(MAAT, sapply(seq_along(lake_summaries), function(i) lake_summaries[[i]]["mean", cxv[i, yrpick]]))

#### LOAD PeRL DATA ####

fnlist <- read.csv("deltas/PeRL/muster_2019_regions.csv", header = F,
               stringsAsFactors = F)[, 1]


v<- lapply(fnlist, function(x) {
  fn <- list.files("deltas/PeRL/waterbodies", paste0(x, ".*shp$"))
  if(length(grep("che0012", x)) > 0) {fn = fn[2]}
  dat <- st_read(paste0("deltas/PeRL/waterbodies/", fn))
  

  return(cbind(dat$AREA, dat$PERIMETER))
  })

##### ANALYSIS OF PERL DATA ########

a <- 100 # replicate muster, use 5600 otherwise (6*900)
v_area <- lapply(v, function(x) {
  x <- x[x[, 1] >= a, 1]
  x
})

b <- 1e6

# fnlist$V1[c(1, 4, 5, 6)]


bins <- seq(2, 8, by = 0.1)

hist_real <- lapply(v_area, function(x) {
  hist(log10(x), breaks = bins, plot = T)
})

real_lnormobj <- lapply(v_area, function(x) {
  x<-round(x)
  jaj <- dislnorm$new(x)
  jaj$setXmin(100)
  jaj$setPars(estimate_pars(jaj))
  jaj$internal$gof <- get_distance_statistic(jaj)
  return(jaj)
})

CCDF_real <- lapply(real_lnormobj, function(x) cbind(unique(x$dat), dist_data_cdf(x, xmax = max(x$dat)+1, lower_tail = F)))

## add kurtosis
trimmer <- c(5400, 1e6) # 0 means no minimum. 5400 is to match the one from GSW
trunc_moment_tab_real <- data.frame(t(sapply(v_area, function(x) {
  # trimmer <- quantile(x, c(0, .99))
  # log10()
  x <- x[(x>= trimmer[1]) & (x<= trimmer[2])]/1e6
  c(mean = mean(x), var = var(x), skew = skewness(x, T), kurt = kurtosis(x),
    l.skew = EnvStats::skewness(x, method = 'l.moment'),
    l.kurtosis = EnvStats::kurtosis(x, method = 'l.moment'))
})))

moment_tab_real <- data.frame(t(sapply(v_area, function(x) {
  x<-x/1e6
  c(mean = mean(x), var = var(x), skew = skewness(x, T), kurt = kurtosis(x),
    l.skew = EnvStats::skewness(x, method = 'l.moment'),
    l.kurtosis = EnvStats::kurtosis(x, method = 'l.moment'))
})))


Pupper_real <- sapply(v_area, function(x) {
  sum(x <= b)/length(x)
})


### estimate 

color_scale_vals <- hcl.colors(length(v_area), palette = "Zissou1")
color_scale_ranks <- round(rank(moment_tab_real$mean))

color_scale_vals <- rep('black', nrow(trunc_moment_tab_real))
GSW_col <- '#5200DF'

GSW <- T # adds GSW data to plot

filename = paste0("deltas/PeRL/plots/realtruncmeanvarskew.png")

png(file = filename, width = 800, height = 400, pointsize = 14)

par(mar = c(5.1, 5.1, 3.9, 2.1), mfrow = c(1, 2))


plot(trunc_moment_tab_real[, c(1, 2)], xlab = "Mean",
     ylab = "Variance", main = "Var. vs. Mean", cex = 2,
     cex.lab = 1.6, cex.main = 1.6, cex.axis = 1.6, pch = 15, 
     col = color_scale_vals[color_scale_ranks],
     ylim = c(0, 0.08))
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))



dfmod <- trunc_moment_tab_real
dfmod <- dfmod[order(dfmod$mean), ]
dfmod <- dfmod[-which((dfmod$mean > 0.2) & (dfmod$var > 0.04)), ] # mean > 0.2 & var >0.04 for x in [5400, 1e6]
mod <- lm(var ~ mean, data = data.frame(mean = dfmod$mean, var = dfmod$var))
lines(dfmod$mean, mod$fitted.values, lwd = 2, lty = 2)


coef1 <- round(mod$coefficients[1])
coef2 <- round(mod$coefficients[2]*100)/100
r2 <- round(summary(mod)$r.squared*100)/100
pval <- round(summary(mod)$coefficients[2, 4])


legend("topleft",
       legend = bquote(atop("Var" %prop% .(coef2)*"Mean",
                            "R"^2 == .(r2)*"," ~ p == .(pval) )),
       box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
       cex = 1.6, text.font = 2)

### GSW

if(GSW) {
  
  ### Use empirical moments
  # dfm <- data.frame(do.call(rbind, lapply(seq_along(lake_summaries), function(i) {
  #   lake_summaries[[i]][, cxv[i, yrpick]][c("mean", "var")]
  # })))
  ### Use LN-based estimates
  dfm <- data.frame(sweep(theo_moment_ln_perm[, 1:2], 2, c(1e6, 1e12), "/"))
  
  ### computations
  colnames(dfm) <- c("mean", "var")
  dfm <- dfm[order(dfm$mean), ]
  points(dfm, pch = 17, cex = 2, col = GSW_col)
  mod2 <- lm(var ~ mean, data = dfm)
  lines(dfm$mean, mod2$fitted.values, col = GSW_col)
  
  
  coef1_GSW <- round(mod2$coefficients[1])
  coef2_GSW <- round(mod2$coefficients[2]*100)/100
  r2_GSW <- round(summary(mod2)$r.squared*100)/100
  pval_GSW <- round(summary(mod2)$coefficients[2, 4])
  
  legend("bottomright",
         legend = bquote(atop("Var" %prop% .(coef2_GSW)*"Mean",
                              "R"^2 == .(r2)*"," ~ p == .(pval_GSW) )),
         box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = GSW_col,
         cex = 1.6, text.font = 2)
  
  

}

plot(trunc_moment_tab_real[, c(1, 3)], xlab = "Mean", ylab = "Skewness coef.",  main = "Skewness Coef. vs. Mean",
     cex = 2, cex.lab = 1.6, cex.main = 1.6, cex.axis = 1.6, pch = 15, col = color_scale_vals[color_scale_ranks])
grid(lwd = 2, col = rgb(0.05, 0.05, 0.05, 0.5))

dfmod <- trunc_moment_tab_real
dfmod <- dfmod[order(dfmod$mean), ]
dfmod <- dfmod[-which((dfmod$mean < 0.05) & (dfmod$skew <= 3)), ]
mod <- lm(skew ~ mean, data = data.frame(mean = 1/sqrt(dfmod$mean), skew = dfmod$skew))
lines(dfmod$mean, mod$fitted.values)

coef1 <- round(mod$coefficients[1])
coef2 <- round(mod$coefficients[2]*100)/100
r2 <- round(summary(mod)$r.squared*100)/100
pval <- round(summary(mod)$coefficients[2, 4]*1e4)/1e4

legend("topright",
       legend = bquote(atop("Sk. Coef. " %prop% frac(.(coef2), sqrt("Mean")),
                            "R"^2 == .(r2_GSW)*"," ~ p == .(pval) )),
       box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = 'black',
       cex = 1.6, text.font = 3)


if(GSW) {
  ### Use empirical moments
  # dfm <- data.frame(do.call(rbind, lapply(seq_along(lake_summaries), function(i) {
  #   lake_summaries[[i]][, cxv[i, yrpick]][c("mean", "skew")]
  # })))
  ### 
  dfm <- data.frame(theo_moment_ln_perm[, c(1, 3)])
  dfm[, 1] <- dfm[, 1]/1e6
  ### Computations
  colnames(dfm) <- c("mean", "skew")
  dfm <- dfm[order(dfm$mean), ]
  points(dfm, pch = 17, cex = 2, col = GSW_col)
  mod2 <- lm(skew ~ mean, data = data.frame(mean = 1/sqrt(dfm$mean), skew = dfm$skew))
  lines(dfm$mean, mod2$fitted.values, col = GSW_col)
  
  
  
  coef1_GSW <- round(mod2$coefficients[1])
  coef2_GSW <- round(mod2$coefficients[2]*100)/100
  r2_GSW <- round(summary(mod2)$r.squared*100)/100
  pval_GSW <- round(summary(mod2)$coefficients[2, 4])
  
  
  legtrac <- legend("bottomleft",
         legend = bquote(atop("Sk. Co."%prop% frac(.(coef2_GSW), sqrt("Mean")),
                              "R"^2 == .(r2_GSW)*"," ~ p == .(pval_GSW) )),
         box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = GSW_col,
         cex = 1.2, text.font = 3, trace = F, plot = T, inset = c(-0.1, 0.05))
  # legend(x = legtrac$text$x, y = legtrac$text$y,
  #                   legend = bquote(atop("Skew" %prop% frac(.(coef2_GSW), sqrt("Mean")),
  #                                        "R"^2 == .(r2_GSW) ~ p == .(pval_GSW) )),
  #                   box.lwd = 0, box.col = "white", bg = "white", bty = "n", text.col = GSW_col,
  #                   cex = 1.4, text.font = 3, trace = F, plot = T)
}


dev.off()

