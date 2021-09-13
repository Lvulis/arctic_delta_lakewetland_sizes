### Useful functions used throughotu the analysis
library(compiler)

# Spatial functions based on sf and raster
clip <- function(x,shape) {
  # Arguments:
  #   x: Input raster
  #   shape: polygon (in same projection)
  # Returns:
  #   clipped raster
  
  return(mask(crop(x, shape), shape))
}

clipper <- function(x, shape, cr, ext) {
  # Perform polygon clipping using raster on an existing matrix through transform
  # to and from raster class
  # Arguments:
  #  x: matrix-stored raster to be clipped
  #  shape: sf polygon 
  #  cr: PROJ4 string
  #  ext: extent 
  # Returns: 
  #  Clipped matrix x
  x <- raster(t(x))
  crs(x) <- cr
  extent(x) <- ext
  x <- clip(x, shape)
  x <- t(as.matrix(x))
  return(x)
}

idExBase <- function(id, nr, nc) {
  # Retrieve x and y indices of an N by M matrix
  # Arguments:
  #   A: Rectangular matrix
  # Returns:
  #   NxM by 2 matrix containing indicse of every point
  idx <- id %% nr
  if(length(idx)==1) {
    if(idx==0) {
      idx <- nr
    }
  } else {
    idx[idx==0] <- nr
  }
  idy <- (id-idx)/nr+1
  return(cbind(idx, idy))
}

idEx <- function(G) {
  # Retrieve x and y indices of an N by M matrix
  # Arguments:
  #   A: Rectangular matrix
  # Returns:
  #   NxM by 2 matrix containing indicse of every point
  nr2 <- nrow(G)
  nc2 <- ncol(G)
  id2 <- 1:length(G)
  out <- idExBase(id2, nr2, nc2)
  return(out)
}

idIn <- function(idx, idy, nro, nco) {
  # function to convert x, y indices to their running index value
  # Arguments:
  #  idx: row index
  #  idy: column index
  #  nro: number of rows in the matrix
  #  nco: number of columns in the matrix
  # Returns:
  #  vector of running indices
  ifelse(((0 < idx) & (idx <= nro)) & ((0 < idy) & (idy <= nco)),
         round((idy-1)*nro + idx),
         NA)

}

edgeId <- function(X) {
  # Get index of edge pixels
  # Arguments:
  #  X: a matrix
  # Returns:
  #  vector of edge indices
  n <- nrow(X)
  m <- ncol(X)
  return(unique(c(1:n, n*(0:(m-1))+1, (n*(m-1)+1):(n*m), (1:m)*n)))
}

### Functions for labelled images
measureProps <- function(x, pixres = 1) {
  # Use skimage measure (calling reticulate) to compute feature properties
  # Arguments:
  #  x: Labelled input integer matrix (no NA)
  #  pixres: spatial resolution of image (default 1 )
  # Return:
  #  prop: matrix containing object properties (n x p) n objects and p proprties (nxp)
  outtie <- measure$regionprops(x)
  oid <- sort(unique(as.vector(x)))
  oid <- oid[oid>0]
  prop <- t(sapply(outtie, function(y) c(y$area, y$perimeter)))
  colnames(prop) <- c("s.area", "s.perimeter")
  rownames(prop) <- oid
  prop[, "s.area"] <- prop[, "s.area"]*pixres^2
  prop[, "s.perimeter"] <- prop[, "s.perimeter"]*pixres
  return(prop)
}


#### Useful math to have on hand
erfc <- function(x) {
  # Calculate complementary error function of x
  2 * pnorm(-sqrt(2) * x)
}


kurtosis <- function(x) {
  # Estimate the kurtosis of a distribution
  # Arguments:
  #  x: vector of observations (integer or double)
  # Returns:
  #  kurtosis ofthe data
  n <- length(x)
  m4 <- (x - mean(x))^4
  m2 <- (x - mean(x))^2
  return(((sum(m4) / n) / ((sum(m2) / n)^2)) - 3 )
}

skewness <- function(x, FP = T) {
  # Compute skewness coefficient with optional fisher-pearson correction
  # Arguments:
  #  x: vector of observations (integer or double)
  # Returns:
  #  skewness of the data
  n<-length(x)
  m3 <- (x - mean(x)) ^ 3
  sk <- (sum(m3) / n) / (sd(x) ^ 3)
  if(FP == T) {
    sk <- (n ^ 2 / ((n - 1) * (n - 2))) * sk
  }
  return(sk)
}

### NEW EMP_CDF TO HAVE STEP FUNCTION INCLUDED...HOW TO DO THIS:

emp_cdf <- function(x){
  cd <- ecdf(x)

  return( cbind(sort(unique(x)), cd(sort(unique(x)))))
}
emp_ccdf <- function(x) {
  CDF <- emp_cdf(x)
  return(cbind(CDF[, 1], 1 - CDF[, 2] + min(CDF[, 2])))
}

# This is the recommended way for certain resamplings
resample <- function(x, ...) x[sample.int(length(x), ...)]

log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceil(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}

adjusted_dist_cdf <- function(x, cut = FALSE, length.out = 100) {
  scale = 1
  xmin = x$getXmin()
  x_values = x$dat
  
  x_axs = lseq(xmin, max(x_values), length.out)
  if (is(x, "discrete_distribution")) 
    x_axs = unique(round(x_axs))
  y = dist_cdf(x, x_axs, F)
  if (!cut) {
    if (is(x, "discrete_distribution")) {
      x$setXmin(1)
    } else {x$setXmin(0)}
    d_cdf = dist_data_cdf(x, lower_tail = F)
    if (is(x, "discrete_distribution")) {
      dif = x$internal[["values"]] - xmin
      upper = which(dif > 0)[1]
      lower = max(upper - 1, 1)
      x_dif = x$internal[["values"]][lower] - x$internal[["values"]][upper]
      y_dif = d_cdf[lower] - d_cdf[upper]
      scale = d_cdf[lower] + y_dif * (xmin - x$internal[["values"]][lower])/x_dif
    } else {
      dif = x$internal[["dat"]] - xmin
      upper = which(dif >= 0)[1]
      lower = max(upper - 1, 1)
      x_dif = x$internal[["dat"]][lower] - x$internal[["dat"]][upper]
      y_dif = d_cdf[lower] - d_cdf[upper]
      scale = d_cdf[lower] + y_dif * (xmin - x$internal[["dat"]][lower])/x_dif
    }
    x$setXmin(xmin)
  }
  y = y * scale
  return(data.frame(x = x_axs, y = 1- y))
}


### Analytical lognormal moments for the Muster et al. check


lnorm_trunc_crxn <- function(mu, sigma, r, xmin = 0, xmax = +Inf) {
  # Correction factor for LN truncated moments
  # Can derive from Jawitz et al. 2004
  # Arguments:
  #  mu: meanlog parameter (gamma in paper)
  #  sigma: sdlog paramter (beta in paper)
  #  xmin: left truncation point
  #  xmax: right truncation point
  # Returns:
  #  constant corresponding to a correction factor
  a0 <- (log(xmin) - mu)/sigma
  b0 <- (log(xmax) - mu)/sigma
  (pnorm(r*sigma - a0) - pnorm(r*sigma - b0)) / (pnorm(b0) - pnorm(a0)) 
}

lnorm_moment <- function(mu, sigma, r, xmin = 0, xmax = +Inf) {
  # Analytical raw moments of the truncated lognormal lognormal distribution
  # Arguments:
  #  mu: meanlog parameter (gamma in paper)
  #  sigma: sdlog paramter (beta in paper)
  #  xmin: left truncation point
  #  xmax: right truncation point
  # Returns:
  #  constant corresponding to a correction factor
  mom <- exp(r*mu + r^2*sigma^2/2) * lnorm_trunc_crxn(mu, sigma, r, xmin = xmin, xmax = xmax)
  return(mom)
}

# Mean, Variance, and Skew of truncated lognormal
lnorm_meanx <- function(mu, sigma, xmin = 0, xmax = +Inf, left = T) {
  lnorm_moment(mu = mu, sigma = sigma, r = 1, xmin = xmin, xmax = xmax)
}
lnorm_varx <- function(mu, sigma, xmin = 0, xmax = +Inf) {
  lnorm_moment(mu, sigma, r = 2, xmin = xmin, xmax = xmax) - lnorm_moment(mu, sigma, r = 1, xmin = xmin, xmax = xmax)^2
}
lnorm_skewx <- function(mu, sigma, xmin = 0, xmax = +Inf) {
  meanx <- lnorm_moment(mu, sigma, r = 1, xmin = xmin, xmax = xmax)
  mom2 <- lnorm_moment(mu, sigma, r = 2, xmin = xmin, xmax = xmax)
  mom3 <- lnorm_moment(mu, sigma, r = 3, xmin = xmin, xmax = xmax)
  (mom3 - (3 * meanx * (mom2 - meanx^2)) - meanx^3 ) / (mom2 - meanx^2)^(3/2)

}

lines_dist_local <- function(x, xscale = 1, yscale = 1, cut = FALSE, draw = TRUE, length.out = 100, ...) {
  scale = 1
  xmin = x$getXmin()
  x_values = x$dat
  x_axs = lseq(xmin, max(x_values), length.out)
  if (is(x, "discrete_distribution")) 
    x_axs = unique(round(x_axs))
  y = dist_cdf(x, x_axs, FALSE)
  if (!cut) {
    if (is(x, "discrete_distribution")) 
      x$setXmin(1)
    else x$setXmin(0)
    d_cdf = dist_data_cdf(x, lower_tail = FALSE)
    if (is(x, "discrete_distribution")) {
      dif = x$internal[["values"]] - xmin
      upper = which(dif > 0)[1]
      lower = max(upper - 1, 1)
      x_dif = x$internal[["values"]][lower] - x$internal[["values"]][upper]
      y_dif = d_cdf[lower] - d_cdf[upper]
      scale = d_cdf[lower] + y_dif * (xmin - x$internal[["values"]][lower])/x_dif
    }
    else {
      dif = x$internal[["dat"]] - xmin
      upper = which(dif >= 0)[1]
      lower = max(upper - 1, 1)
      x_dif = x$internal[["dat"]][lower] - x$internal[["dat"]][upper]
      y_dif = d_cdf[lower] - d_cdf[upper]
      scale = d_cdf[lower] + y_dif * (xmin - x$internal[["dat"]][lower])/x_dif
    }
    x$setXmin(xmin)
  }
  if (is.nan(scale)) 
    scale = 1
  y = y * scale
  if (draw) 
    lines(x_axs*xscale, y*yscale, ...)
  invisible(data.frame(x = x_axs*xscale, y = y*scale))
}

keepLargest <- function(img) {
  # Keep only largest element of an image. Relies on python
  # Arguments:
  #  img: integer-valued binary matrix
  # Returns:
  #  Same sized matrix only containing largest cluster (label id still there)
  cncmp <- measure$label(img, connectivity = 2)
  mode(cncmp) <- 'integer'
  props <- measureProps(cncmp)
  tokp <- names(which.max(props[, 's.area']))
  clean_mask <- rmObjects(cncmp, setdiff(rownames(props), tokp), reenumerate = F)
  return(clean_mask)
}

splitObjects = function(x) {
  # From EBImage: get R indices of each object in x matrix 
  # Arguments: 
  #  x: Integer valued labelled image
  # Returns:
  #  list w/ each entry corresponding to indices of a unique label in x
  z = which(as.integer(x)>0L)
  split(z, x[z])
}


bin_thresh <- function(x, thresh) {
  # Binary threshold a numeric/integer vector or matrix
  # Arguments:
  #  x: Matrix
  #  thresh: threshold value
  # Returns:
  #  x: Binary matrix 
  x[x < thresh] <- 0L
  x[x >= thresh] <- 1L
  mode(x) <- 'integer'
  x
}

colorProp <- function(x, newLab) {
  # Color labelled objects on X with values given by newLab.
  # Can use to color labels (e.g. objects) by their properties
   # Arguments:
  #  x: Matrix
  #  newLab: named vector where names = labels and values = update
  # Returns:
  #  x: Binary matrix
  xs = splitObjects(x)
  if (length(xs) == 0) 
    return(NULL)
  for(obj in names(xs)) {
    print(obj)
    x[xs[[obj]]] <- newLab[obj]
  }
  return(x)
}



diladd <- function(x, dil) {
  # compute added dilation features on labelled (essentially grayscale) img
  # Has a failsafe if two objects crash into one another
  # Arguments:
  #  x: labelled matrix
  #  dil: dilation kernel
  # Returns:
  #  dilated matrix
  domaintemp <- dilate(x, dil) - x
  oldarea <- which(x>0)
  oldarea2 <- oldarea[domaintemp[oldarea] > 0]
  domaintemp[oldarea2] <- x[oldarea2] 
  return(domaintemp)
}

empiricalQQ <- function(x, y) {
  # Construct empirical QQ plot of two univariate distributions
  # Arguments:
  #  x: vector of observations
  #  y: vector of observations
  CDF_1 <- emp_cdf(x)
  out2 <- quantile(y, prob = CDF_1[, 2])
  cbind(CDF_1[, 1], out2)
}

sumStat <- function(x) {
  # Compute mean, median, sd, and skewness of a vector of data x
  c(mean = mean(x), median = median(x), sd = sd(x), skew = skewness(x))
}



### Implementations of poweRlaw's functions for changing
### parallelization procedure
sample_p_helper = function(i, m, x_lower) {
  ## Total sample size
  N = get_n(m)
  ntail_prop = get_ntail(m, prop = TRUE)
  
  ## Proportion to sample
  n1 = sum(runif(N) > ntail_prop) # less than xmin
  
  # q should be of length N
  c(sample(x_lower, n1, replace = TRUE), #less than xmin
    dist_rand(m, N - n1))
}

bootstrap_p_helper = function(i, m, x_lower, xmins, pars, xmax, distance) {
  
  q = sample_p_helper(i, m, x_lower)
  m_cpy = m$getRefClass()$new(q)
  
  est = estimate_xmin(m_cpy, xmins = xmins, pars = pars, xmax = max(m_copy)+1, distance = distance)
  ## Remove the character now, since we will change to data frame.
  est["distance"] = NULL
  unlist(est)
}