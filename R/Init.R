## --- Set initial parameters for MLE code


#' Initialise parameter vector for maximum likelihood algorithm.
#' 
#' Find initial estimates for model parameters using smoothing splines.
#'
#' @param ModelInfo Model information object from set_model_info().
#' @param Dat Data frame containg x and y data.
#' 
#' @return Names vector containing initial parameters estimates.
#' 
#' @keywords initial parameter estimates
#'
#' @examples
#'
#' ## Simulate data and get initial model parameter estimates
#' Models <- set_models (mean_fun="gaussian", err_dist="zip")
#' ModelInfo <- set_model_info (Models[1,])
#' Pars   <- create_default_par_list (models=Models)
#' Data   <- create_simulated_datasets (Pars, seed=12345)
#' theta  <- init_mle (ModelInfo, Dat=data.frame(x=Data$x, y=Data$gaussian_zip))
#' print (rbind (Pars[[1]]$theta, theta))
#' 
#' @export
#'
init_mle <- function (ModelInfo, Dat) {
  ## --- Find parameters to initialise mle
  ## --- Lower value of Spar = less smoothing for smoothing spline
  
  ## --- Set negative log-likelihood
  NLL <- set_nll (ModelInfo=ModelInfo)
  
  ## --- Summarise data set
  DF <- init_data_summary (ModelInfo, Dat)
  
  ## --- Loop until negative log-likehood is finite
  Stop <- FALSE
  spar <- NULL
  while (!Stop) {
    
    ## --- Estimate mean function via spline
    MF <- init_spline_mean_function (ModelInfo, DF, spar)
    spar <- MF$spar
    
    ## --- Estimate error distribution parameters
    thetaE <- init_error_par (ModelInfo, DF, MF)
    
    ## --- Estimate mean function parameters
    thetaM <- init_mean_par (ModelInfo, DF, MF)
    
    ## --- Combine parameters
    theta <- c(thetaM, thetaE)
    
    ## --- Transform parameter to be bounded
    u.theta <- make_unbounded (ModelInfo, theta)
    
    ## --- Find negative log-likelihood
    nll <- NLL (u.theta, ModelInfo, Dat)
    
    ## --- Stop if negative log-likelihood is finite
    if (is.finite(nll)) {
      ## Model parameters legal
      Stop <- TRUE
    } else {
      ## Reduce smoothing span
      spar <- 0.8 * spar
    }
  }
  
  ## --- Return parameters
  return (theta)
}


init_data_summary <- function (ModelInfo, Dat) {
  ## ---  Summarise the data for the initialisation algorithm
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Classify data
  err_class <- classify_dataset (Dat)

  ## Grab data
  y <- Dat$y
  x <- Dat$x

  ## Number of data points / positive data points
  nn <- length(y)
  n0 <- sum(y==0)
  n.pos <- sum(y>0)
  
  ## Combine y,x into data frame
  Dat <- data.frame(y=y, x=x)  ## All data, zeroes included
  
  ## Store summary results
  DF <- list()
  DF$Dat <- Dat               ## Data set to use
  DF$err_class <- err_class   ## Type of data
  DF$nn <- nn                 ## Sample size
  DF$n0 <- n0                 ## Number of zeros 
  DF$n.pos <- n.pos           ## Number of positive counts 
  DF$n.pos <- sum(y > 0)      ## Number of positive data points
  
  ## Return results
  return (DF)
}


init_spline_mean_function <- function (ModelInfo, DF, spar=NULL) {
  ## --- Estimate the mean function via a spline
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x

  ## --- Spline fit depends on number of positive y values
  if (DF$n.pos < 1) {

    ## --- Stop if not data points positve
    stop ("No positive data points (response)")

  } else if (DF$n.pos == 1) {
    
    ## --- Only one positive data point
    ## Fit a gaussian curve centered at m, height H
    ## The spread, s, is set to the closest point to
    ## m (provided it is > 0.01).
    m  <- mean(x[which(y>0)])
    H  <- max(y[which(y>0)])
    d  <- abs (x-m)
    s  <- max (min (d[d>0.01]), 0.01)
    mu <- H*exp(-0.5*((x-m)/s)^2)
    spar <- NULL
    
  } else {
    
    ## --- More than one positive data point
    ## Fit a smoothing spline to estimate mu
    if (is.null(spar)) {
      ## Default smoothing span
      S <- stats::smooth.spline(x, y)
      spar <- S$spar
    } else {
      ## Use supplied smoothing span
      S <- stats::smooth.spline(x, y, spar=spar)
    }
    ## Create mean function using spline (force to be non-negative)
    P <- stats::predict(S, data.frame(x=DF$Dat$x))
    mu <- P$y[,1]
    mu[mu < 0] <- 0
    
    ## --- Maximum height
    H <- max(mu)
    ## --- Mode location
    m <- mean(x[mu==max(mu)])
    ## --- Conservative estimate of standard deviation
    s <- diff(range(x[y>0]))/4
  }

  ## Store parameters
  MF <- list()
  MF$mu <- mu
  MF$H <- H
  MF$m <- m
  MF$s <- s
  MF$spar <- spar
  
  ## Return mean function spline object
  return (MF)
}


## --- Error distributions

init_error_par <- function (ModelInfo, DF, MF) {
  ## --- Initialise error distribution parameters

  ## --- Grab mu
  mu <- MF$mu
  
  ## phi : Dispersion paramter
  phi <- init_estimate_phi (ModelInfo, DF, mu)

  ## sigma : sd/var parameter
  sigma <- init_estimate_sigma (ModelInfo, DF, mu)

  ## pi : Zero spike paramter
  pi <- init_estimate_pi (ModelInfo, DF, mu, sigma, phi)

  ## g0, g1 : zero spike / mean linked parameters
  g0 <- init_estimate_linked_intercept (ModelInfo, DF, pi)
  g1 <- init_estimate_linked_slope     (ModelInfo, DF, g0)

  ## If g0 is calculated then pi is not needed
  if (!is.null(g0)) { pi <- NULL }
  
  ## rho : Tweedie exponent
  rho <- init_estimate_rho (ModelInfo, DF)

  ## Error parameter estimates
  thetaE <- c(pi=pi, g0=g0, g1=g1, phi=phi, rho=rho, sigma=sigma)
  
  ## Return parameters
  return (thetaE)
}

init_estimate_phi <- function (ModelInfo, DF, mu) {
  ## --- Estimate the dispersion parameter
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Models with dispersion parameter
  if ( (err_dist=="negbin")  | (err_dist=="zinb"    ) |
       (err_dist=="zinbl" )  | (err_dist=="zinbl.mu") |
       (err_dist=="tweedie") |
       (err_dist=="ziig")    | (err_dist=="ziigl")    | (err_dist=="ziigl.mu") |
       (err_dist=="zig")     | (err_dist=="zigl")     | (err_dist=="zigl.mu")  ) {
    
    ## --- Estimate variance around mean
    SDat <- stats::na.omit (data.frame(x=DF$Dat$x, R=(DF$Dat$y-mu)^2))
    VFit <- with (SDat, stats::smooth.spline (y=R, x=x))
    P    <- stats::predict(VFit, data.frame(x=DF$Dat$x))
    Vy   <- P$y[,1]
    Vy[Vy < 0] <- 0
    
    ## --- Estimate phi
    if (DF$n.pos==1) {
      ## Set phi to 1 if only one positive y
      phi <- 1
    } else {
      phi <- stats::median(abs(Vy - mu)/mu^2, na.rm=TRUE)
      if (is.na(phi)) { phi <- 0 }
    }

  } else {
    ## No dispersion parameter
    phi <- NULL
  }
    
  ## Return phi
  return (phi)
}

init_estimate_sigma <- function (ModelInfo, DF, mu) {
  ## --- Estimate sigma parameter for Gaussian/Beta error distribution

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## --- Default value
  sigma <- NULL
  
  ## --- Gaussian error : sigma is standard deviation
  if (err_dist == "gaussian") {
    sigma <- stats::sd (DF$Dat$y - mu)
  }

  ## --- Beta distribution : sigma is variance 
  if ( (err_dist == "tab") | (err_dist == "zitab") ) {
    ## Use method of moments
    Var <- stats::var (DF$Dat$y - mu)
    v <- Var / (mu*(1-mu) - Var)
    v <- v[v>0]
    if (length (v) > 0) {
      ## MOM
      sigma <- mean(v)
    } else {
      ## Use average variance if MOM is invalid
      sigma <- mean (Var)
    }
  }
  
  ## Return sigma
  return (sigma)
}

init_estimate_pi <- function (ModelInfo, DF, mu, sigma, phi) {
  ## --- Estimate the zero spike proportion parameter
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Grab sample sizes
  nn <- DF$nn
  n0 <- DF$n0
  
  ## --- Estimate pi
  
  ## Default estimated proportion of zeros
  pi <- NULL
  
  ## --- Subtract expected number of zero from gaussian from number of zero

  ## --- ZITAB
  if (err_dist == "zitab") {
    ## Estimate delta
    delta <- ModelInfo$delta
    
    ## Grab beta parameters
    Alpha <- mu / sigma
    Beta  <- (1 - mu) / sigma
    
    ## Calculate pi using MOM
    p0  <- mean (stats::pbeta (delta, shape1=Alpha, shape2=Beta))
    np0 <- nn * p0
    pi  <- (n0 - np0) / (nn - np0)
  }
  
  ## --- ZINB
  if ( (err_dist == "zinb") | (err_dist == "zinbl") | (err_dist == "zinbl.mu") ) {
    ## Probability of zero : prob^(1/phi)
    prob <- (1/(1 + mu*phi))
    
    ## Estimated zeros not from spike
    ns <- sum(prob^(1/phi))
    
    ## Adjust for expected number of zeros
    pi <- n0/nn
    if (ns <= n0) { pi <- (n0 - ns)/nn }
  }
  
  ## --- ZIP
  if ( (err_dist == "zip") | (err_dist == "zipl") | (err_dist == "zipl.mu") ) {
    ## Estimated zeros not from spike
    ns <- sum(exp(-mu))
    
    ## Adjust for expected number of zeros
    pi <- n0/nn
    if (ns <= n0) { pi <- (n0 - ns)/nn }
  }

  ## --- ZIIG (zero inflated inverse-gaussian)
   if ( (err_dist == "ziig") | (err_dist == "ziigl") | (err_dist == "ziigl.mu") )  {
    ## No need to adjust for zeroes not from spike
    pi <- n0/nn
  }
  
  ## --- ZIG
  if ( (err_dist == "zig") | (err_dist == "zigl") | (err_dist == "zigl.mu") )  {
    ## No need to adjust for zeroes not from spike
    pi <- n0/nn
  }
  
  ## Make sure estimate is not negative/too small
  if (!is.null(pi)) {
    if (pi < 0.01) { pi <- 0.01 }
  }
  
  ## Return pi
  return (pi)
}

init_estimate_linked_intercept <- function (ModelInfo, DF, pi) {
  ## --- Estimate the zero spike proportion intercept parameter, gamma 0

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Convert pi to anti-logit scale if linked model
  if ( (err_dist=="zipl")  | (err_dist=="zipl.mu")  |
       (err_dist=="zinbl") | (err_dist=="zinbl.mu") |
       (err_dist=="ziigl") | (err_dist=="ziigl.mu") |
       (err_dist=="zigl")  | (err_dist=="zigl.mu")  ) {
    ## Zero-linked model parameters
    g0 <- log (pi/(1 - pi))
  } else {
    g0 <- NULL
  }
  
  ## Return gamma_0
  return (g0)
}

init_estimate_linked_slope <- function (ModelInfo, DF, g0) {
  ## --- Estimate the zero spike proportion slope parameter, gamma 1

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Initialise g1 to 0 unless g0 is null
  if (is.null(g0)) {
    g1 <- NULL
  } else {
    g1 <- 0
  }
  
  ## Return gamma_1
  return (g1)
}

init_estimate_rho <- function (ModelInfo, DF) {
  ## --- Set rho (tweedie exponent)

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Set rho to 1.1 if tweedie distribution
  if (err_dist == "tweedie" ) {
    rho <- 1.1
  } else {
    rho <- NULL
  }

  ## Return rho
  return (rho)
}


## --- Mean functions

init_mean_par <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters

  ## Grab mean function
  mean_fun <- ModelInfo$mean_fun

  ## --- Constant
  if (mean_fun == "constant") {
    thetaM <- init_mean_constant (ModelInfo, DF)
  }

  ## --- Uniform
  if (mean_fun == "uniform") {
    thetaM <- init_mean_uniform (ModelInfo, DF)
  }

  ## --- Gaussian
  if (mean_fun == "gaussian") {
    thetaM <- init_mean_gaussian (ModelInfo, DF, MF)
  }
  
  ## --- Mixture gaussian
  if ( (mean_fun == "mixgaussian.equal") | (mean_fun == "mixgaussian") ) {
    thetaM <- init_mean_mixgaussian (ModelInfo, DF, MF)
  }

  ## --- Beta
  if (mean_fun == "beta") {
    thetaM <- init_mean_beta (ModelInfo, DF, MF)
  }

  ## --- Sech
  if ( (mean_fun == "sech") | (mean_fun == "sech.p1") | (mean_fun == "sech.r0p1") ) {
    thetaM <- init_mean_sech (ModelInfo, DF, MF)
  }

  ## --- Modskurt
  if (mean_fun == "modskurt") {
    thetaM <- init_mean_modskurt (ModelInfo, DF, MF)
  }
  
  ## --- HOF-II
  if (mean_fun == "hofII") {
    thetaM <- init_mean_hofII (ModelInfo, DF, MF)
  }

  ## --- HOF IV, IVb, V, Vb
  if ( (mean_fun == "hofIV") | (mean_fun == "hofIVb") |
       (mean_fun == "hofV" ) | (mean_fun == "hofVb" ) ) {
    thetaM <- init_mean_hof (ModelInfo, DF, MF)
  }

  ## Return mean function parameters
  return (thetaM)
}

init_mean_constant <- function (ModelInfo, DF) {
  ## --- Estimate mean function parameters : Constant

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Estimate H
  H <- mean (DF$Dat$y)
  H <- init_adjust_H (err_dist, H)
  
  ## Store mean parameters
  thetaM <- c(H=H)

  ## Return mean parameters
  return (thetaM)
}  

init_mean_uniform <- function (ModelInfo, DF) {
  ## --- Estimate mean function parameters : Uniform

  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Fit poisson/bernoulli-uniform
  mle <- mle_uniform_bernoulli (ModelInfo, DF$Dat)
  
  ## Estimate H, c, and d
  H <- as.list(mle)$H
  H <- init_adjust_H (err_dist, H)
  c <- as.list(mle)$c
  d <- as.list(mle)$d
  
  ## Store mean parameters
  thetaM <- c(H=H, c=c, d=d)
  
  ## Return mean parameters
  return (thetaM)
}

init_mean_gaussian <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Gaussian
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Grab parameters from spline fit
  H <- MF$H
  H <- init_adjust_H (err_dist, H)
  m <- MF$m
  s <- MF$s
  
  ## Store mean parameters
  thetaM <- c(H=H, m=m, s=s)
  
  ## Return mean parameters
  return (thetaM)
}

init_mean_mixgaussian <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Mixture gaussian
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  mean_fun <- ModelInfo$mean_fun

  ## Height
  H <- max(MF$H)
  ## Mixture proportion
  a <- 0.5

  ## Place mode location equidistant along range of x where y>0
  m1 <- as.numeric(stats::quantile (min(x[y>0]):max(x[y>0]), c(0.35)))
  m2 <- as.numeric(stats::quantile (min(x[y>0]):max(x[y>0]), c(0.65)))
  
  ## Set standard deviations
  if (mean_fun == "mixgaussian.equal") {
    ## Component standard deviations equal
    s  <- MF$s; s1 <- NULL; s2 <- NULL
  } else {
    ## Component standard deviation different
    s  <- NULL; s1 <- MF$s; s2 <- MF$s
  }
  
  ## Store mean parameters
  thetaM <- c(H=H, a=a, m1=m1, m2=m2, s=s, s1=s1, s2=s2)

  ## Return mean parameters
  return (thetaM)
}
  
init_mean_beta <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Beta
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Fit poisson/bernoulli-uniform
  mle <- mle_uniform_bernoulli (ModelInfo, DF$Dat)
  
  ## Estimate H, c, and d
  H <- as.list(mle)$H
  H <- init_adjust_H (err_dist, H)
  c <- as.list(mle)$c
  d <- as.list(mle)$d
    
  ## Extend initial estimates c and d to avoid mean function = 0 when counts positive
  c <- c - diff(range(DF$Dat$x))/20
  d <- d + diff(range(DF$Dat$x))/20
  
  ## Remove observation outside of [c,d]
  Y <- y[(x>=c) & (x<=d)]
  X <- x[(x>=c) & (x<=d)]
  ## Standardise X to 0 - 1
  X <- (X - c) / (d - c)
  
  ## Use method of moments to find a and b from beta
  ## Check that there is more than one positive data point
  if (sum(Y>0) > 1) {
    
    ## Weighted mean
    wm <- sum(X*Y)/sum(Y)
    
    ## Weighted variance
    wv <- sum (Y*(X-mean(X))^2)/(sum(Y))
    
    ## Method of moments for a and b
    a <- wm * ( wm*(1-wm)/wv - 1 )
    b <- (1-wm) * ( wm*(1-wm)/wv - 1 )
    
  } else {
    ## Lack of data - make it uniform
    a <- 1
    b <- 1
  }
  
  ## If u-shaped or j-shaped set to uniform
  if ((a < 1) | (b < 1)) {
    a <- 1
    b <- 1
  }
  
  ## Make sure a & b are not too large
  if (a > 5) { a <- 5 }
  if (b > 5) { b <- 5 }
  
  ## Convert to mean and shape parameter
  u <- a / (a + b)
  v <- 1 / (a + b)
  
  ## Store mean parameters
  thetaM <- c(H=H, c=c, d=d, u=u, v=v)

  ## Return mean parameters
  return (thetaM)
}

init_mean_sech <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Sech
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Grab mean function
  mean_fun <- ModelInfo$mean_fun

  ## --- Weighted variance of x
  
  ## Mean and variance of x (weighted by y)
  wm <- sum(x*y)/sum(y)
  wv <- sum (y*(x-mean(x))^2)/(sum(y))
    
  ## --- Estimate parameters from spline
  H <- MF$H
  H <- init_adjust_H (err_dist, H)
  m <- MF$m
  s <- MF$s

  ## Make sure s is not zero if only one data point
  s <- max(sqrt(wv), min(diff(unique(sort(DF$Dat$x)))))
  if ( is.na(s) | is.nan(s) ) { s <- 1 }
  if ( sum(y>0)==1 ) { s <- sqrt(mean((x-m)^2))/6 }
  r <- 0
  p <- 1
  
  ## Store mean parameters
  if (mean_fun == "sech"     ) { thetaM <- c(H=H, m=m, s=s, r=r, p=p) }
  if (mean_fun == "sech.p1"  ) { thetaM <- c(H=H, m=m, s=s, r=r) }
  if (mean_fun == "sech.r0p1") { thetaM <- c(H=H, m=m, s=s) }
  
  ## Return mean parameters
  return (thetaM)
}

init_mean_modskurt <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Modskurt
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Grab mean function
  mean_fun <- ModelInfo$mean_fun

  ## --- Weighted variance of x
  
  ## Mean and variance of x (weighted by y)
  wm <- sum(x*y)/sum(y)
  wv <- sum (y*(x-mean(x))^2)/(sum(y))
    
  ## --- Estimate parameters from spline
  H <- MF$H
  H <- init_adjust_H (err_dist, H)
  m <- MF$m
  s <- MF$s

  ## Make sure s is not zero if only one data point
  s <- max(sqrt(wv), min(diff(unique(sort(DF$Dat$x)))))
  if ( is.na(s) | is.nan(s) ) { s <- 1 }
  if ( sum(y>0)==1 ) { s <- sqrt(mean((x-m)^2))/6 }
  q <- 0.5
  p <- 1
  b <- 2
  
  ## Store mean parameters
  thetaM <- c(H=H, m=m, s=s, q=q, p=p, b=b)
  
  ## Return mean parameters
  return (thetaM)
}

init_mean_hofII <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Hof-II
  ## mu = H * (1/( 1 + exp(-w0*(x-m)) ))

  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist

  ## Grab mu
  mu <- MF$mu
  
  ## --- Estimate parameters from spline
  H <- max(MF$H)
  H <- init_adjust_H (err_dist, H)

  ## Find x location closest to H/2
  e <- abs(mu - H/2)
  m <- x[which (e==min(e))[1]]

  ## The slope of the hofII curve at x=m equals H*w0/4
  ## We estimate the slope at x=m by considering points
  ## at x= m - (1/w0), and x= m + (1/w0)
  ## Or at where -w0(x-m) = {-1, 1}
  
  ## Y values at  m - (1/w0), m + (1/w0)
  UQ <- H * (1/(1 + exp(1)))
  LQ <- H * (1/(1 + exp(-1)))
  ## Find x locations closest to where mu = UQ
  e  <- abs(mu - UQ)
  xUQ <- DF$Dat$x[which (e==min(e))[1]]
  ## Find x locations closest to where mu = LQ
  e  <- abs(mu - LQ)
  xLQ <- DF$Dat$x[which (e==min(e))[1]]
  ## w0 is estimated to be slope / (H/4)
  w0 <- ((UQ-LQ)/(xUQ-xLQ)) / (H/4)
  
  ## Store mean parameters
  thetaM <- c(H=H, m=m, w0=w0)
  
  ## Return mean parameters
  return (thetaM)
}

init_mean_hof <- function (ModelInfo, DF, MF) {
  ## --- Estimate mean function parameters : Hof-IV / Hof-IVb / Hof-V / Hof-Vb
  
  ## --- Grab data
  y <- DF$Dat$y
  x <- DF$Dat$x
  
  ## Grab error distribution
  err_dist <- ModelInfo$err_dist
  
  ## Grab mean function
  mean_fun <- ModelInfo$mean_fun
  
  ## --- Estimate parameters from spline
  H <- max(MF$H)
  H <- init_adjust_H (err_dist, H)
  m <- MF$m
  w <- 4*exp(-.5)/MF$s
  
  ## Store mean parameters
  if (mean_fun == "hofIV"  ) { thetaM <- c(H=H, m=m, w=w, k=1)        }
  if (mean_fun == "hofIVb" ) { thetaM <- c(H=H, m=m, w=w)             }
  if (mean_fun == "hofV"   ) { thetaM <- c(H=H, m=m, w1=w, w2=w, k=1) }
  if (mean_fun == "hofVb"  ) { thetaM <- c(H=H, m=m, w1=w, w2=w)      }
  
  ## Return mean parameters
  return (thetaM)
}

init_adjust_H <- function (err_dist, H) {
  ## --- Adjust H for models where y is [0,1] so it is not on the boundary
  if ( (err_dist == "bernoulli") | (err_dist == "tab") | (err_dist == "zitab") ) {
    if (H > 0.9) { H <- 0.9 }
  }
  
  ## Return H
  return (H)
}

MOM.zinb <- function (Dat) {
  ## Fit a zero-inflated negative binomial distribution
  ## via method-of-moments

  ## Grab count 
  y <- Dat$y

  ## Calculate mean, variance, and proportion of zeros
  ym <- mean (y)
  s2 <- stats::var (y)
  p0 <- mean (y==0)

  ## MOM equation for phi 
  RootPhi <- function (phi) {
    ## Calcualate constants
    k0 <- ((s2 - ym - ym^2*phi)/(s2 - ym*(1-ym)))
    k1 <- ym^2*(1 + phi)/(s2-ym*(1-ym))
    k2 <- ((ym+(s2+ym^2)*phi)/(ym*(1 + phi)))^(-1/phi)
    ## Calculate error
    err <- p0 - k0 - k1*k2
    return (err)
  }

  ## MOM estimates
  phi <- stats::uniroot(RootPhi, lower = -0.5, upper = 1E50)$root
  mu  <- (s2 - ym*(1-ym))/(ym*(1 + phi))
  pi  <- 1 - ym/mu

  ## Return parameter estimates
  theta <- c(pi=pi, phi=phi, mu=mu)
  return (theta)
}


MOM.nb <- function (Dat) {
  ## Fit a zero-inflated negative binomial distribution
  ## via method-of-moments
  
  ## Grab count 
  y <- Dat$y

  ## Calculate mean, variance, and proportion of zeros
  ym <- mean (y)
  s2 <- stats::var (y)
  
  ## MOM estimates
  phi <- (s2 - ym)/ym^2
  mu  <- ym

  ## Return parameter estimates
  theta <- c(phi=phi, mu=mu)
  return (theta)
}
