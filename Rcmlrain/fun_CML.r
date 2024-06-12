# Title: R functoins for CML processing
# Version: 0.0
# Date: 2017-09-03
# Author: Martin Fencl <martin.fencl@cvut.cz>
# Maintainer: Martin Fencl <martin.fencl@cvut.cz>
# Description: Provides set of tools for processing cml data. Designed for 
# ICUD CML workshop in Prague COngress Centre (10th September 2017)
# License: MIT
# Last upadate: 2022/02 (not exactly when and what)




## =======================
## Preprocessing raw data:
## =======================

## Set of functions for processing of raw data in a form of non-regular time series

identify_peaks <- function(tl, tsh, report = T) {
  ## Identifies position of suspicious sudden peaks in the vector of total losses
  ##
  ## Arguments:
  ## tl - vector or time series with total loss or attenuation values
  ## tsh - threshold value in dB to classify data point as sudden peak
  ## report - logical scalar. If True, number of identified peaks is printed
  ## Returns: id.peak - indexes of data points being identified as sudden peaks
  
  # first derivative
  tl <- as.numeric(tl)
  dtl <- tl[-1]-tl[-length(tl)]
  
  # detect local extremes
  local_min <- c(F, abs(dtl[-1] > 0) + (dtl[-length(dtl)] < 0) == 2, F)
  local_max <- c(F, abs(dtl[-1] < 0) + (dtl[-length(dtl)] > 0) == 2, F)
  local_ext <- (local_min + local_max) > 0
  
  # quantify magnitude of peaks and troughs
  dtl2 <- apply(data.frame(abs(c(dtl,0)), abs(c(0,dtl))), 1, min) * local_ext
  
  # identify position of peaks/troughs exceeding threshold value
  id.peak <- which(abs(dtl2) > tsh) 
  
  if(report==T){
    print(paste("Number of sudden peaks:", length(id.peak)))    
  }
  
  
  return(id.peak)
}

#--------------------------------------

zoo_aggreg_by <- function(ts_zoo, step, fun, align = 'center',
                          insert.missing = T, ...){

    ## Aggregates a regular or non-regular time series to the regular time series of defined time step
    ##
    ## parameters:
    ## ts_zoo - zoo time series
    ## step - time step to which agregate the series (in minutes)
    ## step - time step to which agregate the series (in minutes)
    ## 
    ## returns:
    ## ag_zoo - aggregated regular time series of defined time step.
    
    require(zoo)
    tim <- index(ts_zoo)
    
    t_num <- as.numeric(tim)
    
    if(align == 'left'){
        t_num <- floor(t_num/(step*60))*step*60    
    }else if(align == 'center'){
        t_num <- round(t_num/(step*60), 0)*step*60
    }else if(align == 'right'){
        t_num <- ceiling(t_num/(step*60))*step*60
    }else{stop('Error in match.arg(align) : \n 
                arg should be one of center, left, right)')}
    
    
    tim_agr <-  as.POSIXct(t_num, origin="1970-01-01 00:00:00 UTC") #time indexes for aggregtion
    
    ag_zoo <- aggregate(ts_zoo, list(tim_agr), fun, ...)
    
    if(is.regular(ag_zoo, strict = T) == F) {
        if(insert.missing == T){
            ag_zoo <- insert_missing_records(ag_zoo, step)
            print('missing time steps where inserted as NA values')  
        } else {
            warning("new time series is not strictly regular!")
        }
    }
    
    
    return(ag_zoo)
    
}

#--------------------------------------



## =======================
## CML rainfall retrieval:
## =======================

#---------------------------------

# ---- dry-wet classification ----

#---------------------------------

drywet_schleiss <- function (tl,  q = .94, width = 15, align = 'left', returnSD = F,
                             partial = T, na.rm = T, tsh = NULL) {
  
  #' Dry/wet weather classification based on standard deviation moving window
  #' according to Schleiss et al. 2009.
  #' 
  #'@param tl time series with total losses (dB) or specific losses (dB/km)
  #'@param q quantile defining threshold for sds correpsonding to dry and wet weather, i.e.
  #'ratio of dry data points to all data points.
  #'@param width size of moving sd window in data points.
  #'Size of moving window should consider local rainfall climatology. Schleiss et al. (2019)
  #'recommends size corresponding to 15 - 30 min for Paris.
  #'@param align alignment of moving winodw ('left', 'center', 'right')
  #'@param partial logicla or numeric. See zoo.rollapply.
  #'@param na.rm logical. If TRUE then removes NAs when calculating sd within a window
  #'Returns list with sd time series and time series of classifier (F = dry, T = wet)
  
  require(zoo)
  
  tl <- as.zoo(tl)
  
  sd_tl <- rollapply(tl, width = width, sd, by = 1, align = align,
                      partial = partial, na.rm = na.rm)
  
  if (is.null(tsh)) {tsh <- quantile(sd_tl, q, na.rm = T)}
  
  wet <- sd_tl > tsh
  
  if (returnSD == T) {
    out <- list('sd' = sd_tl, 'wet' = wet)
  } else {
    out <- wet
  }
  
  
  return (out)
  
}

#----------------------------------

# ---- baseline identification ----

#----------------------------------

get_baseline <- function(tl, method = 'moving_quantile', ...) {
  # calculate baseline from total loss using different methods
  
  if (!(method %in% c('schleiss', 'fenicia', 'moving_quantile', 'median'))) {
    stop(paste('wrong method value supplied!'))
  }
  
  if (length(which(is.na(tl))) == length(tl)) {
    warning('Baseline set to NA (tl contains NAs only)')
    B <- tl
    return(B)
    stop()
  }  
  
  if (method == 'schleiss') {
    B <- baseline_schleiss(tl, ...)
  }
  
  if (method == 'moving_quantile') {
    B <- baseline_Qsmoothing (tl, q = .5, win = 7 * 24, aggFun = 'mean')
  }
  
  if (method == 'fenicia') {
    B <- baseline_fenicia(tl, m = 3e-3)
  }
  
  if (method == 'median') {
    B <- baseline_median(tl)
  }
  
  return(B)

}

#--------------------------------------


baseline_median <- function (tl, ...) {
  # Calculate constant baseline using median
  B <- tl
  B[] <- rep(median(tl, na.rm = T), length(tl))
  return(B)
}

#--------------------------------------

baseline_fenicia <- function(tl, m = 3e-3){
    
    ## Baseline model for CML rainfall estimation
    ## for details see eq. 3 and 17 in Fenicia et al 2012, Microwave links
    ## for rainfall estimation in an urban environment: 
    ## Insights from an experimental setup in Luxembourg-City)
    ## Arguments: tl - vector with total path loss [dB]
    ##         m   - filter parameter [-]
    ## Returns: vector with baseline values (dB)
    
  
  # identify NA values
    na.ids <- which(is.na(tl))
    if(length(na.ids) == length(tl)){
        stop('tl must not have only NA values!')
    }
  
  # filter out NA values
    B1 <- rep(NA, length(tl))
    
    if(length(na.ids) == 0){
        B <- tl
    }else{
        B <- tl[-na.ids]
    }
    
  # estimate baseline
    for(i in 1:length(B)){
        if(i == 1){
            b0 <- B[1]
        }else{
            B[i] <- min((1-m)*b0 + m*B[i], B[i])
            b0 <- B[i]
        }
    }
    
  # return back intitally filtered NA values
    if(length(na.ids) == 0){
        B1 <- B
    }else{
        B1[-na.ids] <- B    
    }
    
    return(B1)    
} 

#--------------------------------------

baseline_schleiss <- function (tl, wet = NULL, win = 6 * 3600,
                               approxMethod = "linear", ...) {
    
    ## Baseline estimation algorithm of Schleiss:
    ## Coded by Marc Schleiss, EPFL-LTE, 14th May 2013
    
    ## Arguments:
    ## tl = time series with CML attenuations measurements (in dB) corresponding to tim
    ## wet = state vector (0=dry ; 1=rainy) corresponding to tim
    ## w    = safety window (in seconds) before the antennas are declared dry
    
    ## Output:
    ## tabB = vector with baseline attenuations (in dB) corresponding to tim
    
  if (is.null(wet)) {
    warning('dry/wet classification state vector is missing and will be calculated using drywet_schleiss function with default parameters')  
    wet <- drywet_schleiss(tl)     
  }
    
    
  ## get time as numeric vector
    tl0 <- tl
    tim <- as.numeric(index(tl))
    tl <- coredata(tl)
    ## Local variables:
    Ntim <- length(tim)
    Ntl <- length(tl)
    Nwet <- length(wet)
    notNA <- which(!is.na(tl))
    tabB  <- rep(NA, Ntim)
    
    ## Basic checks on input parameters:
    if (Ntim == 0) {stop("tim is empty")}
    if (Ntim != Ntl) {stop("tim and tl must have the same number of elements")}
    if (Ntim != Nwet) {stop("tim and wet must have the same number of elements")}
    if (any(is.na(tim))) {stop("NA values are not allowed in tim")}
    if (any(tl < 0,na.rm = TRUE)) {warning("there were negative attenuation values in tl")}
    if (is.na(win)) {stop("NA value not allowed in window size")}
    if (win < 0) {stop("negative values are not allowed for window size")}
    
    ## define dry and wet periods
    id.dry <- intersect(which(wet == 0), notNA)
    id.wet <- intersect(which(wet == 1), notNA)
    Ndry   <- length(id.dry)
    Nwet   <- length(id.wet)
    #print("OK 1")
    ## Determine the time to the previous wet periods
    time2prev.wet <- rep(NA, Ntim)
    if (Ndry > 0 && Nwet > 0) {
        for(i in id.dry){
            J <- which(id.wet < i)
            nJ <- length(J)
            if (nJ == 0) {next}
            time2prev.wet[i] <- tim[i] - tim[id.wet[J[nJ]]]
        }
    }
    #print("OK 2")
    
    ## Estimate baseline attenuation during periods with wet antennas
    ## The baseline is obtained by linearly interpolating the attenuation values during periods with dry antennas
    id.dry.antenna <- which(time2prev.wet > win)
    #print("OK 3")
    if (length(id.dry.antenna) == 0) {stop("could not find any period with dry antenna")}
    #print("OK 4")
    tabB[id.dry.antenna] <- tl[id.dry.antenna]
    #print("OK 5")
    id.wet.antenna <- setdiff(notNA, id.dry.antenna)
    #print("OK 6")
    if (length(id.wet.antenna) > 0) {
        x <- tim[id.dry.antenna] / 3600
        y <- tl[id.dry.antenna]
        xout <- tim[id.wet.antenna] / 3600
        tabB[id.wet.antenna] <- approx(x = x, y = y,xout = xout,
                                       method = approxMethod)$y
    }
    
    B <- zoo(tabB, index(tl0))
    
    return(B)    
}

#--------------------------------------

baseline_Qsmoothing <- function (tl, q = .5, win = 7 * 24,
                                 aggFun = 'mean', ...) {
  #' Estimate basilne from sub-hurly total losses using moving quantile window.
  
  ## window range (in each step). 
  #
  #'@param tl - zoo series of tl, one channel
  #'@param q - quantile (0-1) of tl hourly subset (within moving window) used for
  #            basline.
  #'@param win - smoothing window size in hours (default is one week)
  #'@return zoo time series corresponding to input TL time series with baseline values 
    
  #'@details The function uses aggregation to hourly time step to reduce computational cost.  
  #' The window length should be set up considering typical duratio of rain
  #' events to ensure sufficienlty high ratio of dry weather records is within
  
  
  tl_hourly <- zoo_aggreg_by(tl, 60, align = 'right', fun = aggFun, na.rm = T)
  b_hourly <- rollapply(tl_hourly, win, quantile, probs = q, na.rm = T,
                   align = 'center', by = 1, partial = T)
  b_hourly2 <- rollapply(b_hourly, win, mean, na.rm = T, align = 'center',
                    by = 1, partial = T)
  b <- tl
  b[] <- NA
  b[index(b_hourly2)] <- b_hourly2
  b <- na.approx(b, method = 'linear', rule = 2, f = .5)
  
  return (b)
  
}

#--------------------------------------

# ---- Wet antenna correction ---------

#--------------------------------------


get_WAA <- function (A, method = 'constant', ...) {
  # wrapper FUnction returning wet antenna attenuation
  
  if (!(method %in% c('kharadly_1', 'kharadly_2', 'leijense',
                      'schleiss', 'constant', 'pastorek',
                      'kharadly-pastorek'))) {
    stop(paste('wrong method value supplied!'))
  }
  
  if (method == 'schleiss') {
    WAA <- WAA_schleiss(A, ...)
  }
  
  if (method == 'leijense') {
    WAA <- WAA_leijense(A, ...)
  }
  
  if (method == 'kharadly-pastorek') {
    WAA <- WAA_kharadly_pastorek(A, ...)
  }
  
  if (method == 'kharadly_1') {
    WAA <- WAA_kharadly_1(A, ...)
  }
  
  if (method == 'kharadly_2') {
    WAA <- WAA_kharadly_2(A, ...)
  }
  
  if (method == 'constant') {
    WAA <- WAA_constant(A, ...)
  }
  
  if (method == 'pastorek') {
    WAA <- waa_pastorek(A, ...)
  }  
  return(WAA)
}


#--------------------------------------

WAA_constant <- function(A, p = 1.5, ...) {
  # Get constant WAA
  return(zoo(p, index(A)))
}


#--------------------------------------

WAA_kharadly_pastorek <- function(A, L = 1, p = c(6, .125), ...){
  ## Single-frequency Wet antenna attenuation model (Kharadly 2001)
  ## Ipnuts: 
  ##       A   - attenuation of CML for which WAA is estimated
  ##       p    - vector with model parameters c(C, d)
  ## Returns: Aw - wet antenna attenuation of one antenna
  
  Cc <- p[1]   #maximal expected attenuation caused by WA effect
  d <- p[2] #empirical parameter
  
  #filter out negative values
  A[which(A < 0)] <- 0
  
  Aw <- Cc*(1 - exp(-d * A / L))
  
  
  return(Aw)
  
}

#--------------------------------------

WAA_kharadly_1 <- function(A, p = c(6, .125), ...){
  ## Single-frequency Wet antenna attenuation model (Kharadly 2001)
  ## Ipnuts: 
  ##       A   - attenuation of CML for which WAA is estimated
  ##       p    - vector with model parameters c(C, d)
  ## Returns: Aw - wet antenna attenuation of one antenna
  
  Cc <- p[1]   #maximal expected attenuation caused by WA effect
  d <- p[2] #empirical parameter
  
  # filter out negative values
  A[which(A < 0)] <- 0
  
  Aw <- Cc*(1 - exp(-d * A))
  
  
  return(Aw)
  
}

#--------------------------------------
WAA_kharadly_2 <- function(A, A2, p, ...){
  ## Dual-frequency Wet antenna attenuation model (Kharadly 2001)
  ## Ipnuts: 
  ##       A   - attenuation of CML for which WAA is estimated
  ##       A2  - attenuation of CML of same (similar) path and
  ##                    different frequency
  ##       p    - vector with model parameters c(Sp, gamma), where Sp
  ##                    is ratio of path attenuations (based on ITU) and gamma
  ##                    ratio of WAA
  ## Returns: Aw - wet antenna attenuation of one antenna
  
  Sp <- p[1]   #maximal expected attenuation caused by WA effect
  gamma <- p[2] #empirical parameter
  
  #filter out negative values
  A[which(A < 0)] <- 0
  
  Aw <- (A2 - Sp * A) / (gamma - Sp)
  
  return(Aw)
}

#--------------------------------------

WAA_leijense <-function(R, fr, thickpars, refra, ...){
  
  ## Calculates wet antena attenuation according to Leijnse 2007
  ## Last update: 2015/07/14
  ##
  ## Arguments:
  ##      R           - vector with rain rates [mm/h]
  ##      fr          - NWL frequency [GHz]
  ##      thickpars   - vector with gamma and delta parameter to calculate
  ##                     thickness of a water film from rain rate
  ##          refra   - refractive index of water, vector of two elements 
  ##                      with real and imaginary part of a complex number
  ## Returns:         - Wet antenna attenuation [dB]
  
  
  if(missing(thickpars)){
    gamma <- 2.06*10^-5
    delta <- 0.24
  }else{
    gamma <- thickpars[1]
    delta <- thickpars[2]
  }
  
  if(missing(refra)){
    refra <- fun_refra(fr)
  }
  
  
  if(length(which(R<0)) > 0){R[which(R<0)] <- 0}
  
  R0 <- R[which(R==0)]
  R1 <- R[which(R > 0)]
  
  if(length(R1) > 0){
    
    f<-fr*10^9  #frequency of the link [Hz]
    
    #Thicknes of the drop layer on the antenna
    l <- gamma*R1^delta
    
    #Calculate wet antenna attenuation
    
    j <- complex(imaginary=1)
    c <- 2.99*10^8
    expo <- j*2*pi*f/c
    
    
    #Segelstein, D., 1981: "The Complex Refractive Index of Water", M.S.Thesis,University of Missouri--Kansas City
    mh2o<-complex(real=refra[1],imaginary=refra[2])     
    mair <- 1
    mant <- 1.73+0.014*j
    lant <- 0.001
    
    
    x1 <- (mair+mh2o)*(mh2o+mant)*(mant+mair)*exp(-expo*(mant*lant+mh2o*l))
    x2 <- (mair-mh2o)*(mh2o-mant)*(mant+mair)*exp(-expo*(mant*lant-mh2o*l))
    x3 <- (mair+mh2o)*(mh2o-mant)*(mant-mair)*exp(expo*(mant*lant-mh2o*l))
    x4 <- (mair-mh2o)*(mh2o+mant)*(mant-mair)*exp(expo*(mant*lant+mh2o*l))
    
    y1 <- (mair+mant)^2*exp(-expo*mant*lant)
    y2 <- (-(mair-mant)^2*exp(expo*mant*lant))
    
    frac <- (x1+x2+x3+x4)/(2*mh2o*(y1+y2))
    Aw1 <- 10*log(abs(frac)^2,10)
    
    if(length(R0) > 0){
      Aw <- rep(NA, length(R0) + length(R1))
      Aw[which(R==0)] <- 0
      Aw[which(R>0)] <- Aw1
    }else{Aw <- Aw1}        
    
  }else{
    if(length(R0) > 0){Aw <- rep(0, length(R0))}else{Aw <- NULL}
  }
  
  
  return(Aw)
}

#--------------------------------------

WAA_schleiss <- function (A, wet = NULL, tauW = 15, Wmax = 2.3, w0 = 0, ...) {
  
  ## Dynamic wet-antenna attenuation model
  ## Coded by Marc Schleiss, EPFL-LTE, 14th May 2013
  ## References: "Quantification and modeling of wet-antenna attenuation for commercial microwave links"
  ## by Schleiss, M., J. Rieckermann and A. Berne, IEEE Geosci. Remote Sens. Lett., in press.
  
  ## Arguments:
  ## A = time series with CML attenuations measurements (in dB) corresponding to tim, after removal of the baseline
  ## wet = state vector (0=dry ; 1=rainy) corresponding to tim
  ## tauW = average antenna wetting time (in minutes)
  ## Wmax = maximum wet-antenna attenuation (in dB)
  ## w0   = initial wet-antenna attenuation (in dB)
  
  ## Output:
  ## get time as numeric vector
  if (is.null(wet)) {
    warning('dry/wet classification state vector is missing and will be calculated using drywet_schleiss function with default parameters')  
    wet <- drywet_schleiss(A)     
  }
  
  tim <- as.numeric(index(A))
  
  ## Aw = vector with wet-antenna attenuations (in dB) corresponding to tim
  print(paste("execution started at", Sys.time()))
  ## Local variables:
  Ntim <- length(tim)
  NAtt <- length(A)
  Nwet <- length(wet)
  notNA <- which(!is.na(A))
  NNA   <- length(notNA)
  
  ## Basic checks on input parameters:
  if (Ntim == 0) {stop("tim is empty")}
  if (Ntim != NAtt) {stop("tim and A must have the same number of elements")}
  if (Ntim != Nwet) {stop("tim and wet must have the same number of elements")}
  if (any(is.na(tim))) {stop("NA values are not allowed in tim")}
  if (is.na(tauW)) {stop("NA value not allowed for tauW")}
  if (is.na(Wmax)) {stop("NA value not allowed for Wmax")}
  if (tauW <= 0) {stop("tauW must be strictly positive")}
  if (Wmax < 0) {stop("negative values are not allowed for Wmax")}
  if (any(A < 0,na.rm = TRUE)){warning("there were negative attenuation values in A")}
  
  ## Compute wet-antenna attenuation
  Aw <- rep(NA, Ntim)
  if (NNA > 0) {Aw[notNA[1]] <- w0}
  if (NNA == 1) {return(Aw)}
  
  wet.NNA <- wet[notNA]
  tim.NNA <- tim[notNA]
  dt <- tim.NNA[-1] - tim.NNA[-length(tim.NNA)]
  
  for (itr in 2 : NNA) {
    j <- notNA[itr]
    if (wet[j] == 0){
      Aw[j] <- max(min(A[j], Wmax), 0)
      next
    }
    i <- notNA[itr - 1]
    dt <- (tim[j] - tim[i]) / 60
    Aw[j] <- Aw[i] + (Wmax - Aw[i]) * 3 * dt / tauW
    if(Aw[j] > Wmax){Aw[j] <- Wmax}
    if(Aw[j] > A[j]){Aw[j] <- A[j]}
    if(Aw[j] < 0){Aw[j] <- 0}
  }
  print(paste("execution endeded at", Sys.time()))
  
  Aw <- zoo(Aw, index(A))
  return(Aw)
}

#----------------------------

waa_pastorek <- function (A, L, p = c(1.5, 0.5), ...) {
  
  Ai <- A
  Ai[Ai < 0] <- 0
  for (i in 1 : 2) {
    Aw <- 2 * p[1] * (Ai/L)^p[2]
    Ai <- Ai - Aw
    Ai[Ai < 0] <- 0
  }
  return(Aw)
}



#-----------------------------------

# ---- specific attenuation  ----

#-----------------------------------

get_specificAtt <- function (tl, B, Aw, L, wet = NULL) {
  # get specific raindrop path attenuation
  # Arguments:
  # tl - numeric vector of total losses (dB)
  # B - numeric vector of baseline values (dB)
  # Aw - wet antenna attenuation
  # L - length of the CML (km)
  # wet - state vector with dry/wet time steps (False) wet (True)
  
  if (is.null(wet)) {wet <- T}
  
  k <- (tl - B - Aw) * wet / L
  k[k < 0] <- 0
  
  return(k)
}


#-----------------------------------

# ---- attenuation to rain rate ----

#-----------------------------------

get_ITU_pars <- function(freq, pol, conv = T, digits = 3){
    ## returns R-k power law parameters for given frequency and polarization. 
    ## Parameters are based on Rec. ITU-R P.838-3.
    ##
    ## Arguments:   freq  - Frequency in GHz 
    ##           pol   - polarization ("V" or "H")
    ##           conv - indicating if to convert original ITU parameters to 
    ##                  (alpha & beta) parameters for rainfall estimation or keep
    ##                  original ITU values (a & b)
    ##          digits - integer indicating number of decimal places to round the parameters, value NULL does no rounding
    ##
    ## Returns:  data frame with R-k power law parameters.
    
    
    #check inputs  
    if(length(freq) != length(pol)){stop("frequency vector does not have same length as polarization vector")}
    if(length(pol) != length(which(is.na(match(pol, c("V", "H")))==F))){
        stop("Polarizaton has to be character vector with only V or H symbol for vert. resp. horizon. polarization")
    }
    if(min(freq) < 1 || max(freq) > 1000){stop("Frequency is out of ITU rec. range")}
    
    #load ITU data (can be later extended to 1 GHz)
    Fr <- c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,120,150,200,300,400,500,600,700,800,900,1000)
    alphaH <- c(54069.00,18815.58,6585.12,2888.68,1349.47,595.63, 301.99,192.83,144.51,116.29,95.95,68.41,51.96,41.09,33.35,27.72,23.55,20.39,17.89,15.87,14.19,12.78,11.56,10.51,9.6,8.79,8.08,7.45,6.88,6.38,5.92,5.51,5.14,4.8,4.5,4.22,3.96,3.73,3.52,3.32,3.14,2.98,2.83,2.69,2.56,2.44,2.33,2.22,2.13,2.04,1.95,1.88,1.8,1.74,1.67,1.61,1.56,1.51,1.46,1.41,1.37,1.33,1.29,1.25,1.22,1.18,1.15,1.12,1.1,1.07,1.05,1.02,1,0.98,0.96,0.94,0.92,0.9,0.89,0.87,0.86,0.84,0.83,0.81,0.8,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.68,0.67,0.66,0.66,0.65,0.64,0.64,0.63,0.55,0.49,0.46,0.46,0.48,0.5,0.52,0.54,0.57,0.59,0.6)
    betaH <- c(1.03,0.98,0.94,0.89,0.81,0.70,0.62,0.59,0.59,0.61,0.63,0.68,0.72,0.76,0.80,0.82,0.85,0.86,0.88,0.89,0.9,0.91,0.92,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.06,1.08,1.09,1.1,1.11,1.12,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.21,1.22,1.23,1.24,1.24,1.25,1.26,1.27,1.27,1.28,1.29,1.29,1.3,1.31,1.31,1.32,1.32,1.33,1.34,1.34,1.35,1.35,1.36,1.36,1.37,1.37,1.38,1.38,1.38,1.39,1.39,1.4,1.4,1.41,1.41,1.41,1.42,1.42,1.42,1.43,1.43,1.43,1.44,1.44,1.44,1.45,1.45,1.45,1.45,1.46,1.46,1.46,1.46,1.47,1.51,1.54,1.57,1.59,1.6,1.6,1.6,1.59,1.58,1.57,1.56)
    alphaV <- c(178137.64,54317.86,16439.15,6340.75,2970.44,1540.13,781.02,393.46,229.07,161.39,127.54,85.17,60.89,48.56,39.99,32.85,27.25,23.08,19.97,17.6,15.72,14.19,12.89,11.77,10.79,9.91,9.13,8.43,7.79,7.21,6.69,6.21,5.78,5.38,5.02,4.7,4.39,4.12,3.87,3.64,3.43,3.24,3.06,2.89,2.74,2.61,2.48,2.36,2.25,2.15,2.05,1.97,1.88,1.81,1.74,1.67,1.61,1.55,1.5,1.45,1.4,1.36,1.32,1.28,1.24,1.2,1.17,1.14,1.11,1.08,1.06,1.03,1.01,0.99,0.97,0.95,0.93,0.91,0.89,0.87,0.86,0.84,0.83,0.82,0.8,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.68,0.67,0.66,0.65,0.65,0.64,0.64,0.63,0.55,0.49,0.46,0.46,0.48,0.5,0.53,0.55,0.57,0.59,0.6)
    betaV  <- c(1.16,1.12,1.05,0.99,0.94,0.88,0.80,0.71,0.65,0.63,0.64,0.68,0.72,0.78,0.82,0.86,0.89,0.92,0.94,0.96,0.97,0.99,1,1.01,1.02,1.02,1.03,1.04,1.05,1.05,1.06,1.07,1.08,1.09,1.1,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.21,1.22,1.23,1.24,1.25,1.26,1.26,1.27,1.28,1.28,1.29,1.3,1.31,1.31,1.32,1.32,1.33,1.34,1.34,1.35,1.35,1.36,1.36,1.37,1.37,1.38,1.38,1.39,1.39,1.39,1.4,1.4,1.41,1.41,1.41,1.42,1.42,1.42,1.43,1.43,1.43,1.44,1.44,1.44,1.45,1.45,1.45,1.45,1.46,1.46,1.46,1.46,1.47,1.47,1.47,1.47,1.48,1.48,1.51,1.55,1.58,1.6,1.6,1.59,1.59,1.58,1.58,1.57,1.57)
    
    # interpolate parameters to match the link frequencies
    cml.pars <- matrix(NA, nrow=length(pol), ncol=2)
    colnames(cml.pars) <- c("alpha", "beta")
    
    for(i in 1:length(pol)){
        if(pol[i]=="H"){
            cml.pars[i ,1] <- approx(Fr,  alphaH, freq[i])$y    
            cml.pars[i ,2] <- approx(Fr,  betaH, freq[i])$y  
        }
        
        if(pol[i]=="V"){
            cml.pars[i ,1] <- approx(Fr,  alphaV, freq[i])$y    
            cml.pars[i ,2] <- approx(Fr,  betaV, freq[i])$y  
        }
    }
    
    if(conv == F){
        cml.pars[ ,2] <- 1/cml.pars[ ,2]
        cml.pars[ ,1] <- cml.pars[ ,1]^(-cml.pars[ ,2])
        colnames(cml.pars) <- c("a", "b") 
    }
    
    if (!is.null(digits)) {cml.pars <- round(cml.pars, digits)}
      
    return(as.data.frame(cml.pars))
}

#--------------------------------------

kRmodel <- function(k, p){
  ## function to get rain rate from attenuation
  ## Intputs: k - speicific attenutation [dB/km]
  ##          p - model pamaeters, c("alpha", "beta")
  ## Output: r.mod - modeled rainfall intensity [mm/h]
  
  k[k < 0] <- 0
  r.mod <- p[1] * k ^ p[2]
  return(r.mod)
}

#--------------------------------------


## =========================
## Optimize model parameters
## =========================

fit_kRmodel <- function(k, R, alpha.lim, beta.lim, logtransform = F)
    # function to fit power law k-R model: R = alpha*k^beta;
    # Arguments:
    # k - vector with specific attenuation
    # R - vector with reference rain rates
    # alpha.lim, beta.lim - limits of slope (alpha) and power (beta) parameters
    # and their initial values: c(min, max, ini))
    # logtransform - logaritmic transformation to reduce heteroscedasticity of
    # data
{    
    
    #p.ini <- c(mean(p1.lim), mean(p2.lim)) #
    k <- as.numeric(k)
    R <- as.numeric(R)
    p.ini <- c(alpha.lim[3], beta.lim[3])
    p <- optim(par = p.ini, fn = minfun_kRmodel, gr = NULL,
               lower = c(alpha.lim[1], beta.lim[1]),        
               upper = c(alpha.lim[2], beta.lim[2]),       
               method = "L-BFGS-B",
               k = k,
               r = R,
               logtransform)$par
    
    return(p)
}

#--------------------------------------

minfun_kRmodel <- function(par, k, r, logtransform = F)
    #function, which is minimized in fitting
    # here: LS criterion
    # alpha <- pars[1]
    # beta <- pars[2]
{    
    r[which(r < 0)] <- 0
    r.mod <- kRmodel(k, par)
    cost.val <- costfun_kRmodel(r.mod, r, logtransform)
    return(cost.val)
    
}

#--------------------------------------

costfun_kRmodel <- function(r.mod, r, logtransform = F){
    
    if (logtransform == F) {
        cost <- sum((r.mod - r)^2, na.rm=T)
    } else if (logtransform == T) {
        cost <- sum((log(r.mod + 1) - log(r + 1))^2, na.rm = T)
    } else {
        stop ('logtransform is not logical!')
    }
    
    return (cost)
}

#--------------------------------------

      

## ==============
## CML adjusting
## ==============

# Functions for adjusting CMLs as proposed in Fencl et al. (2017), 
# Gauge-adjusted rainfall estimates from commercial microwave links, HESS

dyn_cal_power <- function(A, R, L, beta, w, r.thr, al.lim, aw.lim, b.lim, wei)
    ## function to dynamicaly calibrate simpified k-R relation (R=alpha*k^beta-aw)
    ## alpha and aw are calibratio parameters.    
    ## Arguments: 
    ##  A - rain attenuation vector [dB] (aggregated to arbitrary time step) 
    ##  R - rain rate vector [mm/h] corresponding to A 
    ##  L - CML length [m]
    ##  beta - beta parameter (e.g. taken form ITU) 
    ##  w - length of calibration window  (number of time steps)
    ##  r.thr - minimal rain rate to include time step for calibration
    ##  al.lim, aw.lim, b.lim - vector with upper, lower limit and initial value
    ##  for alpha resp. aw, resp. beta
    ##  wei - weights binding the value close to initial one (vector of three
    ##        elements, only third elements intended for beta is currently
    ##        implemented)
    ##  Returns:
    ##  p 

{
    #simple argument test
    if(length(A)!=length(R)){stop("A and R must have same length")}
    if(missing(L)==T){
        warning("L is missing, A is considered to correspond to specific attenuation!")
        L <- 1000
    }
    
    if(missing(beta)==T){
        warning("beta is missing! beta is set to unity!")
        beta <- 1
    }
    
    if(missing(w)==T){
        warning("w is missing, w is set to 3!")
        w <- 3
    }
    
    if(missing(r.thr)==T){
        warning("r.thr is missing, no time steps will be excluded from optimization!")
        r.thr <- 0
    }
    
    if(missing(al.lim)==T){
        warning("al.lim is missing, al.lim is set to c(0,15)!")
        al.lim <- c(0,15)
    }
    
    if(missing(aw.lim)==T){
        warning("aw.lim is missing, aw.lim is set to c(0,10)!")
        aw.lim <- c(0,10)
    }
    
    
    #treat NA values
    k0 <- 1000*A/L
    k0[which(k0 < 0)] <- 0
    #k0 <- k0^beta
    
    id.na <- which(is.na(k0+R)==T)
    if(length(id.na)== 0){id.na <- 10^8}
    
    r1 <- R[-id.na]
    k1 <- k0[-id.na]
    
    p <- matrix(NA, length(r1), 3)
    
    if(w > length(which(k1 >= r.thr))){
        warning("w is longer than calibration set, w is set to length of calibration set!")
        w <- length(which(r1 >= r.thr))
    }
    
    #par estimation for each time step
    id.out <- which(r1 < r.thr)  #rain rates lower than threshold
    if(length(id.out)== 0){id.out <- 10^8}
    r.w1 <- r1[-id.out][1:w]
    k.w1 <- k1[-id.out][1:w]
    first1 <- seq(1, length(r1), 1)[-id.out][w]  #first  
    
    p.j <- fit_kR_power(k.w1, r.w1, al.lim, aw.lim, b.lim, wei) #6000/cml.info2$length[i]))
    #f.nam <- paste("opt", zero_before(first1,3), ".png", sep="")
    #plot_costF(k.w1, r.w1, al.lim, aw.lim, 20, p.j, f.nam)
    
    p[first1, ] <- p.j
    
    for(j in (first1+1):length(r1)){
        if(r1[j] < r.thr){   #if rain rate which was excluded (r < r.thr) next
            p[j, ] <- p.j
        }else{
            #fill calibration window with next value
            r.w1 <- c(r.w1[-1], r1[j])
            k.w1 <- c(k.w1[-1], k1[j])
            p.j <- fit_kR_power(k.w1, r.w1, al.lim, aw.lim, b.lim, wei) #6000/cml.info2$length[i]))
            #f.nam <- paste("opt", zero_before(j,3), ".png", sep="")
            #plot_costF(k.w1, r.w1, al.lim, aw.lim, 20, p.j, f.nam)
            p[j, ] <- p.j
        }
        
    }
    
    #set parameters of initial time steps to value equal to first estimated pair of parameters
    for(j in 1:first1){p[j, ] <- p[first1,]}
    p0 <- matrix(NA, length(k0), 3)
    p0[-id.na, ] <- p
    
    return(p0)       
}

#--------------------------------------

fit_kR_power <- function(k, ref, p1.lim, p2.lim, p3.lim, wei)
    #function to fit power law k-R model: R = alpha*(k - kw)^beta;
    # where kw = specific wet antenna attenuation
    
    # k - vector with specific attenuation
    # ref - vector with reference rain rates
    # p1.lim - limits of slope (alpha) parameter - c(min, max, ini))
    # p2.lim - limits of intersect (kw) parameter - c(min, max, ini))
    # p3.lim - limits of power (beta) - c(min, max, ini))
    # wei
    
    
{    
    
    #p.ini <- c(mean(p1.lim), mean(p2.lim)) #
    p.ini <- c(p1.lim[3], p2.lim[3], p3.lim[3])
    
    p <- optim(par= p.ini, fn=minfun_power,
               lower= c(p1.lim[1], p2.lim[1], p3.lim[1]),        
               upper= c(p1.lim[2], p2.lim[2], p3.lim[2]),       
               method= "L-BFGS-B",
               x = k,
               r = ref,
               w = wei,
               beta0=p3.lim[3])$par
    
    return(p)
}

#--------------------------------------

minfun_power <- function(pars, x, r.ref, w, beta0)
    #function, which is minimized in fitting
    # here: LS criterion
    #alpha <- pars[1]
    #kw  <- pars[2]
    #beta <- pars[3]
    #w   - weight to bind beta close to its initial value (beta0)
    #beta0 - initial value of beta
{    
    r.ref[which(r.ref < 0)] <- 0
    r.mod <- kRmodel_power(x, pars)
    cost.val <- cost_fun_power(r.mod, r.ref, pars, w, beta0)
    return(cost.val)
    
}

#--------------------------------------

kRmodel_power <- function(k, p){
    ## function to get rain rate from attenuation
    ## Intputs: k - speicific attenutation [dB/km]
    ##          p - model pamaeters, c("alpha", "kw", "beta")
    ## Output: r.mod - modeled rainfall intensity [mm/h]
    
    if(class(p)=="numeric" || class(p)=="integer"){p <- matrix(p, 1,3)}
    kw <- (k-p[, 2])*((k-p[, 2]) >= 0)  #correct for WAA and set negative values to zero
    r.mod <- p[ ,1]*kw^p[ ,3]
    r.mod[which(r.mod < 0)] <- 0
    return(r.mod)
}

#--------------------------------------

cost_fun_power <- function(mod, ref, p, w, beta0){
    sqrt(sum((mod - ref)^2, na.rm=T)/length(ref))/max(c(ref,1), na.rm=T) + w[3]*(p[3]-beta0)^2
}

#--------------------------------------

disaggreg_parameters <- function(p.mtx, p.idx, pout.idx, method="constant")
    
    ## function to disaggregate parameter arrays or matrixes
    ## (function was designerd for dyn_cal  disaggregation
    ## input: p.mtx   - matrix with parameters
    ##        p.idx - time stamps (or other index) coresponding to parameters (p.mtx rows)
    ##        pout.idx - time stamps to which p.mtx should be disaggregated
    ##        method - method to use for extrapolation: "constant" or "linear" (see ?approx)
    ## Returns: pmtx.out - table with dissagregated parameters with each row corresponding to pout.idx
    
{
    
    #enable proceed p.mtx in a numeric vector class
    if(is.null(nrow(p.mtx))){
        warning("p.mtx is a vector")
        n.col <- 1
        n.row <- length(p.mtx)
    }else{
        n.col <- ncol(p.mtx)
        n.row <- nrow(p.mtx)
    }    
    
    ## basic parameter check
    if(length(p.idx) != n.row){
        stop("length of input index vector(pin.idx) do not
             match number of rows in the parameter tabele (p.mtx)!")
    }
    
    
    ## merge and extrapolate
    p0 <- zoo(p.mtx, p.idx)
    p.out <- zoo(rep(0, length(pout.idx)), pout.idx)
    pmtx.out <- zoo(matrix(NA, length(pout.idx), n.col), pout.idx)
    
    p.mgd <- merge.zoo(p.out, p0)
    t.mgd <- index(p.mgd)
    
    for(i in 2:ncol(p.mgd)){
        
        x <- apply(p.mgd[ ,c(1,i)], 1 ,sum)
        x <- zoo(x, t.mgd)
        
        #NA check 
        if(length(which(is.na(x))==T) == length(x)){
            warning(paste("p.mtx[", i-1, ", ] do not match any index. NAs itroduced!", sep=""))
        }else{
            x <- na.approx(x, method=method, rule=2, f=1)
            pmtx.out[ ,i-1] <- x[pout.idx]
        }
    }
    
    return(pmtx.out)
    }

#--------------------------------------



## ===================================
## CML RAINFALL SPATIAL RECONSTRUCTION 
## ===================================


cml.to.points <- function(c.tab, Lmax){
    
    #------
    # Schematize CMLs as single points (in their centre) or as multiple points
    # equaly spaced along the CML path
    #
    # Arguments:    c.tab - data.frame with CML cartesian coordinates
    #                    (Ax, Ay, Bx, By) [m]
    #            Lmax - maximal length [m] of segment represented by one point
    # Returns:   list containg: K - table with point coordinates and CML id to
    #            which they belong to, K.dist - square matrix with distances
    #            between points (cols and rows coorespond to rows of K table)
    #         
    # Last modified: 2017/08/14 
    #------
    
    # CML position length
    xA <- c.tab[ ,1]
    yA <- c.tab[ ,2]
    xB <- c.tab[ ,3]
    yB <- c.tab[ ,4]
    Lcml <- sqrt((xA - xB)^2 + (yA - yB)^2) # CML length
    
    # calculated number of points representing each CML
    
    NL <- length(c.tab[,1]) # Number of CMLs 
    NK <- ceiling(Lcml/Lmax) # vector with count of points representing each CML
    
    # define matrix to store info about points
    
    K <- matrix(NA, ncol=3, nrow=sum(NK))
    colnames(K) <- c("X", "Y", "cml_id")
    K <- as.data.frame(K)
    K$cml_id <-  rep(1:NL, times = NK) # assign CML ids to each point
    # calculate positions of points
    
    for (i in 1:NL){
        id <- which(K$cml_id == i) #(sum(NK[1:i-1])+1):(sum(NK[1:i])) 
        # calculate X coordinates
        K$X[id] <- c(xA[i] + seq(from = 1, by = 2, length = NK[i])
                     * (xB[i] - xA[i])/(2 * NK[i])) 
        # calculate Y coordinates
        K$Y[id] <- c(yA[i] + seq(from = 1, by = 2, length = NK[i])
                     * (yB[i] - yA[i])/(2 * NK[i])) 
    }
    
    return(K)
    
}

#--------------------------------------

K.distance.mtx <- function(K){
    
    #------
    # function to calculate distances between CML points 
    #
    # Arguments:    K - table with coordinates of points representing CML as 
    #                returned by fun cml.to.points
    # Returns:   K.dist - square matrix with distances between points,
    #                     cols and rows coorespond to rows of K table
    #         
    # Last modified: 2017/08/14 
    #------
    
    K.dist <- matrix(NA, ncol=nrow(K), nrow=nrow(K))
    for (i in 1:nrow(K)){
        K.dist[1:nrow(K),i] <- sqrt((K[i,2] - K[1:nrow(K),2])^2 +
                                        (K[i,3] - K[1:nrow(K),3])^2)
    }
    
    return(K.dist)
}

#--------------------------------------

distribute.cmlR.1D <- function(R, K, K.dist, q.var, z=5, n.iter = 10, rad=3000,
                               tune = F)
{
    
    #------
    # Distribute path-averaged CML rainfall along points representing CMLs
    #
    # Arguments:   
    # R  - path averaged rainfall intensity of CMLs
    # K     - table with coordinates of points representing CML and id of CML to
    #         which the points belong to
    # K.dist - square matrix with distances between points
    # q.var  - varinace due to CML quantization (and WAA) for each CML
    # z - coefficient adjusting ration between distance and q.var
    #     weights
    # n.iter - number of iterations when esitmating R distribution
    # rad - radius of influence (decorrelarion distance) [m]
    # tune - if tune = T then results for each iteration are returned in a form
    #        of matrix else final iterated rainfall is
    #                   returnd in a form of K table
    # Returns:   K - K table with two more columns corresponding to initial and
    #                esitmated (distributed) rainfall intensities
    #            OR (if tune = T)
    #            rain - matrix with iterated rainfall intensitis where rows
    #                   correspond to itteratons and columns to CML points     
    #         
    # Last modified: 2017/08/14 
    #------
    
    # 
    NP <- nrow(K)  # number of points
    NL <- nlevels(as.factor(K$cml_id)) #number of CMLs
    
    
    # transform varaibles belonging to CMLs to match CML points (K table)
    q.var <- q.var[K$cml_id]
    R <-  R[K$cml_id]
    
    # add to K rable columns with path-averaged and distributed rainfall
    # distrib. rainfall is initialy set as initial rainfall
    K <- cbind(K, data.frame('R'= R, 'rij' = R)) 
    
    # calculate weighting matrix (accounting for distance + CML accuracy)
    W <- z/((K.dist/1000)^2)
    + 1/matrix(q.var, ncol = NP, nrow = NP, byrow=T)
    
    #------
    # Estimate rainfall distribution
    
    rain <- matrix(NA, ncol=NP, nrow=n.iter) # mtx for est. R in each iteration
    
    for(loop in 1:n.iter){ 
        # loop 20 times to iterate areal rainfall
        
        for(i in 1:NL){ # loop to list all the links (NL = number of links)
            K.est <- K[which(K$cml_id == i), 4:5]
            K.ref <- K[which(K$cml_id != i), 4:5]
            
            #assign rain rate from K.ref (using distance weighted mean)
            NKi <- nrow(K.est) # number of points belonging to CML 
            
            # loop to estimate rain rate at all points of particular link
            for(j in 1:NKi){                                               
                
                # use only info from points of neighboring links)
                dis <- K.dist[which(K$cml_id != i), which(K$cml_id == i)[j]] 
                
                # take into account only points in a vicinity (3 km)
                ind <- which(dis < rad) 
                W.j <- W[which(K$cml_id != i), which(K$cml_id == i)[j]] # weights
                
                # estimate rain rate
                K.est$rij[j] <- sum(W.j[ind] * K.ref$rij[ind]) / sum(W.j[ind]) 
            }
            
            # optimze CML estimated rain rate so that mean rain rate of all
            # points of blonging to particular CML corresponds to the CML
            # path-averaged rain rate
            
            rij <-  K.est$rij + K.est$R - mean(K.est$rij)
            rij[which(rij < 0)] <- 0 # treat negative rainfall intensities
            K$rij[which(K$cml_id == i)] <- rij # place data into K data.frame
        }
        rain[loop,] <- K$rij # "rain" matrix (for tunning the algorithm)
        
    }  #end optimization loop
    
    # return
    if(tune == F){
        return(K)    
    }else{
        return(rain)    
    }
    
}

#--------------------------------------

Rextrapolate_to_2Dgrid <- function(K2, x.lim, y.lim, g.size, q.var, z, rad = 3)                                  {
    
    #------
    # extrapolate rain rates at CML points to the regular 2D grid
    #
    # reconstructs the rainfall distribution from CML path averaged rainfall
    #
    # Arguments:    K2 - data.frame CML point info and rainfall
    #        x.lim - Xcoordinates of field edges
    #        y.lim - Ycoordinates of field edges
    #        g.size - size of a grid cell [m]
    #       q.var  - varinace due to CML quantization (and WAA) for each CML
    #            z - coefficient adjusting ratio between distance and q.var
    #                weights
    # rad - radius of influence (decorrelarion distance) [m]
    # Returns:   R2D - data.frame with two columns with X and Y ccordinates
    #                  and third column with rainfall intensities
    #
    # Last modified: 2017/08/15 by Martin Fencl
    #------
    
    #create grid
    #define size and resolution
    g.st <- c("X"=x.lim[1],"Y"=y.lim[1])              # coordinates of top/left edge
    g.dim <- c(x.lim[2]-x.lim[1], y.lim[2]-y.lim[1])  # dimensions in m
    print(g.dim)
    #g.size <- 500                                     # resolution in m
    
    #x and y coordinates of cell centres
    x.cor <- g.st[1] + seq(from=0,to=g.dim[1]-g.size, by=g.size) + g.size/2
    y.cor <- g.st[2] + seq(from=0,to=g.dim[2]-g.size, by=g.size) + g.size/2
    
    #define grid matrix
    G <- matrix(NA, ncol=3, nrow=length(x.cor)*length(y.cor))
    colnames(G) <- c("X", "Y", "R")
    G <- as.data.frame(G)
    
    #cell centres coordinates
    k <- 0
    #assign to the each point of grid respective coordinates
    for (j in y.cor[length(y.cor):1]){
        for (i in x.cor){
            k <- k + 1
            G[k, 1:2] <- c(i,j)
        }
    }
    
    #assign rain rate to the grid points
    for(i in 1:length(G[,1])){
        
        # Distance weight
        dis <- sqrt((G$X[i]-K2$X)^2 + (G$Y[i]-K2$Y)^2)
        W.d <- round(((1-dis/rad)/(dis/rad))^2, 2) #distance weight
        
        # CML accuracy weight - (inverse variation of noise)
        W.n  <- 1/q.var[K2$cml_id]
        
        # Total weight
        W.j <- z*W.d + W.n        
        W.j[which(dis > rad)] <- 0 # zero Weight to points further than rad
        
        # if grid point is at the same position as data point assign the
        # data point value, lse use weighted mean
        if (length(which(dis == 0)) == 1) {G$R[i]<- K2$rij[which(dis == 0)]  
        }else{
            G$R[i] <- sum(W.j * K2$rij, na.rm = T) / sum(W.j)                            
        }
        
    }
    
    
    return(G)     #data.frame with coordinates and rain rates
    
}  


#--------------------------------------

select.best.link <- function(cell.id, e.link, var1, var2){
    #------
    # select most relevant link from all links belogning to one cell
    #
    #  Arguments:
    # cell.id   - id of cell to which more links belong
    # e.link - data.frame with all link ids and respective cell ids
    # var1   - variance of link due to siganl quantization and wet antenna
    #          (vector with values for all links)
    # var2   - expected variance along the link path due to rainfall spatial
    #          variability
    #          (vector with values for all links)
    # Returns:
    # integer - id of selected link
    #------  
    
    cml.id <- e.link[which(e.link[ ,2] == cell.id), 1] 
    vars <- var1[cml.id] + var2[cml.id]
    vars.id <- which(vars  == min(vars))
    
    return(cml.id[vars.id])
}

#--------------------------------------

R.along.cml.cell <- function(cell, dis0, L2, w.rain=1/5){  
    #------
    # reconstructs the rainfall distribution along particular CMLs by
    # joined analysis of nearby CMLs
    #
    # Arguments:  cell - data.frame with CML info schematized by square cells 
    #          dis0 - square matrix, where [n , m] element represent distances 
    #                 between cell[n, ] and cell[m, ].
    #          L2 - vector with amounts of cells representing respective CML
    #          w.rain - weight reflecting rainfall variablity along the CML
    #                   which is assumed to be proportional to the rainfall
    #                   intensity    
    #         
    # Last modified: 2017/08/14 by Martin Fencl
    #------
    
    # calculate weights
    
    var <- cell[,5] + cell[,6] + (cell[,3]*dis0^2*w.rain)^2 # calculate variance matrix
    
    # ---------------------------
    
    # iterate (20 x) to estimate rainfall distribution
    
    NL <- length(L2)      # Number of CMLs
    
    for(loop in 1:20){                 # loop 20 times to iterate areal rainfall
        
        # -------------
        for (i in 1:NL){                   # loop to list all the links (NL=number of links)
            if(L2[i] <= 1){next}             # if CML is represented by a single cell
            
            id.cml <- which(cell[ ,2] == i) # 
            
            r.est <- cell[id.cml, 3:4]         # rainfall along CML i
            r.ref <- cell[-id.cml, 3:4]        # rainfall along all the other CMLs
            
            
            # Loop over segments of CML i, estimate rainfall intens. for a each segment 
            for(s in 1:length(id.cml)){         
                dis <- dis0[-id.cml, id.cml[s]] # distances to segments of neigboring CMLs
                id.s <- which(dis < 3) # limit the selection to vicinity of 3 km
                
                # if no neighboring segment found increase dist. to 5 km
                if(length(id.s)==0){id.s <- which(dis < 5)}  
                if(length(id.s)==0){next} # if no eighboring segment found skip estimation
                
                # calculate weights
                var.s <- var[-id.cml, id.cml[s]][id.s]      # select values from variance matrix
                W.s <- (var.s*sum(1/var.s))^-1
                r.est[s, 2] <- sum(W.s* r.ref[id.s , 2])   # estimate rain rate
            }
            
            
            # optimze CML estimated rainfall given the constraint that mean rainfall
            # intensity of all segments belonging to particular CML corresponds to the
            # CML's path-averaged rainfall intensiy)
            
            rij <-  r.est[ ,2] + r.est[, 1] - mean(r.est[, 2])
            rij[which(rij < 0)] <- 0                      # treat negative rain rates
            
            
            cell[id.cml, 4] <- rij                        # place data into cell data.frame
        }
        
        # -------------
    }  #end optimization loop
    
    return(cell)
    
}

#--------------------------------------

Rcell2reg.grid <- function (cell, dis0, rad1 = 3, rad2 = 5){
    
    #------
    # extrapolate rain rates along CMLs to the regular 2D grid
    #
    # Arguments:  cell - data.frame with CML info schematized by square cells 
    #          dis0 - distance matrix of all cells
    #      
    #         
    # Last modified: 2017/08/14 by Martin Fencl
    #------
    
    rad1 <- 3  # radius of influence (decorr. distance) [km]
    var <- cell[,5]+cell[,6]+(cell[,4]*dis0^2/5)^2  # calculate variance matrix  
    id.cml <- which(is.na(cell[,2]) == F)               # identify CML cells
    
    for(i in 1:length(cell[,1])){
        
        if(is.na(cell[i, 2])==F){next} # skip cells belonging to CMLs
        
        dis <- dis0[id.cml, i] # distance of cell i to the all CML cells
        
        id.r <- which(dis < rad1) # select cells within radius rad1 (3 km)
        
        # if no CML cell in  increase radius to rad 2 (5 km)
        if(length(id.r)==0){id.r <- which(dis < rad2)}  
        if(length(id.r)==0){next} #if no CML cell within radius rad2, skip
        
        # estimate weights
        var.i <- var[id.cml, i][id.r]   # select values from variance matrix
        W.i <- (var.i*sum(1/var.i))^-1
        cell[i, 4] <- sum(W.i* cell[id.cml[id.r] ,4])   # estimate rain rate
    }
    
    return(cell[ , 4])
}  

#--------------------------------------

get_cml_accuracy <- function(Fr, Pol, length, std = .75) {
    # The accuracy is used as weight in the reconstruction algorithm.
    # Signal quantization + wet antenna effect is assumed to cause normally 
    # distributed error of +- 1.5 dB in estimated baseline (N(0, 0.75)) which is
    # propagated to estimated rain rainfall intensity [mm/h]. Since
    # rainfall is calculated from specific attenuation, longer links will provide
    # more accurate path-averaged rainfall esitmates. Uncertainty due to k-R power
    # law approximation is not assumed here.
    
    p <- get_ITU_pars(Fr, Pol)
    
    var1 <- (p$alpha * (std / length) ^ p$beta)^2
    return(var1)
}

## ================
## Analyze results:
## ================

#------------------------------

# ---- Performance metrics ----

#------------------------------

#--------------------------------------

summary_stats <- function (rRef, rEst) {
  # Returns vector with summary statistics (RMSE, rel. error and R^2)
  # Arguments:
  # Rref - vector with reference rainfall (along link path)
  # Rref - vector reference rainfall (along link path)
  # Returns: vector with summary stats
  # Details: (NA are values are omitted)
  metrics <- c('rmse' = rmse(rRef, rEst),
               'rel_error' = rel_bias(rRef, rEst),
               'R2' = cor(rRef, rEst, use = 'na.or.complete')^2)
  1+1
  return(metrics)
}



rmse <- function(x,y){sqrt(mean((x - y)^2, na.rm=T))}

#--------------------------------------

rel_bias <- function(x,y){
  ## relative bias (rel. error of total cumulative quantities)
  ## Arguments: x - reference
  ##         y - evaluated variable
  ## Returns: relaitve bias
  
  dif <- y - x
  nas <- which(is.na(dif)==T)
  if(length(nas) > 0){
    err <- sum((y-x)[-nas])/sum(x[-nas])    
  }else{
    err <- sum(y-x)/sum(x)    
  }
  return(err)
}

#--------------------------------------

nse <- function(x,y){
  ## Nash-Sutcliffe Efficiency index
  ## Arguments: x - reference
  ##         y - evaluated variable
  ## Returns: Nash-Sutcliffe index
  
  dif <- y - x
  nas <- which(is.na(dif)==T)
  if(length(nas) > 0){
    nse <- 1 - sum(((x-y)[-nas])^2)/sum((x[-nas] - mean(x[-nas]))^2)    
  }else{
    nse <- 1 - sum((x-y)^2)/sum((x - mean(x))^2)    
  }
  return(nse)
}



#-----------------------------

# ---- Plotting functions ----

#-----------------------------

plot2 <- function(..., h = seq(0, 50, 5), v = seq(0, 50, 5), g.col = '#00000044'){
  # plot graph with grid lines
  
  plot(...)
  
  abline(h = h, lty=3, col = g.col) #grid
  abline(v = v, lty=3, col = g.col) #grid
}

#--------------------------------------

polyline <- function(x,y, ..., from='zero', border=NA){
  ## Function to plot time series as solid polygon
  ## Arguments: x, ... - Arguments of the plot function
  ##         from - indicates the baseline of the  polygon, 'min' sets the
  ##                baseline to the minimal value of y, 'zero ' sets it to 0 
  
  if(from == 'zero'){
    poly = data.frame(x[c(1 ,1:length(x), length(x), 1)],
                      c(0, y[1:length(y)], 0, 0))
  }else if(from == 'min'){
    y.min = min(y, na.rm=T)
    poly = data.frame(x[c(1 ,1:length(x), length(x), 1)],
                      c(y.min, y[1:length(y)], y.min, y.min))
  }else{stop('from has to be \'zero\' or \'min\'')}
  
  polygon(poly, ..., border=border)
}

#--------------------------------------

plot_scatter = function(x,y, y.intersp = 2, cex.leg = 1.2, pos = 'topleft', ...){
  ## Function to plot nice scatter plot with legend 
  
  plot2(x,y, ...)
  co = round(cor(x,y, use= 'na.or.complete'), 2)
  l = lm(y ~ x)
  abline(l, col=2)
  l = round(l[[1]], 2)
  legend(pos, c(paste('r =', co), paste('interc. =', l[1]), paste('slope =', l[2])),
         y.intersp = y.intersp, cex = cex.leg)
}

#--------------------------------------

plot_3_scatters = function(z,...){
  ## Function to plot three scatter plots in row showing correlation between
  ## three variables
  ## Arguments: z - data.frame with three variables
  
  layout(matrix(1:3, 1,3))
  par(mar=c(4,4,.5,.5))
  
  plot_scatter(z[ ,1], z[ ,2], xlab = colnames(z)[1], ylab = colnames(z)[2], ...)
  plot_scatter(z[ ,1], z[ ,3], xlab = colnames(z)[1], ylab = colnames(z)[3],...)
  plot_scatter(z[ ,2], z[ ,3], xlab = colnames(z)[2], ylab = colnames(z)[3],...)
}

#--------------------------------------

tim.colors.palette <- function() {
  
  # function returns tim.colors() palette from fields package
  
  return(
    c(
      "#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", 
      "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", 
      "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", 
      "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", 
      "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
      "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", 
      "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", 
      "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", 
      "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", 
      "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
      "#AF0000", "#9F0000", "#8F0000", "#800000"
    )
  )
  
}

#--------------------------------------

zero_before <- function(x, dig){
  ## function to convert numbers to strings and add zeros before them 
  ##
  ## Arguments: x   - integer
  ##         dig - number of digits
  ##
  ## Returns: number as character with zeros before
  
  if(round(x, 0) != x){stop("input has to be an integer")}
  if(round(dig, 0) != dig){stop("number of digits has to be an integer")}
  
  n.zer <- dig - floor(log(x, 10)) - 1
  if(n.zer < 0){
    n.zer <- 0
    warning("Input is longer than defined number of digits!")
  }
  
  return(paste(paste(rep("0", n.zer), sep="", collapse=""), x, sep=""))
}

#--------------------------------------

zero_before2 <- function(x, dig){
  ## function to convert numbers to strings and add zeros to enable alphabetical order.
  ##
  ## Inputs: x   - integer
  ##         dig - number of digits
  ##
  ## Outputs: vector as character with zeros before 'shorter' numerals
  
  #if(round(x, 0) != x){stop("input has to be an integer")}
  if(round(dig, 0) != dig){stop("number of digits has to be an integer")}
  
  # treat zero and negative values
  sZero <- which(x == 0)
  sNeg <-  which(x < 0)
  x[sZero] <- 1
  x[sNeg] <- abs(x[sNeg])
  
  # use log10 to estimate number of digits
  nDig <- floor(log(x, 10))
  nDig[nDig < 0] <- 0
  nDig <- nDig + 1
  
  xOut <- rep('', length(x))
  for (i in 1 : dig) {
    zeros <- paste0(rep('0', dig - i), collapse = '')
    xOut[nDig == i] <- paste0(zeros, x[nDig == i])
  }
  xOut[nDig > dig] <- x
  # treat zero and negative values
  xOut[sZero] <- paste0(rep("0", dig), collapse="")
  xOut[sNeg] <- paste0('-', xOut[sNeg])
  
  return(xOut)
}

#--------------------------------------

## ========
## Helpers:
## ========

insert_missing_records <- function(zoo_ser, step){
  ## inserts time steps and assign to them NA value to make the input time series regular
  ##
  ## Arguments:
  ## dat  - data.frame of n columns with time stamp in the first column
  ## step - time step which should have the resulting time series (in minutes)
  ##
  ## returns:
  ## dat2 - data.frame with consistent time series
  
  require(zoo)
  
  dat <- data.frame('time' = index(zoo_ser), 'data'= coredata((zoo_ser)))
  
  dat2 <- reg_series(dat, step)
  zoo2_ser <- zoo(dat2[ , -1], dat2[ , 1])
  colnames(zoo2_ser) <- colnames(zoo_ser)
  
  return(zoo2_ser)
}

#--------------------------------------

reg_series <- function(dat, step){
  ## function to insert time steps to make time series regular
  ##
  ## Arguments:
  ## dat  - data.frame of n columns with time stamp in the first column
  ## step - time step which should have the resulting time series (in minutes)
  ##
  ## returns:
  ## dat2 - data.frame with consistent time series
  tim0 <- dat[,1]
  tim <- seq(tim0[1],tim0[length(tim0)], by=step*60)
  id.add <- which(is.element(tim, tim0)==T)
  dat2 <- as.data.frame(matrix(NA, ncol=ncol(dat), nrow=length(tim)))
  dat2[ ,1] <- tim
  dat2[id.add, -1] <- dat[ ,-1]
  
  # name colnames according to original data.frame
  colnames(dat2) <- colnames(dat)
  
  return(dat2)
}

#--------------------------------------

fun_refra <- function(fr){
    
    # Function to get refractive index
    # Arguments:   fr - frequency in GHz
    # Returns:  refr.out - refractive index, vector with real (first element)
    #                     and imaginary (second element) part of a refractive
    #                     index
    
    f <- fr*10^9
    
    c <- 2.99*10^8
    wl <- (c/f)*10^6   #wavelength in micrometers
    
    refr <- c(
        1.000E-02,0.968416,1.745E-03,
        1.099E-02,0.964778,2.370E-03,
        1.199E-02,0.960953,3.146E-03,
        1.300E-02,0.956954,4.072E-03,
        1.400E-02,0.952792,5.174E-03,
        1.600E-02,0.944124,7.958E-03,
        1.799E-02,0.934744,1.164E-02,
        2.000E-02,0.924583,1.636E-02,
        2.198E-02,0.913973,2.227E-02,
        2.399E-02,0.902694,2.950E-02,
        2.600E-02,0.890837,3.818E-02,
        2.799E-02,0.878766,4.851E-02,
        2.999E-02,0.866493,6.064E-02,
        3.199E-02,0.854141,7.461E-02,
        3.396E-02,0.842171,9.074E-02,
        3.597E-02,0.830720,1.093E-01,
        3.802E-02,0.819753,1.303E-01,
        3.999E-02,0.809997,1.534E-01,
        4.198E-02,0.802291,1.798E-01,
        4.395E-02,0.797737,2.088E-01,
        4.603E-02,0.797007,2.414E-01,
        4.797E-02,0.805579,2.766E-01,
        5.000E-02,0.820742,2.998E-01,
        5.200E-02,0.830957,3.153E-01,
        5.395E-02,0.835240,3.310E-01,
        5.598E-02,0.835295,3.498E-01,
        5.794E-02,0.831628,3.739E-01,
        5.998E-02,0.830901,4.119E-01,
        6.194E-02,0.840575,4.558E-01,
        6.397E-02,0.866994,5.033E-01,
        6.607E-02,0.903527,5.355E-01,
        6.808E-02,0.941801,5.634E-01,
        6.998E-02,0.981692,5.791E-01,
        7.194E-02,1.020921,5.859E-01,
        7.396E-02,1.049744,5.805E-01,
        7.603E-02,1.068724,5.859E-01,
        7.798E-02,1.087685,5.981E-01,
        7.998E-02,1.111682,6.135E-01,
        8.204E-02,1.140628,6.292E-01,
        8.395E-02,1.173382,6.453E-01,
        8.590E-02,1.214969,6.573E-01,
        8.790E-02,1.259495,6.573E-01,
        8.995E-02,1.302663,6.528E-01,
        9.204E-02,1.346760,6.439E-01,
        9.397E-02,1.387639,6.292E-01,
        9.594E-02,1.425425,6.050E-01,
        9.795E-02,1.455868,5.752E-01,
        1.000E-01,1.476628,5.430E-01,
        1.021E-01,1.493473,5.185E-01,
        1.040E-01,1.506677,4.929E-01,
        1.059E-01,1.516305,4.707E-01,
        1.079E-01,1.523589,4.485E-01,
        1.099E-01,1.528933,4.303E-01,
        1.119E-01,1.535363,4.148E-01,
        1.140E-01,1.543211,3.988E-01,
        1.159E-01,1.548070,3.826E-01,
        1.180E-01,1.553435,3.705E-01,
        1.199E-01,1.560870,3.596E-01,
        1.219E-01,1.570304,3.490E-01,
        1.239E-01,1.584638,3.387E-01,
        1.259E-01,1.606068,3.220E-01,
        1.279E-01,1.626822,2.876E-01,
        1.300E-01,1.633849,2.392E-01,
        1.321E-01,1.619420,1.870E-01,
        1.340E-01,1.586268,1.489E-01,
        1.361E-01,1.536403,1.333E-01,
        1.380E-01,1.496271,1.422E-01,
        1.400E-01,1.471129,1.678E-01,
        1.419E-01,1.461485,1.927E-01,
        1.439E-01,1.460977,2.167E-01,
        1.459E-01,1.469275,2.409E-01,
        1.479E-01,1.489551,2.641E-01,
        1.500E-01,1.521276,2.772E-01,
        1.521E-01,1.559942,2.772E-01,
        1.542E-01,1.596861,2.581E-01,
        1.560E-01,1.620422,2.338E-01,
        1.581E-01,1.641473,2.022E-01,
        1.600E-01,1.650184,1.670E-01,
        1.622E-01,1.652917,1.351E-01,
        1.641E-01,1.653100,1.039E-01,
        1.660E-01,1.647245,7.241E-02,
        1.679E-01,1.635062,3.998E-02,
        1.698E-01,1.605555,3.998E-03,
        1.722E-01,1.568183,2.004E-03,
        1.742E-01,1.549412,1.182E-03,
        1.750E-01,1.543062,8.391E-04,
        1.799E-01,1.513343,5.995E-05,
        1.849E-01,1.491881,1.250E-06,
        1.901E-01,1.475183,3.622E-07,
        1.950E-01,1.462543,1.850E-07,
        2.000E-01,1.451724,1.101E-07,
        2.051E-01,1.442296,6.711E-08,
        2.099E-01,1.434685,3.844E-08,
        2.148E-01,1.427828,1.999E-08,
        2.198E-01,1.421603,1.270E-08,
        2.249E-01,1.415921,1.158E-08,
        2.301E-01,1.410702,1.101E-08,
        2.350E-01,1.406358,1.071E-08,
        2.399E-01,1.402321,1.049E-08,
        2.449E-01,1.398535,9.904E-09,
        2.500E-01,1.394993,9.307E-09,
        2.553E-01,1.391674,8.606E-09,
        2.600E-01,1.388881,7.994E-09,
        2.649E-01,1.386239,7.444E-09,
        2.698E-01,1.383726,6.852E-09,
        2.748E-01,1.381341,6.292E-09,
        2.799E-01,1.379072,5.791E-09,
        2.851E-01,1.376902,5.405E-09,
        2.897E-01,1.375086,4.795E-09,
        2.951E-01,1.373098,4.403E-09,
        2.999E-01,1.371437,4.148E-09,
        3.048E-01,1.369839,3.826E-09,
        3.097E-01,1.368287,3.546E-09,
        3.148E-01,1.366812,3.325E-09,
        3.199E-01,1.365376,3.190E-09,
        3.251E-01,1.363990,3.082E-09,
        3.304E-01,1.362616,2.984E-09,
        3.350E-01,1.361513,2.883E-09,
        3.396E-01,1.360441,2.766E-09,
        3.451E-01,1.359231,2.653E-09,
        3.499E-01,1.358224,2.528E-09,
        3.548E-01,1.357247,2.420E-09,
        3.597E-01,1.356295,2.316E-09,
        3.648E-01,1.355370,2.217E-09,
        3.698E-01,1.354470,2.117E-09,
        3.750E-01,1.353594,2.031E-09,
        3.802E-01,1.352740,1.940E-09,
        3.846E-01,1.352046,1.840E-09,
        3.899E-01,1.351231,1.761E-09,
        3.954E-01,1.350438,1.663E-09,
        3.999E-01,1.349793,1.580E-09,
        4.046E-01,1.349159,1.489E-09,
        4.102E-01,1.348417,1.422E-09,
        4.150E-01,1.347811,1.339E-09,
        4.198E-01,1.347219,1.258E-09,
        4.246E-01,1.346636,1.169E-09,
        4.295E-01,1.346066,1.088E-09,
        4.345E-01,1.345505,1.018E-09,
        4.395E-01,1.344956,9.393E-10,
        4.446E-01,1.344418,8.685E-10,
        4.498E-01,1.343889,8.087E-10,
        4.550E-01,1.343368,7.795E-10,
        4.603E-01,1.342858,7.600E-10,
        4.645E-01,1.342455,7.495E-10,
        4.699E-01,1.341961,7.291E-10,
        4.753E-01,1.341475,7.011E-10,
        4.797E-01,1.341093,7.092E-10,
        4.853E-01,1.340620,7.158E-10,
        4.898E-01,1.340248,7.342E-10,
        4.955E-01,1.339791,7.849E-10,
        5.000E-01,1.339430,9.243E-10,
        5.047E-01,1.339073,1.078E-09,
        5.105E-01,1.338635,1.267E-09,
        5.152E-01,1.338288,1.461E-09,
        5.200E-01,1.337944,1.570E-09,
        5.248E-01,1.337607,1.640E-09,
        5.297E-01,1.337273,1.757E-09,
        5.346E-01,1.336943,1.887E-09,
        5.395E-01,1.336615,2.098E-09,
        5.445E-01,1.336292,2.269E-09,
        5.495E-01,1.335972,2.442E-09,
        5.546E-01,1.335656,2.659E-09,
        5.598E-01,1.335344,2.869E-09,
        5.649E-01,1.335035,3.132E-09,
        5.702E-01,1.334729,3.434E-09,
        5.754E-01,1.334425,3.844E-09,
        5.794E-01,1.334200,4.434E-09,
        5.848E-01,1.333902,5.221E-09,
        5.902E-01,1.333609,6.365E-09,
        5.957E-01,1.333316,7.723E-09,
        5.998E-01,1.333100,9.634E-09,
        6.053E-01,1.332813,1.132E-08,
        6.095E-01,1.332598,1.238E-08,
        6.152E-01,1.332317,1.330E-08,
        6.194E-01,1.332106,1.399E-08,
        6.252E-01,1.331826,1.472E-08,
        6.295E-01,1.331619,1.502E-08,
        6.353E-01,1.331345,1.552E-08,
        6.397E-01,1.331144,1.570E-08,
        6.457E-01,1.330877,1.606E-08,
        6.501E-01,1.330683,1.674E-08,
        6.546E-01,1.330490,1.777E-08,
        6.607E-01,1.330238,1.940E-08,
        6.653E-01,1.330052,2.031E-08,
        6.699E-01,1.329869,2.098E-08,
        6.745E-01,1.329690,2.177E-08,
        6.808E-01,1.329452,2.300E-08,
        6.855E-01,1.329278,2.471E-08,
        6.902E-01,1.329106,2.653E-08,
        6.950E-01,1.328938,2.963E-08,
        6.998E-01,1.328769,3.348E-08,
        7.047E-01,1.328603,4.100E-08,
        7.096E-01,1.328440,4.998E-08,
        7.145E-01,1.328279,5.995E-08,
        7.194E-01,1.328120,7.291E-08,
        7.244E-01,1.327963,9.137E-08,
        7.295E-01,1.327808,1.150E-07,
        7.345E-01,1.327652,1.348E-07,
        7.396E-01,1.327502,1.458E-07,
        7.447E-01,1.327350,1.530E-07,
        7.499E-01,1.327201,1.559E-07,
        7.551E-01,1.327055,1.580E-07,
        7.603E-01,1.326909,1.580E-07,
        7.656E-01,1.326764,1.570E-07,
        7.691E-01,1.326667,1.527E-07,
        7.745E-01,1.326524,1.478E-07,
        7.798E-01,1.326382,1.409E-07,
        7.852E-01,1.326244,1.339E-07,
        7.907E-01,1.326104,1.282E-07,
        7.943E-01,1.326012,1.258E-07,
        7.998E-01,1.325874,1.250E-07,
        8.054E-01,1.325739,1.270E-07,
        8.091E-01,1.325648,1.330E-07,
        8.147E-01,1.325512,1.448E-07,
        8.204E-01,1.325379,1.621E-07,
        8.241E-01,1.325290,1.819E-07,
        8.299E-01,1.325157,2.041E-07,
        8.356E-01,1.325025,2.243E-07,
        8.395E-01,1.324937,2.459E-07,
        8.453E-01,1.324805,2.690E-07,
        8.492E-01,1.324718,2.929E-07,
        8.551E-01,1.324590,3.153E-07,
        8.590E-01,1.324502,3.348E-07,
        8.650E-01,1.324373,3.546E-07,
        8.710E-01,1.324244,3.748E-07,
        8.750E-01,1.324159,3.907E-07,
        8.790E-01,1.324074,4.053E-07,
        8.851E-01,1.323946,4.234E-07,
        8.892E-01,1.323859,4.403E-07,
        8.954E-01,1.323732,4.622E-07,
        8.995E-01,1.323648,4.862E-07,
        9.057E-01,1.323520,5.150E-07,
        9.099E-01,1.323434,5.699E-07,
        9.141E-01,1.323351,6.696E-07,
        9.204E-01,1.323222,8.304E-07,
        9.247E-01,1.323138,1.060E-06,
        9.290E-01,1.323054,1.368E-06,
        9.354E-01,1.322926,1.771E-06,
        9.397E-01,1.322842,2.169E-06,
        9.441E-01,1.322757,2.557E-06,
        9.506E-01,1.322630,2.932E-06,
        9.550E-01,1.322546,3.190E-06,
        9.594E-01,1.322462,3.358E-06,
        9.661E-01,1.322333,3.464E-06,
        9.705E-01,1.322249,3.502E-06,
        9.750E-01,1.322165,3.480E-06,
        9.795E-01,1.322080,3.418E-06,
        9.840E-01,1.321994,3.336E-06,
        9.908E-01,1.321866,3.253E-06,
        9.954E-01,1.321780,3.131E-06,
        1.000E+00,1.321695,3.000E-06,
        1.009E+00,1.321521,2.688E-06,
        1.021E+00,1.321303,2.352E-06,
        1.030E+00,1.321128,2.001E-06,
        1.040E+00,1.320952,1.690E-06,
        1.050E+00,1.320775,1.419E-06,
        1.059E+00,1.320596,1.299E-06,
        1.069E+00,1.320416,1.259E-06,
        1.079E+00,1.320233,1.329E-06,
        1.089E+00,1.320051,1.499E-06,
        1.099E+00,1.319865,1.708E-06,
        1.109E+00,1.319678,2.038E-06,
        1.119E+00,1.319488,2.628E-06,
        1.130E+00,1.319296,3.869E-06,
        1.140E+00,1.319103,5.951E-06,
        1.151E+00,1.318909,9.306E-06,
        1.159E+00,1.318763,1.069E-05,
        1.169E+00,1.318566,1.120E-05,
        1.180E+00,1.318366,1.160E-05,
        1.191E+00,1.318162,1.181E-05,
        1.199E+00,1.318008,1.199E-05,
        1.211E+00,1.317799,1.191E-05,
        1.219E+00,1.317641,1.179E-05,
        1.230E+00,1.317427,1.160E-05,
        1.239E+00,1.317263,1.139E-05,
        1.250E+00,1.317042,1.100E-05,
        1.259E+00,1.316873,1.079E-05,
        1.271E+00,1.316645,1.090E-05,
        1.279E+00,1.316470,1.139E-05,
        1.291E+00,1.316233,1.221E-05,
        1.300E+00,1.316052,1.400E-05,
        1.309E+00,1.315868,1.639E-05,
        1.321E+00,1.315618,1.912E-05,
        1.330E+00,1.315425,2.251E-05,
        1.340E+00,1.315228,2.849E-05,
        1.349E+00,1.315031,4.047E-05,
        1.361E+00,1.314760,4.505E-05,
        1.371E+00,1.314547,5.804E-05,
        1.380E+00,1.314329,7.802E-05,
        1.390E+00,1.314104,1.060E-04,
        1.400E+00,1.313871,1.530E-04,
        1.409E+00,1.313671,2.540E-04,
        1.419E+00,1.313518,3.197E-04,
        1.429E+00,1.313373,3.538E-04,
        1.439E+00,1.313220,3.629E-04,
        1.449E+00,1.313055,3.637E-04,
        1.459E+00,1.312888,3.604E-04,
        1.469E+00,1.312715,3.387E-04,
        1.479E+00,1.312525,3.018E-04,
        1.489E+00,1.312318,2.659E-04,
        1.500E+00,1.312093,2.248E-04,
        1.510E+00,1.311852,1.958E-04,
        1.521E+00,1.311604,1.741E-04,
        1.531E+00,1.311352,1.602E-04,
        1.542E+00,1.311097,1.441E-04,
        1.549E+00,1.310923,1.348E-04,
        1.560E+00,1.310659,1.240E-04,
        1.570E+00,1.310387,1.140E-04,
        1.581E+00,1.310114,1.071E-04,
        1.589E+00,1.309928,9.940E-05,
        1.600E+00,1.309642,9.347E-05,
        1.611E+00,1.309352,8.804E-05,
        1.622E+00,1.309055,8.310E-05,
        1.629E+00,1.308855,8.096E-05,
        1.641E+00,1.308548,7.903E-05,
        1.648E+00,1.308341,7.591E-05,
        1.660E+00,1.308021,7.398E-05,
        1.671E+00,1.307672,7.404E-05,
        1.679E+00,1.307435,7.495E-05,
        1.690E+00,1.307073,7.601E-05,
        1.698E+00,1.306829,7.743E-05,
        1.710E+00,1.306453,8.050E-05,
        1.722E+00,1.306070,8.410E-05,
        1.730E+00,1.305809,8.900E-05,
        1.742E+00,1.305413,9.510E-05,
        1.750E+00,1.305142,1.000E-04,
        1.762E+00,1.304727,1.051E-04,
        1.770E+00,1.304442,1.120E-04,
        1.778E+00,1.304155,1.219E-04,
        1.791E+00,1.303718,1.330E-04,
        1.799E+00,1.303418,1.359E-04,
        1.811E+00,1.302947,1.371E-04,
        1.820E+00,1.302616,1.380E-04,
        1.828E+00,1.302269,1.418E-04,
        1.841E+00,1.301709,1.552E-04,
        1.849E+00,1.301291,1.861E-04,
        1.862E+00,1.300633,3.205E-04,
        1.871E+00,1.300214,5.209E-04,
        1.879E+00,1.299860,7.224E-04,
        1.888E+00,1.299545,9.221E-04,
        1.901E+00,1.298998,1.161E-03,
        1.910E+00,1.298791,1.678E-03,
        1.919E+00,1.298793,1.827E-03,
        1.932E+00,1.298681,1.922E-03,
        1.941E+00,1.298590,1.909E-03,
        1.950E+00,1.298472,1.848E-03,
        1.959E+00,1.298308,1.717E-03,
        1.968E+00,1.298051,1.548E-03,
        1.982E+00,1.297607,1.402E-03,
        1.991E+00,1.297292,1.250E-03,
        2.000E+00,1.296913,1.101E-03,
        2.009E+00,1.296499,9.904E-04,
        2.018E+00,1.296066,8.888E-04,
        2.028E+00,1.295606,8.050E-04,
        2.042E+00,1.294919,7.392E-04,
        2.051E+00,1.294457,6.742E-04,
        2.061E+00,1.293973,6.206E-04,
        2.070E+00,1.293476,5.725E-04,
        2.080E+00,1.292966,5.294E-04,
        2.089E+00,1.292438,4.884E-04,
        2.099E+00,1.291899,4.643E-04,
        2.109E+00,1.291353,4.403E-04,
        2.118E+00,1.290795,4.176E-04,
        2.128E+00,1.290221,3.970E-04,
        2.138E+00,1.289634,3.826E-04,
        2.148E+00,1.289033,3.705E-04,
        2.158E+00,1.288418,3.587E-04,
        2.168E+00,1.287787,3.506E-04,
        2.178E+00,1.287139,3.434E-04,
        2.188E+00,1.286474,3.395E-04,
        2.198E+00,1.285790,3.379E-04,
        2.208E+00,1.285087,3.387E-04,
        2.218E+00,1.284365,3.410E-04,
        2.228E+00,1.283619,3.458E-04,
        2.239E+00,1.282852,3.571E-04,
        2.249E+00,1.282064,3.739E-04,
        2.259E+00,1.281256,3.898E-04,
        2.270E+00,1.280421,4.081E-04,
        2.280E+00,1.279561,4.293E-04,
        2.291E+00,1.278675,4.506E-04,
        2.301E+00,1.277755,4.686E-04,
        2.312E+00,1.276797,4.918E-04,
        2.317E+00,1.276305,5.114E-04,
        2.328E+00,1.275295,5.430E-04,
        2.339E+00,1.274257,5.995E-04,
        2.350E+00,1.273184,6.365E-04,
        2.360E+00,1.272062,6.852E-04,
        2.371E+00,1.270902,7.427E-04,
        2.382E+00,1.269682,7.921E-04,
        2.388E+00,1.269059,8.488E-04,
        2.399E+00,1.267787,9.095E-04,
        2.410E+00,1.266450,9.904E-04,
        2.421E+00,1.265059,1.071E-03,
        2.432E+00,1.263577,1.150E-03,
        2.438E+00,1.262824,1.250E-03,
        2.449E+00,1.261297,1.348E-03,
        2.460E+00,1.259683,1.472E-03,
        2.472E+00,1.257969,1.580E-03,
        2.477E+00,1.257102,1.709E-03,
        2.489E+00,1.255347,1.810E-03,
        2.500E+00,1.253465,1.900E-03,
        2.512E+00,1.251445,1.953E-03,
        2.518E+00,1.250383,1.990E-03,
        2.529E+00,1.248127,2.017E-03,
        2.541E+00,1.245672,2.069E-03,
        2.553E+00,1.243014,2.142E-03,
        2.564E+00,1.240167,2.269E-03,
        2.570E+00,1.238670,2.311E-03,
        2.576E+00,1.237089,2.338E-03,
        2.582E+00,1.235426,2.387E-03,
        2.588E+00,1.233679,2.425E-03,
        2.594E+00,1.231834,2.476E-03,
        2.606E+00,1.227769,2.575E-03,
        2.612E+00,1.225483,2.703E-03,
        2.618E+00,1.223082,2.977E-03,
        2.624E+00,1.220535,3.302E-03,
        2.630E+00,1.218078,4.016E-03,
        2.636E+00,1.215699,4.363E-03,
        2.649E+00,1.209954,4.828E-03,
        2.655E+00,1.206519,5.368E-03,
        2.661E+00,1.202951,6.278E-03,
        2.667E+00,1.199289,7.325E-03,
        2.673E+00,1.195340,8.547E-03,
        2.679E+00,1.191390,1.049E-02,
        2.685E+00,1.188087,1.270E-02,
        2.698E+00,1.179962,1.451E-02,
        2.704E+00,1.174582,1.640E-02,
        2.710E+00,1.168874,1.861E-02,
        2.716E+00,1.160993,2.050E-02,
        2.723E+00,1.152876,2.817E-02,
        2.729E+00,1.149520,3.800E-02,
        2.742E+00,1.142068,4.622E-02,
        2.748E+00,1.136183,5.480E-02,
        2.754E+00,1.132860,6.483E-02,
        2.761E+00,1.131711,7.444E-02,
        2.767E+00,1.132778,8.352E-02,
        2.780E+00,1.130913,9.285E-02,
        2.786E+00,1.127959,1.020E-01,
        2.793E+00,1.127558,1.119E-01,
        2.799E+00,1.129478,1.210E-01,
        2.812E+00,1.128413,1.312E-01,
        2.818E+00,1.125532,1.422E-01,
        2.825E+00,1.125351,1.541E-01,
        2.831E+00,1.127523,1.670E-01,
        2.838E+00,1.133346,1.798E-01,
        2.851E+00,1.142386,1.940E-01,
        2.858E+00,1.145545,2.060E-01,
        2.864E+00,1.152284,2.182E-01,
        2.871E+00,1.162372,2.290E-01,
        2.884E+00,1.178446,2.392E-01,
        2.891E+00,1.185419,2.493E-01,
        2.897E+00,1.195889,2.581E-01,
        2.904E+00,1.208002,2.647E-01,
        2.917E+00,1.229654,2.715E-01,
        2.924E+00,1.240033,2.759E-01,
        2.931E+00,1.252073,2.798E-01,
        2.938E+00,1.263935,2.804E-01,
        2.951E+00,1.285942,2.824E-01,
        2.958E+00,1.297762,2.817E-01,
        2.965E+00,1.307891,2.785E-01,
        2.979E+00,1.326310,2.759E-01,
        2.985E+00,1.334533,2.721E-01,
        2.999E+00,1.352917,2.721E-01,
        3.048E+00,1.411876,2.398E-01,
        3.097E+00,1.452013,1.918E-01,
        3.148E+00,1.466753,1.348E-01,
        3.199E+00,1.461522,9.243E-02,
        3.251E+00,1.449409,6.106E-02,
        3.304E+00,1.432585,3.688E-02,
        3.350E+00,1.417064,2.611E-02,
        3.396E+00,1.404875,1.949E-02,
        3.451E+00,1.393260,1.321E-02,
        3.499E+00,1.384213,9.393E-03,
        3.548E+00,1.376092,6.789E-03,
        3.597E+00,1.368863,5.150E-03,
        3.648E+00,1.362546,4.234E-03,
        3.698E+00,1.356937,3.596E-03,
        3.750E+00,1.351891,3.402E-03,
        3.802E+00,1.347393,3.402E-03,
        3.846E+00,1.343958,3.530E-03,
        3.899E+00,1.340174,3.800E-03,
        3.954E+00,1.336658,4.157E-03,
        3.999E+00,1.333929,4.600E-03,
        4.046E+00,1.331403,5.067E-03,
        4.102E+00,1.328504,5.621E-03,
        4.150E+00,1.326183,6.220E-03,
        4.198E+00,1.323997,6.883E-03,
        4.246E+00,1.321906,7.600E-03,
        4.295E+00,1.319948,8.449E-03,
        4.345E+00,1.318113,9.307E-03,
        4.395E+00,1.316398,1.030E-02,
        4.446E+00,1.314920,1.140E-02,
        4.498E+00,1.313587,1.238E-02,
        4.550E+00,1.312483,1.361E-02,
        4.603E+00,1.311785,1.472E-02,
        4.645E+00,1.311588,1.548E-02,
        4.699E+00,1.311451,1.570E-02,
        4.753E+00,1.311148,1.552E-02,
        4.797E+00,1.310657,1.499E-02,
        4.853E+00,1.309721,1.441E-02,
        4.898E+00,1.308720,1.370E-02,
        4.955E+00,1.307228,1.312E-02,
        5.000E+00,1.305885,1.241E-02,
        5.047E+00,1.304258,1.180E-02,
        5.105E+00,1.301965,1.111E-02,
        5.152E+00,1.299910,1.061E-02,
        5.200E+00,1.297550,1.011E-02,
        5.248E+00,1.294933,9.904E-03,
        5.297E+00,1.292117,9.790E-03,
        5.346E+00,1.289015,9.881E-03,
        5.395E+00,1.285729,1.030E-02,
        5.445E+00,1.282194,1.078E-02,
        5.495E+00,1.278291,1.158E-02,
        5.546E+00,1.273883,1.258E-02,
        5.598E+00,1.268802,1.418E-02,
        5.649E+00,1.262994,1.659E-02,
        5.702E+00,1.256584,2.031E-02,
        5.754E+00,1.248370,2.482E-02,
        5.794E+00,1.242239,3.295E-02,
        5.848E+00,1.234896,4.323E-02,
        5.902E+00,1.229289,6.220E-02,
        5.957E+00,1.231892,8.646E-02,
        5.998E+00,1.242862,1.069E-01,
        6.053E+00,1.268459,1.250E-01,
        6.095E+00,1.295314,1.309E-01,
        6.152E+00,1.330121,1.172E-01,
        6.194E+00,1.341605,8.786E-02,
        6.252E+00,1.339863,6.947E-02,
        6.295E+00,1.335754,5.699E-02,
        6.353E+00,1.329242,4.952E-02,
        6.397E+00,1.325038,4.485E-02,
        6.457E+00,1.320468,4.176E-02,
        6.501E+00,1.317726,3.925E-02,
        6.546E+00,1.314837,3.731E-02,
        6.607E+00,1.311404,3.563E-02,
        6.653E+00,1.309021,3.450E-02,
        6.699E+00,1.306716,3.371E-02,
        6.745E+00,1.304521,3.310E-02,
        6.808E+00,1.301901,3.272E-02,
        6.855E+00,1.300125,3.242E-02,
        6.902E+00,1.298382,3.220E-02,
        6.950E+00,1.296751,3.212E-02,
        6.998E+00,1.295193,3.197E-02,
        7.047E+00,1.293609,3.190E-02,
        7.096E+00,1.292093,3.197E-02,
        7.145E+00,1.290696,3.205E-02,
        7.194E+00,1.289296,3.205E-02,
        7.244E+00,1.287944,3.220E-02,
        7.295E+00,1.286624,3.220E-02,
        7.345E+00,1.285242,3.227E-02,
        7.396E+00,1.283912,3.242E-02,
        7.447E+00,1.282606,3.249E-02,
        7.499E+00,1.281248,3.257E-02,
        7.551E+00,1.279895,3.272E-02,
        7.603E+00,1.278508,3.279E-02,
        7.656E+00,1.277123,3.302E-02,
        7.691E+00,1.276220,3.310E-02,
        7.745E+00,1.274794,3.325E-02,
        7.798E+00,1.273363,3.348E-02,
        7.852E+00,1.271952,3.371E-02,
        7.907E+00,1.270543,3.395E-02,
        7.943E+00,1.269613,3.410E-02,
        7.998E+00,1.268163,3.426E-02,
        8.054E+00,1.266657,3.450E-02,
        8.091E+00,1.265652,3.466E-02,
        8.147E+00,1.264125,3.490E-02,
        8.204E+00,1.262564,3.514E-02,
        8.241E+00,1.261488,3.530E-02,
        8.299E+00,1.259903,3.563E-02,
        8.356E+00,1.258240,3.579E-02,
        8.395E+00,1.257072,3.604E-02,
        8.453E+00,1.255384,3.637E-02,
        8.492E+00,1.254220,3.654E-02,
        8.551E+00,1.252405,3.688E-02,
        8.590E+00,1.251193,3.714E-02,
        8.650E+00,1.249353,3.748E-02,
        8.710E+00,1.247433,3.783E-02,
        8.750E+00,1.246095,3.809E-02,
        8.790E+00,1.244791,3.844E-02,
        8.851E+00,1.242789,3.880E-02,
        8.892E+00,1.241424,3.916E-02,
        8.954E+00,1.239322,3.952E-02,
        8.995E+00,1.237862,3.988E-02,
        9.057E+00,1.235657,4.035E-02,
        9.099E+00,1.234142,4.072E-02,
        9.141E+00,1.232659,4.110E-02,
        9.204E+00,1.230259,4.148E-02,
        9.247E+00,1.228589,4.196E-02,
        9.290E+00,1.226967,4.234E-02,
        9.354E+00,1.224439,4.293E-02,
        9.397E+00,1.222699,4.333E-02,
        9.441E+00,1.220909,4.373E-02,
        9.506E+00,1.218113,4.434E-02,
        9.550E+00,1.216115,4.475E-02,
        9.594E+00,1.214136,4.537E-02,
        9.661E+00,1.211068,4.600E-02,
        9.705E+00,1.208909,4.664E-02,
        9.750E+00,1.206729,4.718E-02,
        9.795E+00,1.204471,4.784E-02,
        9.840E+00,1.202228,4.851E-02,
        9.908E+00,1.198600,4.929E-02,
        9.954E+00,1.195932,4.998E-02,
        1.000E+01,1.193164,5.079E-02,
        1.005E+01,1.190334,5.174E-02,
        1.009E+01,1.187365,5.270E-02,
        1.014E+01,1.183900,5.380E-02,
        1.021E+01,1.180893,5.805E-02,
        1.026E+01,1.178360,5.634E-02,
        1.030E+01,1.174182,5.845E-02,
        1.035E+01,1.170827,5.995E-02,
        1.040E+01,1.167354,6.191E-02,
        1.045E+01,1.163960,6.394E-02,
        1.050E+01,1.160584,6.619E-02,
        1.054E+01,1.157248,6.852E-02,
        1.059E+01,1.153843,7.092E-02,
        1.064E+01,1.150368,7.359E-02,
        1.069E+01,1.146959,7.652E-02,
        1.074E+01,1.143601,7.958E-02,
        1.079E+01,1.140345,8.294E-02,
        1.084E+01,1.137372,8.646E-02,
        1.089E+01,1.134419,8.970E-02,
        1.094E+01,1.131445,9.328E-02,
        1.099E+01,1.128640,9.678E-02,
        1.104E+01,1.125466,9.995E-02,
        1.109E+01,1.122010,1.039E-01,
        1.114E+01,1.118841,1.083E-01,
        1.119E+01,1.116059,1.129E-01,
        1.125E+01,1.113289,1.172E-01,
        1.130E+01,1.110334,1.218E-01,
        1.135E+01,1.107674,1.270E-01,
        1.140E+01,1.105361,1.321E-01,
        1.146E+01,1.103057,1.370E-01,
        1.151E+01,1.100705,1.422E-01,
        1.156E+01,1.097503,1.472E-01,
        1.159E+01,1.096584,1.520E-01,
        1.164E+01,1.096068,1.570E-01,
        1.169E+01,1.094339,1.621E-01,
        1.175E+01,1.092339,1.678E-01,
        1.180E+01,1.090622,1.741E-01,
        1.186E+01,1.089062,1.802E-01,
        1.191E+01,1.086474,1.865E-01,
        1.194E+01,1.086163,1.927E-01,
        1.199E+01,1.087480,1.990E-01,
        1.205E+01,1.087926,2.055E-01,
        1.211E+01,1.087993,2.112E-01,
        1.216E+01,1.086723,2.177E-01,
        1.219E+01,1.087212,2.238E-01,
        1.225E+01,1.089721,2.295E-01,
        1.230E+01,1.090913,2.359E-01,
        1.236E+01,1.091270,2.420E-01,
        1.239E+01,1.092375,2.476E-01,
        1.245E+01,1.095643,2.528E-01,
        1.250E+01,1.098011,2.593E-01,
        1.256E+01,1.099603,2.641E-01,
        1.259E+01,1.100816,2.690E-01,
        1.265E+01,1.104624,2.740E-01,
        1.271E+01,1.107403,2.791E-01,
        1.276E+01,1.108999,2.837E-01,
        1.279E+01,1.110319,2.883E-01,
        1.285E+01,1.114243,2.929E-01,
        1.291E+01,1.116753,2.977E-01,
        1.294E+01,1.118262,3.012E-01,
        1.300E+01,1.122067,3.060E-01,
        1.306E+01,1.124841,3.103E-01,
        1.309E+01,1.126485,3.139E-01,
        1.315E+01,1.130583,3.183E-01,
        1.321E+01,1.133825,3.227E-01,
        1.324E+01,1.135773,3.257E-01,
        1.330E+01,1.139515,3.295E-01,
        1.334E+01,1.141428,3.325E-01,
        1.340E+01,1.145850,3.363E-01,
        1.346E+01,1.149628,3.402E-01,
        1.349E+01,1.151643,3.426E-01,
        1.355E+01,1.156338,3.466E-01,
        1.361E+01,1.160150,3.490E-01,
        1.365E+01,1.161869,3.514E-01,
        1.371E+01,1.165763,3.546E-01,
        1.374E+01,1.167947,3.571E-01,
        1.380E+01,1.172049,3.596E-01,
        1.384E+01,1.174089,3.621E-01,
        1.390E+01,1.178513,3.646E-01,
        1.396E+01,1.182458,3.680E-01,
        1.400E+01,1.184740,3.696E-01,
        1.406E+01,1.189086,3.722E-01,
        1.409E+01,1.191399,3.739E-01,
        1.416E+01,1.195603,3.757E-01,
        1.419E+01,1.197623,3.774E-01,
        1.426E+01,1.201594,3.791E-01,
        1.429E+01,1.203552,3.809E-01,
        1.435E+01,1.207465,3.826E-01,
        1.439E+01,1.209428,3.844E-01,
        1.445E+01,1.213645,3.862E-01,
        1.449E+01,1.215328,3.871E-01,
        1.455E+01,1.218762,3.898E-01,
        1.459E+01,1.220973,3.916E-01,
        1.466E+01,1.225566,3.934E-01,
        1.469E+01,1.227627,3.943E-01,
        1.476E+01,1.231631,3.961E-01,
        1.479E+01,1.233597,3.970E-01,
        1.486E+01,1.237500,3.988E-01,
        1.489E+01,1.239445,3.998E-01,
        1.496E+01,1.243348,4.016E-01,
        1.500E+01,1.245318,4.025E-01,
        1.507E+01,1.249380,4.044E-01,
        1.510E+01,1.251704,4.053E-01,
        1.514E+01,1.253631,4.053E-01,
        1.521E+01,1.256977,4.072E-01,
        1.524E+01,1.258880,4.081E-01,
        1.531E+01,1.263173,4.100E-01,
        1.535E+01,1.265082,4.100E-01,
        1.542E+01,1.268440,4.119E-01,
        1.545E+01,1.270391,4.128E-01,
        1.549E+01,1.272559,4.138E-01,
        1.556E+01,1.276473,4.148E-01,
        1.560E+01,1.278233,4.157E-01,
        1.567E+01,1.282639,4.176E-01,
        1.570E+01,1.284709,4.176E-01,
        1.574E+01,1.286576,4.186E-01,
        1.581E+01,1.290576,4.196E-01,
        1.585E+01,1.292723,4.205E-01,
        1.589E+01,1.294706,4.205E-01,
        1.596E+01,1.298872,4.225E-01,
        1.600E+01,1.301310,4.225E-01,
        1.603E+01,1.303145,4.225E-01,
        1.611E+01,1.306556,4.234E-01,
        1.614E+01,1.308540,4.244E-01,
        1.622E+01,1.313112,4.254E-01,
        1.626E+01,1.315327,4.254E-01,
        1.629E+01,1.317122,4.254E-01,
        1.637E+01,1.320901,4.264E-01,
        1.641E+01,1.322675,4.264E-01,
        1.644E+01,1.324631,4.274E-01,
        1.648E+01,1.326773,4.274E-01,
        1.656E+01,1.330870,4.283E-01,
        1.660E+01,1.333056,4.283E-01,
        1.663E+01,1.334869,4.283E-01,
        1.671E+01,1.338869,4.293E-01,
        1.675E+01,1.341074,4.293E-01,
        1.679E+01,1.342949,4.293E-01,
        1.687E+01,1.347481,4.303E-01,
        1.690E+01,1.349696,4.293E-01,
        1.694E+01,1.351233,4.293E-01,
        1.698E+01,1.352834,4.293E-01,
        1.706E+01,1.356772,4.303E-01,
        1.710E+01,1.359301,4.303E-01,
        1.714E+01,1.361083,4.293E-01,
        1.722E+01,1.364655,4.303E-01,
        1.726E+01,1.367219,4.303E-01,
        1.730E+01,1.369211,4.293E-01,
        1.734E+01,1.370751,4.293E-01,
        1.742E+01,1.374519,4.293E-01,
        1.746E+01,1.376758,4.293E-01,
        1.750E+01,1.378598,4.283E-01,
        1.754E+01,1.380029,4.283E-01,
        1.762E+01,1.383660,4.283E-01,
        1.766E+01,1.385875,4.283E-01,
        1.770E+01,1.387734,4.274E-01,
        1.774E+01,1.389417,4.274E-01,
        1.778E+01,1.390838,4.264E-01,
        1.786E+01,1.394313,4.274E-01,
        1.791E+01,1.396377,4.264E-01,
        1.795E+01,1.398169,4.264E-01,
        1.799E+01,1.399826,4.254E-01,
        1.803E+01,1.401123,4.254E-01,
        1.811E+01,1.404604,4.254E-01,
        1.816E+01,1.406786,4.254E-01,
        1.820E+01,1.408657,4.244E-01,
        1.824E+01,1.410419,4.244E-01,
        1.828E+01,1.412092,4.234E-01,
        1.837E+01,1.415276,4.234E-01,
        1.841E+01,1.417548,4.234E-01,
        1.845E+01,1.419809,4.225E-01,
        1.849E+01,1.421557,4.215E-01,
        1.854E+01,1.422820,4.205E-01,
        1.862E+01,1.426178,4.205E-01,
        1.866E+01,1.428308,4.196E-01,
        1.871E+01,1.429982,4.186E-01,
        1.875E+01,1.431240,4.176E-01,
        1.879E+01,1.432797,4.176E-01,
        1.884E+01,1.434643,4.167E-01,
        1.888E+01,1.435881,4.157E-01,
        1.897E+01,1.439563,4.157E-01,
        1.901E+01,1.441618,4.138E-01,
        1.905E+01,1.442846,4.128E-01,
        1.910E+01,1.444197,4.119E-01,
        1.914E+01,1.445486,4.110E-01,
        1.919E+01,1.446666,4.100E-01,
        1.923E+01,1.447502,4.091E-01,
        1.932E+01,1.450255,4.091E-01,
        1.936E+01,1.452188,4.081E-01,
        1.941E+01,1.453825,4.072E-01,
        1.945E+01,1.455604,4.062E-01,
        1.950E+01,1.456898,4.044E-01,
        1.954E+01,1.457713,4.035E-01,
        1.959E+01,1.458719,4.025E-01,
        1.963E+01,1.459690,4.016E-01,
        1.968E+01,1.460391,4.007E-01,
        1.977E+01,1.463349,4.007E-01,
        1.982E+01,1.465400,3.988E-01,
        1.986E+01,1.466543,3.970E-01,
        1.991E+01,1.467000,3.952E-01,
        1.995E+01,1.467249,3.943E-01,
        2.000E+01,1.467642,3.934E-01,
        2.099E+01,1.483693,3.818E-01,
        2.198E+01,1.499422,3.722E-01,
        2.301E+01,1.516402,3.629E-01,
        2.399E+01,1.529309,3.482E-01,
        2.500E+01,1.537967,3.356E-01,
        2.600E+01,1.544080,3.227E-01,
        2.698E+01,1.546670,3.103E-01,
        2.799E+01,1.546272,2.991E-01,
        2.897E+01,1.542658,2.889E-01,
        2.999E+01,1.535500,2.817E-01,
        3.097E+01,1.527225,2.791E-01,
        3.199E+01,1.519076,2.798E-01,
        3.304E+01,1.511879,2.830E-01,
        3.396E+01,1.505906,2.863E-01,
        3.499E+01,1.498932,2.916E-01,
        3.597E+01,1.492960,2.991E-01,
        3.698E+01,1.486740,3.068E-01,
        3.802E+01,1.481006,3.190E-01,
        3.899E+01,1.478232,3.317E-01,
        3.999E+01,1.476571,3.442E-01,
        4.102E+01,1.475642,3.587E-01,
        4.198E+01,1.477194,3.739E-01,
        4.295E+01,1.480747,3.880E-01,
        4.395E+01,1.485266,4.025E-01,
        4.498E+01,1.491543,4.176E-01,
        4.603E+01,1.499424,4.323E-01,
        4.699E+01,1.508821,4.465E-01,
        4.797E+01,1.520272,4.579E-01,
        4.898E+01,1.531473,4.675E-01,
        5.000E+01,1.542270,4.773E-01,
        5.200E+01,1.567492,4.975E-01,
        5.395E+01,1.594131,5.079E-01,
        5.598E+01,1.619157,5.162E-01,
        5.794E+01,1.643739,5.233E-01,
        5.998E+01,1.669053,5.258E-01,
        6.194E+01,1.690223,5.246E-01,
        6.397E+01,1.709762,5.246E-01,
        6.607E+01,1.729441,5.233E-01,
        6.808E+01,1.747333,5.209E-01,
        6.998E+01,1.762824,5.174E-01,
        7.194E+01,1.777162,5.138E-01,
        7.396E+01,1.790800,5.103E-01,
        7.603E+01,1.805539,5.079E-01,
        7.798E+01,1.819110,5.021E-01,
        7.998E+01,1.830882,4.964E-01,
        8.204E+01,1.842330,4.907E-01,
        8.395E+01,1.851943,4.839E-01,
        8.590E+01,1.859854,4.773E-01,
        8.790E+01,1.867327,4.718E-01,
        8.995E+01,1.874242,4.654E-01,
        9.204E+01,1.880545,4.600E-01,
        9.397E+01,1.886330,4.548E-01,
        9.594E+01,1.891384,4.485E-01,
        9.795E+01,1.895435,4.434E-01,
        1.000E+02,1.899131,4.383E-01,
        1.099E+02,1.907505,4.176E-01,
        1.199E+02,1.911671,4.167E-01,
        1.300E+02,1.919973,4.186E-01,
        1.400E+02,1.927412,4.205E-01,
        1.500E+02,1.934154,4.264E-01,
        1.600E+02,1.941655,4.333E-01,
        1.698E+02,1.948419,4.403E-01,
        1.799E+02,1.955736,4.506E-01,
        1.901E+02,1.965156,4.611E-01,
        2.000E+02,1.974559,4.697E-01,
        2.099E+02,1.983438,4.784E-01,
        2.198E+02,1.992287,4.873E-01,
        2.301E+02,2.001418,4.964E-01,
        2.399E+02,2.010446,5.056E-01,
        2.500E+02,2.020318,5.138E-01,
        2.600E+02,2.029224,5.209E-01,
        2.698E+02,2.037243,5.282E-01,
        2.799E+02,2.045135,5.355E-01,
        2.897E+02,2.052476,5.430E-01,
        2.999E+02,2.059773,5.505E-01,
        3.199E+02,2.073976,5.660E-01,
        3.396E+02,2.086956,5.805E-01,
        3.597E+02,2.099543,5.967E-01,
        3.802E+02,2.112811,6.135E-01,
        3.999E+02,2.125742,6.292E-01,
        4.198E+02,2.139507,6.453E-01,
        4.395E+02,2.153213,6.589E-01,
        4.603E+02,2.166254,6.727E-01,
        4.797E+02,2.177335,6.852E-01,
        5.000E+02,2.188736,6.995E-01,
        5.200E+02,2.200349,7.125E-01,
        5.395E+02,2.210869,7.241E-01,
        5.598E+02,2.220374,7.359E-01,
        5.794E+02,2.228339,7.478E-01,
        5.998E+02,2.236685,7.617E-01,
        6.501E+02,2.254575,7.921E-01,
        6.998E+02,2.270109,8.314E-01,
        7.499E+02,2.290196,8.725E-01,
        7.998E+02,2.312599,9.116E-01,
        8.492E+02,2.337241,9.501E-01,
        8.995E+02,2.363856,9.835E-01,
        9.506E+02,2.385313,1.011E+00,
        1.000E+03,2.399111,1.042E+00,
        1.100E+03,2.436760,1.120E+00,
        1.200E+03,2.481153,1.191E+00,
        1.300E+03,2.527536,1.261E+00,
        1.400E+03,2.577344,1.330E+00,
        1.500E+03,2.629097,1.393E+00,
        1.600E+03,2.679108,1.452E+00,
        1.700E+03,2.729264,1.511E+00,
        1.800E+03,2.781861,1.567E+00,
        1.900E+03,2.831974,1.617E+00,
        2.000E+03,2.881863,1.670E+00,
        2.100E+03,2.933900,1.718E+00,
        2.200E+03,2.983258,1.763E+00,
        2.300E+03,3.032401,1.809E+00,
        2.400E+03,3.084049,1.853E+00,
        2.500E+03,3.133464,1.891E+00,
        2.600E+03,3.179887,1.931E+00,
        2.700E+03,3.228984,1.973E+00,
        2.800E+03,3.279470,2.009E+00,
        2.900E+03,3.326631,2.043E+00,
        3.000E+03,3.374610,2.079E+00,
        3.100E+03,3.422465,2.110E+00,
        3.200E+03,3.468221,2.142E+00,
        3.300E+03,3.516889,2.174E+00,
        3.400E+03,3.563346,2.199E+00,
        3.500E+03,3.607096,2.228E+00,
        3.600E+03,3.650102,2.255E+00,
        3.700E+03,3.695213,2.286E+00,
        3.800E+03,3.741930,2.310E+00,
        3.900E+03,3.785136,2.333E+00,
        4.000E+03,3.829496,2.360E+00,
        4.100E+03,3.873564,2.380E+00,
        4.200E+03,3.917021,2.404E+00,
        4.300E+03,3.960586,2.423E+00,
        4.400E+03,4.003601,2.445E+00,
        4.500E+03,4.045111,2.460E+00,
        4.600E+03,4.084851,2.481E+00,
        4.700E+03,4.125763,2.500E+00,
        4.800E+03,4.166540,2.518E+00,
        4.900E+03,4.207585,2.535E+00,
        5.000E+03,4.248425,2.551E+00,
        5.100E+03,4.288766,2.567E+00,
        5.200E+03,4.328263,2.581E+00,
        5.300E+03,4.367201,2.595E+00,
        5.400E+03,4.403706,2.607E+00,
        5.500E+03,4.442443,2.625E+00,
        5.600E+03,4.482558,2.636E+00,
        5.700E+03,4.518750,2.647E+00,
        5.800E+03,4.555811,2.662E+00,
        5.901E+03,4.593558,2.671E+00,
        6.001E+03,4.631138,2.685E+00,
        6.100E+03,4.667698,2.692E+00,
        6.200E+03,4.704528,2.705E+00,
        6.299E+03,4.740422,2.711E+00,
        6.400E+03,4.776470,2.723E+00,
        6.500E+03,4.811400,2.727E+00,
        6.599E+03,4.844068,2.737E+00,
        6.700E+03,4.881015,2.747E+00,
        6.800E+03,4.915201,2.750E+00,
        6.899E+03,4.946351,2.758E+00,
        7.000E+03,4.979800,2.766E+00,
        7.101E+03,5.013994,2.774E+00,
        7.199E+03,5.047771,2.780E+00,
        7.300E+03,5.083439,2.787E+00,
        7.399E+03,5.116001,2.786E+00,
        7.501E+03,5.146161,2.792E+00,
        7.600E+03,5.177179,2.796E+00,
        7.700E+03,5.209531,2.801E+00,
        7.800E+03,5.241539,2.805E+00,
        7.900E+03,5.273172,2.808E+00,
        8.000E+03,5.304929,2.811E+00,
        8.100E+03,5.336323,2.814E+00,
        8.200E+03,5.367389,2.816E+00,
        8.300E+03,5.398286,2.824E+00,
        8.400E+03,5.428878,2.825E+00,
        8.500E+03,5.459208,2.826E+00,
        8.600E+03,5.489262,2.827E+00,
        8.700E+03,5.519027,2.827E+00,
        8.800E+03,5.548489,2.833E+00,
        8.900E+03,5.577699,2.833E+00,
        9.000E+03,5.606586,2.832E+00,
        9.100E+03,5.635201,2.830E+00,
        9.200E+03,5.663535,2.835E+00,
        9.300E+03,5.691521,2.833E+00,
        9.400E+03,5.719272,2.831E+00,
        9.500E+03,5.746661,2.835E+00,
        9.600E+03,5.773802,2.832E+00,
        9.700E+03,5.800631,2.835E+00,
        9.800E+03,5.827179,2.831E+00,
        9.900E+03,5.853423,2.828E+00,
        1.000E+04,5.879378,2.830E+00,
        1.100E+04,6.131865,2.807E+00,
        1.200E+04,6.346035,2.773E+00,
        1.300E+04,6.538143,2.733E+00,
        1.400E+04,6.711149,2.685E+00,
        1.500E+04,6.867192,2.630E+00,
        1.600E+04,7.007965,2.576E+00,
        1.700E+04,7.135674,2.519E+00,
        1.800E+04,7.252419,2.461E+00,
        1.900E+04,7.358822,2.397E+00,
        2.000E+04,7.455943,2.338E+00,
        2.100E+04,7.544423,2.280E+00,
        2.200E+04,7.625553,2.224E+00,
        2.300E+04,7.701126,2.170E+00,
        2.400E+04,7.768902,2.113E+00,
        2.500E+04,7.831158,2.064E+00,
        2.600E+04,7.889643,2.011E+00,
        2.700E+04,7.941322,1.961E+00,
        2.800E+04,7.989355,1.914E+00,
        2.900E+04,8.033791,1.868E+00,
        3.000E+04,8.074469,1.824E+00,
        3.100E+04,8.112180,1.781E+00,
        3.200E+04,8.147128,1.741E+00,
        3.300E+04,8.179843,1.701E+00,
        3.400E+04,8.209818,1.663E+00,
        3.500E+04,8.238281,1.627E+00,
        3.600E+04,8.264599,1.591E+00,
        3.700E+04,8.288448,1.557E+00,
        3.800E+04,8.311297,1.525E+00,
        3.900E+04,8.332788,1.493E+00,
        4.000E+04,8.352700,1.463E+00,
        4.100E+04,8.371979,1.434E+00,
        4.200E+04,8.389151,1.405E+00,
        4.300E+04,8.405156,1.378E+00,
        4.400E+04,8.420858,1.352E+00,
        4.500E+04,8.435660,1.327E+00,
        4.600E+04,8.449912,1.302E+00,
        4.700E+04,8.462521,1.278E+00,
        4.800E+04,8.474673,1.256E+00,
        4.900E+04,8.486470,1.233E+00,
        5.000E+04,8.497290,1.212E+00,
        5.100E+04,8.507814,1.191E+00,
        5.200E+04,8.517490,1.171E+00,
        5.300E+04,8.526995,1.152E+00,
        5.400E+04,8.535947,1.133E+00,
        5.500E+04,8.545560,1.115E+00,
        5.600E+04,8.554652,1.097E+00,
        5.700E+04,8.563363,1.080E+00,
        5.800E+04,8.571605,1.063E+00,
        5.901E+04,8.579613,1.047E+00,
        6.001E+04,8.587191,1.031E+00,
        6.100E+04,8.594362,1.016E+00,
        6.200E+04,8.601342,1.001E+00,
        6.299E+04,8.607948,9.867E-01,
        6.400E+04,8.614382,9.727E-01,
        6.500E+04,8.620473,9.589E-01,
        6.599E+04,8.626323,9.456E-01,
        6.700E+04,8.632030,9.328E-01,
        6.800E+04,8.637431,9.202E-01,
        6.899E+04,8.642623,9.078E-01,
        7.000E+04,8.647693,8.960E-01,
        7.101E+04,8.652568,8.845E-01,
        7.199E+04,8.657181,8.728E-01,
        7.300E+04,8.661691,8.618E-01,
        7.399E+04,8.666029,8.509E-01,
        7.501E+04,8.670273,8.406E-01,
        7.600E+04,8.674287,8.302E-01,
        7.700E+04,8.678217,8.201E-01,
        7.800E+04,8.681999,8.104E-01,
        7.900E+04,8.685638,8.007E-01,
        8.000E+04,8.689204,7.914E-01,
        8.100E+04,8.692637,7.822E-01,
        8.200E+04,8.695941,7.733E-01,
        8.300E+04,8.699181,7.645E-01,
        8.400E+04,8.702404,7.560E-01,
        8.500E+04,8.705507,7.476E-01,
        8.600E+04,8.708553,7.394E-01,
        8.700E+04,8.711486,7.314E-01,
        8.800E+04,8.714366,7.236E-01,
        8.900E+04,8.717139,7.159E-01,
        8.999E+04,8.719808,7.083E-01,
        9.099E+04,8.722432,7.009E-01,
        9.200E+04,8.725011,6.938E-01,
        9.300E+04,8.727494,6.867E-01,
        9.399E+04,8.729884,6.798E-01,
        9.499E+04,8.732236,6.730E-01,
        9.601E+04,8.734549,6.663E-01,
        9.701E+04,8.736776,6.598E-01,
        9.799E+04,8.738919,6.534E-01,
        9.899E+04,8.741029,6.471E-01,
        1.000E+05,8.743107,6.409E-01,
        1.099E+05,8.761392,5.845E-01,
        1.199E+05,8.776133,5.380E-01,
        1.300E+05,8.788091,4.983E-01,
        1.400E+05,8.797836,4.636E-01,
        1.500E+05,8.804337,4.337E-01,
        1.600E+05,8.809677,4.073E-01,
        1.698E+05,8.814082,3.837E-01,
        1.799E+05,8.817867,3.630E-01,
        1.901E+05,8.821129,3.448E-01,
        2.000E+05,8.823834,3.277E-01,
        2.099E+05,8.826184,3.123E-01,
        2.198E+05,8.828230,2.983E-01,
        2.301E+05,8.830103,2.860E-01,
        2.399E+05,8.831656,2.741E-01,
        2.500E+05,8.833091,2.634E-01,
        2.600E+05,8.834348,2.535E-01,
        2.698E+05,8.835449,2.441E-01,
        2.799E+05,8.836475,2.356E-01,
        2.897E+05,8.837375,2.274E-01,
        2.999E+05,8.838141,2.202E-01,
        3.097E+05,8.838809,2.132E-01,
        3.199E+05,8.839434,2.064E-01,
        3.304E+05,8.840020,2.008E-01,
        3.396E+05,8.840491,1.944E-01,
        3.499E+05,8.840974,1.891E-01,
        3.597E+05,8.841394,1.840E-01,
        3.698E+05,8.841790,1.790E-01,
        3.802E+05,8.842166,1.745E-01,
        3.899E+05,8.842491,1.701E-01,
        3.999E+05,8.842801,1.659E-01,
        4.102E+05,8.843095,1.621E-01,
        4.198E+05,8.843349,1.580E-01,
        4.295E+05,8.843591,1.545E-01,
        4.395E+05,8.843823,1.509E-01,
        4.498E+05,8.844044,1.478E-01,
        4.603E+05,8.844254,1.448E-01,
        4.699E+05,8.844436,1.415E-01,
        4.797E+05,8.844609,1.386E-01,
        4.898E+05,8.844776,1.358E-01,
        5.000E+05,8.844936,1.333E-01,
        5.105E+05,8.845089,1.309E-01,
        5.200E+05,8.845220,1.282E-01,
        5.297E+05,8.845345,1.258E-01,
        5.395E+05,8.845467,1.233E-01,
        5.495E+05,8.845583,1.213E-01,
        5.598E+05,8.845696,1.191E-01,
        5.702E+05,8.845804,1.172E-01,
        5.794E+05,8.845896,1.150E-01,
        5.902E+05,8.845997,1.132E-01,
        5.998E+05,8.846082,1.114E-01,
        6.095E+05,8.846164,1.093E-01,
        6.194E+05,8.846244,1.076E-01,
        6.295E+05,8.846321,1.061E-01,
        6.397E+05,8.846396,1.044E-01,
        6.501E+05,8.846468,1.030E-01,
        6.607E+05,8.846538,1.013E-01,
        6.699E+05,8.846596,9.995E-02,
        6.808E+05,8.846662,9.858E-02,
        6.902E+05,8.846716,9.700E-02,
        6.998E+05,8.846769,9.567E-02,
        7.096E+05,8.846821,9.436E-02,
        7.194E+05,8.846871,9.285E-02,
        7.295E+05,8.846920,9.179E-02,
        7.396E+05,8.846967,9.053E-02,
        7.499E+05,8.847013,8.929E-02,
        7.603E+05,8.847058,8.827E-02,
        7.691E+05,8.847094,8.705E-02,
        7.798E+05,8.847137,8.586E-02,
        7.907E+05,8.847178,8.488E-02,
        7.998E+05,8.847212,8.391E-02,
        8.091E+05,8.847244,8.275E-02,
        8.204E+05,8.847283,8.181E-02,
        8.299E+05,8.847314,8.087E-02,
        8.395E+05,8.847344,7.976E-02,
        8.492E+05,8.847374,7.885E-02,
        8.590E+05,8.847403,7.795E-02,
        8.710E+05,8.847437,7.723E-02,
        8.790E+05,8.847459,7.617E-02,
        8.892E+05,8.847485,7.547E-02,
        8.995E+05,8.847512,7.461E-02,
        9.099E+05,8.847537,7.375E-02,
        9.204E+05,8.847563,7.308E-02,
        9.290E+05,8.847582,7.224E-02,
        9.397E+05,8.847606,7.142E-02,
        9.506E+05,8.847630,7.076E-02,
        9.594E+05,8.847648,6.995E-02,
        9.705E+05,8.847671,6.931E-02,
        9.795E+05,8.847688,6.852E-02,
        9.908E+05,8.847710,6.804E-02,
        1.000E+06,8.847727,6.727E-02,
        1.099E+06,8.847889,6.121E-02,
        1.199E+06,8.848013,5.621E-02,
        1.300E+06,8.848109,5.197E-02,
        1.400E+06,8.848183,4.828E-02,
        1.500E+06,8.848243,4.506E-02,
        1.600E+06,8.848291,4.234E-02,
        1.698E+06,8.848330,3.979E-02,
        1.799E+06,8.848364,3.765E-02,
        1.901E+06,8.848392,3.571E-02,
        2.000E+06,8.848415,3.395E-02,
        2.099E+06,8.848435,3.234E-02,
        2.198E+06,8.848452,3.089E-02,
        2.301E+06,8.848467,2.957E-02,
        2.399E+06,8.848480,2.837E-02,
        2.500E+06,8.848491,2.721E-02,
        2.600E+06,8.848501,2.623E-02,
        2.698E+06,8.848510,2.522E-02,
        2.799E+06,8.848518,2.437E-02,
        2.897E+06,8.848524,2.348E-02,
        2.999E+06,8.848531,2.274E-02,
        3.097E+06,8.848536,2.202E-02,
        3.199E+06,8.848541,2.132E-02,
        3.304E+06,8.848546,2.074E-02,
        3.396E+06,8.848550,2.008E-02,
        3.499E+06,8.848554,1.953E-02,
        3.597E+06,8.848557,1.900E-02,
        3.698E+06,8.848560,1.848E-02,
        3.802E+06,8.848563,1.802E-02,
        3.899E+06,8.848566,1.757E-02,
        3.999E+06,8.848568,1.713E-02,
        4.102E+06,8.848570,1.670E-02,
        4.198E+06,8.848572,1.632E-02,
        4.295E+06,8.848574,1.591E-02,
        4.395E+06,8.848576,1.559E-02,
        4.498E+06,8.848577,1.523E-02,
        4.603E+06,8.848579,1.492E-02,
        4.699E+06,8.848580,1.458E-02,
        4.797E+06,8.848581,1.428E-02,
        4.898E+06,8.848582,1.399E-02,
        5.000E+06,8.848584,1.373E-02,
        5.105E+06,8.848585,1.348E-02,
        5.200E+06,8.848585,1.321E-02,
        5.297E+06,8.848586,1.297E-02,
        5.395E+06,8.848587,1.273E-02,
        5.495E+06,8.848588,1.250E-02,
        5.598E+06,8.848589,1.227E-02,
        5.702E+06,8.848589,1.207E-02,
        5.794E+06,8.848590,1.185E-02,
        5.902E+06,8.848591,1.166E-02,
        5.998E+06,8.848591,1.148E-02,
        6.095E+06,8.848592,1.129E-02,
        6.194E+06,8.848592,1.111E-02,
        6.295E+06,8.848593,1.091E-02,
        6.397E+06,8.848593,1.076E-02,
        6.501E+06,8.848594,1.061E-02,
        6.607E+06,8.848594,1.044E-02,
        6.699E+06,8.848594,1.030E-02,
        6.808E+06,8.848595,1.016E-02,
        6.902E+06,8.848595,9.995E-03,
        6.998E+06,8.848595,9.858E-03,
        7.096E+06,8.848596,9.700E-03,
        7.194E+06,8.848596,9.567E-03,
        7.295E+06,8.848596,9.458E-03,
        7.396E+06,8.848596,9.328E-03,
        7.499E+06,8.848597,9.200E-03,
        7.603E+06,8.848597,9.095E-03,
        7.691E+06,8.848597,8.949E-03,
        7.798E+06,8.848597,8.847E-03,
        7.907E+06,8.848598,8.746E-03,
        7.998E+06,8.848598,8.626E-03,
        8.091E+06,8.848598,8.527E-03,
        8.204E+06,8.848598,8.429E-03,
        8.299E+06,8.848598,8.333E-03,
        8.395E+06,8.848598,8.218E-03,
        8.492E+06,8.848598,8.124E-03,
        8.590E+06,8.848599,8.031E-03,
        8.710E+06,8.848599,7.958E-03,
        8.790E+06,8.848599,7.849E-03,
        8.892E+06,8.848599,7.759E-03,
        8.995E+06,8.848599,7.688E-03,
        9.099E+06,8.848599,7.600E-03,
        9.204E+06,8.848599,7.530E-03,
        9.290E+06,8.848599,7.444E-03,
        9.397E+06,8.848600,7.359E-03,
        9.506E+06,8.848600,7.291E-03,
        9.594E+06,8.848600,7.208E-03,
        9.705E+06,8.848600,7.142E-03,
        9.795E+06,8.848600,7.060E-03,
        9.908E+06,8.848600,7.011E-03,
        1.000E+07,8.848600,6.931E-03
    )
    
    
    refr <- matrix(refr, length(refr/3), 3, byrow=T)
    
    refr.out <- c(NA,NA)
    refr.out[1] <- approx(refr[,1], refr[,2], wl)$y
    refr.out[2] <- approx(refr[,1], refr[,3], wl)$y
    
    return(refr.out)
    
}
