#' Calculate weighted response rate over all strata in each wave. The weighted response rate
#' one of indicators to measure nonresponse error. Each wave has N simulated RR.
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param s.w stratum proportion in population
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b  shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return dataframe with empirical rr for N times in each wave
simulated.RR<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  # Response propensity samples from wave 0 to the end
  lst.w<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)

  # Stratum Weight matrix with one row
  w<- t(s.w$stratum.allocation)

  # Calculate
  rr.w<- lapply(lst.w,function(x){
    # transform df to matrix
    x<- as.matrix(x)
    # RR for each iteration
    res<- w%*%x})
  # list to df
  rr<- as.data.frame(do.call(rbind,rr.w))
  return(rr)
}

#' Empirically compute coefficient of variation of response propensity
#' for each wave. The basis of it is on simulated weighted response rate
#' and sampling draws from posterior response propensity.
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat  historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param s.w stratum proportion in population
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return dataframe with N replicated draws with relative to one wave
simulated.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  # Simulated response propesnity list by waves
  rho_sim<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)

  # df RR simulation in waves
  rr_sim<- simulated.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref,N)

  # Stratum weight matrix
  w<- t(s.w$stratum.allocation)

  # identify matrix to expand simulated wave RR so as to have the same rows as simulated Rho in one wave
  I<- matrix(1,nrow = nrow(n.dat),ncol = 1)

  # the number of waves by the name of simulated rho list
  t<- length(names(rho_sim))
  # Do computation by wave indicators
  cv_sim<- lapply(1:t,function(i){
    # simulated rho in wave i
    rho_i<- as.matrix(rho_sim[[i]])

    # simulated rr in wave i
    rr_i<- as.matrix(rr_sim[i,])
    # extend simulated rr
    rr_i_ext<- I%*%rr_i

    # numerator of CV
    x<- sqrt(w%*%((rho_i-rr_i_ext)^2))
    # denominator = rr_i
    res<- as.data.frame(x/rr_i)
    return(res)})

  # lst to df
  cv<- do.call(rbind,cv_sim)
  return(cv)
}
