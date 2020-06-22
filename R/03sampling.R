#' Darw samples from posterior to constrcut empirical response propensity for each stratum
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return a list with wave indicator
rho.sample.w<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref,N=10000){
  # Posterior parameters
  post.para<- pposterior(n.dat,h.dat,s.score,a,b,svy.ref)

  # Divide posterior parameter df by Wave
  lst<- split(post.para,post.para$Wave)
  # Draw sapling from each wave df
  lst.w<- lapply(lst,function(i){
    # Sampling for each row
    lst.g<- lapply(1:max(i$Strata),function(j) stats::rbeta(N, shape1 = i[j,"shape1"], shape2 = i[j,"shape2"],ncp=0))
    # list to df
    res<- as.data.frame(do.call(rbind,lst.g))
    colnames(res)<- c(1:N)
    return(res)})
  return(lst.w)
}



#' Draw samples from power prior to achieve empirical simulation of response propensity.
#' This can be treated as the sampling in wave 0.
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param N the number of observations (replication) with default 10000
#'
#' @return a list
rho.sample.0<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref,N=10000){
  # Power prior parameters
  p.para<- pprior(n.hat,h.dat,s.score,a,b,svy.ref)
  # Draw samples by Strata(row)
  lst.g<- lapply(1:nrow(p.para),function(i) stats::rbeta(N, shape1 = p.para[i,"shape1"], shape2 = p.para[i,"shape2"], ncp = 0))
  #
  res<- as.data.frame(do.call(rbind,lst.g))
  colnames(res)<- c(1:N)
  # put it in a list with name w0
  lst<- list(w0=res)
  return(lst)
}


#' Merge samples from power prior and posterior to one list with indicators from wave 0 to the end
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return a list from wave 0 to the end
rho.sample<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref,N=10000){
  # Sample draws in wave 0
  rho.w0<- rho.sample.0(n.dat,h.dat,s.score,a,b,svy.ref,N)
  # sample draws in new wvaes
  rho.w<- rho.sample.w(n.dat,h.dat,s.score,a,b,svy.ref,N)

  # Combine those list
  lst<- c(rho.w0,rho.w)
}

