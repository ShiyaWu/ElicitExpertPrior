#' Evaluated RR
#'
#' @param n.dat  new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param s.w subgroup proportion
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @return dataframe with empirical rr for N times in each wave
simulated.RR<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  lst.w<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)
  w<- t(s.w$stratum.allocation)
  rr.w<- lapply(lst.w,function(x) w%*%as.matrix(x))
  rr<- as.data.frame(do.call(rbind,rr.w))
  return(rr)
}

#' Estimated CV
#'
#' @param n.dat  new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param s.w subgroup proportion
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return dataframe with N replicated draws with relative to one wave
simulated.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  rho_sim<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)
  rr_sim<- simulated.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref,N)

  w<- t(s.w$stratum.allocation)
  I<- matrix(1,nrow = nrow(n.dat),ncol = 1)

  t<- length(names(rho_sim))
  cv_sim<- lapply(1:t,function(i){
    rho_i<- as.matrix(rho_sim[[i]])
    rr_i<- as.matrix(rr_sim[i,])
    rr_i_ext<- I%*%rr_i

    res<- as.data.frame(sqrt(w%*%((rho_i-rr_i_ext)^2))/rr_i)
    return(res)
    })

  cv<- do.call(rbind,cv_sim)
  return(cv)
}
