#' Empirical samples from response propensity
#'
#' Draw samples from power prior and updated posterior and later used to
#' evaluate quality indicators.
#'
#' @param n.dat new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return a list from wave 0 to the end
rho.sample<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref,N=10000){
  p.para<- pprior(n.hat,h.dat,s.score,a,b,svy.ref)
  res.0<- as.data.frame(do.call(rbind,lapply(1:nrow(p.para),
                                             function(i) stats::rbeta(N, shape1 = p.para[i,"shape1"],
                                                                      shape2 = p.para[i,"shape2"], ncp = 0))))
  colnames(res.0)<- c(1:N)

  post.para<- pposterior(n.dat,h.dat,s.score,a,b,svy.ref)
  lst.w<- lapply(split(post.para,post.para$Wave),function(i){
                   res<- as.data.frame(do.call(rbind,
                                               lapply(1:max(i$Strata),function(j) stats::rbeta(N, shape1 = i[j,"shape1"],
                                                                                                     shape2 = i[j,"shape2"],ncp=0))))
                   colnames(res)<- c(1:N)
                   return(res)
                   }
                 )

  rho.w0<- list(w0=res.0)
  rho.w<- lst.w

  lst<- c(rho.w0,rho.w)
  return(lst)
}

