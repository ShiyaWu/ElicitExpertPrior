#'Given historical survey data sets, construct shape parameters of power prior to new survey
#' and later will be updated by upcoming data from new survey. The prior to each historical surveys
#' is represented by shape parameters a and b. The default values for them is 1. The special case
#' means the uniform priors considered.
#'
#' @param n.hat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return dataframe includes shape parameters of power prior
pprior<- function(n.hat,h.dat,s.score,a=1,b=1,svy.ref){
  # Input reshaped historical survey datasets lst
  h.lst<- redat(h.dat)
  r.dat<- h.lst$r # Respondents
  n.dat<- h.lst$n # Sample size

  if(svy.ref==1){
    # power prior shape parameters
    para_a<- as.matrix(r.dat)%*%t(s.score) + matrix(a,nrow = nrow(n.dat))
    para_b<- as.matrix(n.dat-r.dat)%*%t(s.score) + matrix(b,nrow = nrow(n.dat))
  }

  if(svy.ref==0){
    # no prior information to new survey
    para_a<- matrix(a,nrow = nrow(n.dat))
    para_b<- matrix(b,nrow = nrow(n.dat))
  }

  # additonal information to shape parameters
  beta.para<- data.frame(shape1=para_a,shape2=para_b)

  return(beta.para)
}


#' Posterior shape parameters by updating the prior of wave t to posterior of wave t.
#' The prior of wave t is equivalent to the posterior of wave t-1
#' at the end of data collection in wave t.
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return data frame with shape 1, shape 2, strata and wave
pposterior<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref){
  # power prior parameters
  para.prior<- pprior(n.dat,h.dat,s.score,a,b,svy.ref)
  # shape 1
  p.a<- para.prior$shape1
  # shape 2
  p.b<- para.prior$shape2

  # New survey data sets
  dat.lst<- redat(n.dat)
  dat.r<- dat.lst$r
  dat.n<- dat.lst$n
  # Cumulative sum of success(response) and failure(non-response)
  fun1<- function(dat.x){
    # transfer df to matrix
    dat.x<- as.matrix(dat.x)
    res<- t(apply(dat.x,1,cumsum))
    # transfer matrix to df
    res<- as.data.frame(res)
    return(res)}
  # Cumulative response
  cum.r<- fun1(dat.r)
  # Cumulative non-response
  cum.nr<- fun1(dat.n-dat.r)

  # posterior parameters
  fun2<- function(p.x,p.y,cum.x,cum.y){
    w<- ncol(cum.x)
    lst<- lapply(1:w,function(i){
      x<- cum.x[,i]+p.x
      y<- cum.y[,i]+p.y
      # Add wave & strata group
      res<- data.frame(shape1=x,shape2=y,Strata=c(1:nrow(cum.y)),Wave=rep(paste0("w",i),nrow(cum.x)))
      return(res)})
    # list to df
    para<- do.call(rbind,lst)
    return(para)}
  #
  post.para<-fun2(p.a,p.b,cum.r,cum.nr)

  # return posterior para
  return(post.para)
}
