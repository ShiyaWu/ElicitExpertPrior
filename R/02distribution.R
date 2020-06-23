#' Power prior with Beta form
#'
#' Given historical survey data sets, their posterior distribution
#' serves as the prior distribution to a new survey. The prior will be
#' updated to posterior when new data comes in.
#'
#' @param n.hat new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
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


#' Beta posterior distributions
#'
#' \code{pprior} parameters are updated to \code{pposterior} parameter wave by wave
#' and posterior will serve as the prior for next data collection wave.
#'
#' @param n.dat new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#'
#' @return data frame with shape 1, shape 2, strata and wave
pposterior<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref){

  para.prior<- pprior(n.dat,h.dat,s.score,a,b,svy.ref)
  p.a<- para.prior$shape1
  p.b<- para.prior$shape2


  dat.lst<- redat(n.dat)
  dat.r<- dat.lst$r
  dat.n<- dat.lst$n

  fun1<- function(dat.x){
    dat.x<- as.matrix(dat.x)
    res<- t(apply(dat.x,1,cumsum))
    # transfer matrix to df
    res<- as.data.frame(res)
    return(res)
    }
  cum.r<- fun1(dat.r)
  cum.nr<- fun1(dat.n-dat.r)


  fun2<- function(p.x,p.y,cum.x,cum.y){
    w<- ncol(cum.x)
    lst<- lapply(1:w,function(i){
      x<- cum.x[,i]+p.x
      y<- cum.y[,i]+p.y
      res<- data.frame(shape1=x,shape2=y,Strata=c(1:nrow(cum.y)),Wave=rep(paste0("w",i),nrow(cum.x)))
      return(res)
      })
    para<- do.call(rbind,lst)
    return(para)
    }
  post.para<-fun2(p.a,p.b,cum.r,cum.nr)

  return(post.para)
}
