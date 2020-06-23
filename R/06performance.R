#' Prediction performance
#'
#' Some statistics, including posterior mean, lower(upper) quantile and RMSE.
#'
#' @param simulated.fun simulated indicators
#' @param true.fun true indicators
#' @param n.dat new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score s.score historcial-level similarity scores dependent on feature-level weights
#' @param s.w subgroup proportion
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param adjust CV with bias adjusted is 1 and 0 otherwise
#'
#' @return dataframe with posterior mean, variance, quantiles and resulting RMSE with regards to waves
#' @export
sumstats<- function(adjust="no",simulated.fun,true.fun,n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  theta_sim<- as.matrix(simulated.fun(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  theta_sim_m<- apply(theta_sim,1,mean)
  theta_sim_v<- apply(theta_sim,1,var)
  theta_sim_q<- apply(theta_sim,1,quantile,probs = c(0.025, 0.975))
  res<- data.frame(Wave=paste0("w",0:(nrow(theta_sim)-1)),Mean=theta_sim_m,Variance=theta_sim_v,Quantile1=theta_sim_q[1,],Quantile2=theta_sim_q[2,])

  if(adjust=="no"){theta_true<- true.fun(n.dat,h.dat,s.score,s.w,a,b,svy.ref)}
  if(adjust=="yes"){theta_true<- true.fun(n.dat)}
  row.names(theta_true)<- NULL
  res<- cbind(res,theta_true)

  rmse_lst<- lapply(1:(nrow(res)-1),function(i){
    theta_i<- res[i+1,"True"]
    mean<- res[i,"Mean"]
    variance<- res[i,"Variance"]
    rmse<- sqrt((theta_i-mean)^2+variance)
    return(rmse)
    })
  rmse<- do.call(cbind,rmse_lst)

  res<- cbind(res,RMSE=round(c(0,rmse),digits = 3))
  return(res)
}






#' Partial CV RMSE
#'
#' @param n.dat new survey data sets at wave levels
#' @param h.dat historical survey data sets
#' @param s.score s.score historcial-level similarity scores dependent on feature-level weights
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param N the number of observations (replication) with default 10000
#' @param s.w subgroup proportion
#'
#' @return matrix
partial.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  d_rp<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)
  t_rp<- as.matrix(rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref))

  d_wrr<- simulated.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref,N)
  t_wrr<- t(true.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  w<- sqrt(s.w$stratum.allocation)
  I<- matrix(1,nrow = nrow(n.dat),ncol = 1)
  T.len<- length(names(d_rp))
  G<- nrow(t_rp)

  fun1<- function(i){
    d_rp_i<- as.matrix(d_rp[[i]])
    d_wrr_i_ext<- I%*%as.matrix(d_wrr[i,])

    res<- w*((d_rp_i-d_wrr_i_ext)/d_wrr_i_ext)
    return(res)
    }
  d_cv<- lapply(1:T.len,fun1)
  names(d_cv)<- paste0("w",0:(T.len-1))
  t_cv<- w*((t_rp-(I%*%t_wrr))/(I%*%t_wrr))

  d_prop<- lapply(names(d_cv),function(t){
    cv_t<- d_cv[[t]]
    x<- apply(cv_t,2,which.min)
    df<- as.data.frame(table(x))
    df$prop<- (df$Freq)/N
    return(df)
    })
  t_prop<- lapply(apply(t_cv,2,which.min),function(i) i)

  rmse<- lapply(1:(T.len-1),function(t){
    d_prop_t<- d_prop[[t]]
    t_prop_t<- t_prop[[t+1]]
    J<- rep(0,nrow(d_prop_t))

    if(t_prop_t %in% d_prop_t$x == TRUE){ J[which(d_prop_t$x == t_prop_t)]<- 1}
    if(t_prop_t %in% d_prop_t$x == FALSE){ J<- J}
    bias<- (sum(d_prop_t$prop-J)^2)/G
    variance<- (sum((d_prop_t$prop)*(1-d_prop_t$prop)))/G

    res<- sqrt(bias + variance)
    return(res)
    })
  names(rmse)<- paste0("w",1:(T.len-1))
  rmse<- do.call(rbind,rmse)

  return(rmse)
}
