#' The root of mean square error is derived from true response propensity, posterior mean
#' and posterior variance of indicators, weighted response rate and coefficient of variation
#' of response propensity.
#'
#' @param simulated.fun function generates simulated indicator, simulated.RR or simulated.CV
#' @param true.fun function to generated true indicator, true.RR or true.CV, and even adj.true.CV
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param s.w stratum proportion in population
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param N the number of observations (replication) with default 10000
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param adjust whether true indicator is adjusted or not. The default value is no.
#'
#' @return dataframe with posterior mean, variance, quantiles and resulting RMSE with regards to waves
#' @export
sumstats<- function(adjust="no",simulated.fun,true.fun,n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  # survey quality indicator prediction denoted by theta
  theta_sim<- as.matrix(simulated.fun(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  # expectation
  theta_sim_m<- apply(theta_sim,1,mean)
  # variance
  theta_sim_v<- apply(theta_sim,1,var)
  # quantiles
  theta_sim_q<- apply(theta_sim,1,quantile,probs = c(0.025, 0.975))

  # Summary statistics about measurement: location & spread
  res<- data.frame(Wave=paste0("w",0:(nrow(theta_sim)-1)),Mean=theta_sim_m,Variance=theta_sim_v,Quantile1=theta_sim_q[1,],Quantile2=theta_sim_q[2,])

  # True indicator
  if(adjust=="no"){theta_true<- true.fun(n.dat,h.dat,s.score,s.w,a,b,svy.ref)}
  if(adjust=="yes"){theta_true<- true.fun(n.dat)}
  row.names(theta_true)<- NULL
  # merge it in summary df
  res<- cbind(res,theta_true)

  # RMSE
  rmse_lst<- lapply(1:(nrow(res)-1),function(i){
    # True indicator in wave i
    theta_i<- res[i+1,"True"]
    # Expected indicator in wave i-1
    mean<- res[i,"Mean"]
    # Variance indicator in wave i-1
    variance<- res[i,"Variance"]
    # rmse in wave i
    rmse<- sqrt((theta_i-mean)^2+variance)
    return(rmse)})
  #
  rmse<- do.call(cbind,rmse_lst)
  # Old code calculates the timely rmse instead of cumulative rmse on average
  # If want to output timely rmse, then not cumsum()

  # cumulative sum of rmse and then average by the number up to current wave
  #rmse<- cumsum(rmse)/c(1:(nrow(res)-1))
  # add averaged rmse to res
  res<- cbind(res,RMSE=round(c(0,rmse),digits = 3))

  return(res)
}






#' Partial CV as extra evaluation criterion
#'
#' @param n.dat new survey data sets indexied with wave
#' @param h.dat historical survey data sets
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param N the number of observations (replication) with default 10000
#' @param s.w stratum proportion in population
#'
#' @return matrix
partial.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref,N=10000){
  # RP from wave 0 to final wave
  ## Sampled RP from posteriors
  d_rp<- rho.sample(n.dat,h.dat,s.score,a,b,svy.ref,N)
  ## True RP
  t_rp<- as.matrix(rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref))

  # WRR from wave 0 to final wave
  ## Empirical WRR by Sampled RP
  d_wrr<- simulated.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref,N)
  ## True WRR
  t_wrr<- t(true.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  # Stratum population distribution
  w<- sqrt(s.w$stratum.allocation)
  # Identify matrix to extend WRR
  I<- matrix(1,nrow = nrow(n.dat),ncol = 1)
  # Wave number
  T.len<- length(names(d_rp))
  # Strata numver
  G<- nrow(t_rp)

  # Function to get unconditional partial CV
  fun1<- function(i){
    # Sampled RP in wave i
    d_rp_i<- as.matrix(d_rp[[i]])
    # Empirical WRR in wave i
    d_wrr_i<- as.matrix(d_wrr[i,])
    d_wrr_i_ext<- I%*%d_wrr_i

    # fomula
    res<- w*((d_rp_i-d_wrr_i_ext)/d_wrr_i_ext)
    return(res)}
  # PCV for Sampled RP
  d_cv<- lapply(1:T.len,fun1)
  names(d_cv)<- paste0("w",0:(T.len-1))
  # PCV for True RP
  t_cv<- w*((t_rp-(I%*%t_wrr))/(I%*%t_wrr))


  # Proportion list in waves to derive stratum needed extra effort the most
  d_prop<- lapply(names(d_cv),function(t){
    ## partial CV in wave t
    cv_t<- d_cv[[t]]
    ## for each draw which stratum has largest negative value
    x<- apply(cv_t,2,which.min)
    df<- as.data.frame(table(x))
    ## times out of N draws
    df$prop<- (df$Freq)/N
    return(df)})
  t_prop<- lapply(apply(t_cv,2,which.min),function(i) i)

  # RMSE
  rmse<- lapply(1:(T.len-1),function(t){
    ## Propotion in wave t
    d_prop_t<- d_prop[[t]]
    t_prop_t<- t_prop[[t+1]]
    ## identify vector
    J<- rep(0,nrow(d_prop_t))
    ## if realized stratum exist in sampled stratum in wave t
    if(t_prop_t %in% d_prop_t$x == TRUE){ J[which(d_prop_t$x == t_prop_t)]<- 1}
    if(t_prop_t %in% d_prop_t$x == FALSE){ J<- J}
    ## Bias term
    bias<- (sum(d_prop_t$prop-J)^2)/G
    ## Variance term
    variance<- (sum((d_prop_t$prop)*(1-d_prop_t$prop)))/G
    ## rmse
    res<- sqrt(bias + variance)
    return(res)})
  names(rmse)<- paste0("w",1:(T.len-1))
  # list to df
  rmse<- do.call(rbind,rmse)


  return(rmse)
}
