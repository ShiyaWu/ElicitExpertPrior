#' Calculate true rho based on observed respondents and sample size from wave 0 to the end of new survey data collection
#'
#' @param n.dat new survey data sets indexied with wave
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return data frame with true response propensity from wave 0 to data collection end
rho.true.w<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref){
  # Aggregated historical survey datasets by feature weights
  p.para<- pprior(n.hat,h.dat,s.score,a,b,svy.ref)
  #dat.w0<- data.frame(w0=p.para$shape1/(p.para$shape1+p.para$shape2))
  dat.w0<- data.frame(w0=(p.para$shape1-1)/(p.para$shape1+p.para$shape2-2))

  # New survey data sets with wave
  lst<- redat(n.dat)
  r<- lst$r
  n<- lst$n
  # True rho with wave
  rho<- r/n
  colnames(rho)<- paste0("w",1:ncol(r))

  # Combine w0 with new waves
  res<- cbind(dat.w0,rho)
  return(res)
}

#' Weighted response rate for new survey observations is calculated. This is used as a benchmark
#' in comparison with empirical weighted response rate calculated by simulated posterior response
#' propensity in wach wave.
#'
#' @param n.dat new survey data sets indexied with wave
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param s.w stratum proportion in population
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return dataframe with wave rows and one real rr
true.RR<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref){
  # True observed response propensity
  rho_true<- rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref)

  # Stratum Weight matrix with one row
  w<- t(s.w$stratum.allocation)

  # transform df to matrix
  rho_true<- as.matrix(rho_true)

  # Weighted response rate
  rr_true<- as.data.frame(t(w%*%rho_true))
  colnames(rr_true)<- "True"
  return(rr_true)
}


#' Coefficeint of variatino of response propensity is specified by
#' true weighted response rate and true response propensity
#'
#' @param n.dat new survey data sets indexied with wave
#' @param s.score similarilty scores to weight hostorical surveys. The choice is related with the degree
#' of feauture importance to response, either equal influence or unequal influence.
#' @param s.w stratum proportion in population
#' @param a shape1 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param b shape2 parameter of beta prior to historical survey data sets. The default value is 1.
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return dataframe
true.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref){
  # True response propensity with relative to wave
  rho_true<- as.matrix(rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref))

  # True weighted response rate overall strata
  rr_true<- t(true.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  # the number of waves
  t<- ncol(rho_true)
  # Stratum weight
  w<- t(s.w$stratum.allocation)

  # Identity matrix
  I<- matrix(1,nrow = nrow(rho_true),ncol = 1)
  # numerator of cv
  res<- sqrt(w%*%(rho_true-(I%*%rr_true))^2)
  # denominator = rr_true
  cv<- as.data.frame(t(res/rr_true))
  colnames(cv)<- "True"
  return(cv)
}



#' True coefficient of variation of response propensity is adjusted to produce the new true results.
#'
#' @param n.dat new survey data sets indexied with wave
#'
#' @return dataframe from wave 0 to the end
adj.true.CV<- function(n.dat){
  # Input microdata file containing records of all parents in the sample for a time portion
  # including their stratum number and 0-1 response
  bias.adjust<- function(i){
    # file path
    path<- paste0("D:/R/Project1/Input/tp_",i,".csv")
    # import microdata
    mic_data<- read.csv(path)[,-1]

    # Add resp column
    mic_data$resp<- rep(0,nrow(mic_data))
    # Revise resp column with 1 when code result is 30 and remains 0 otherwise
    mic_data[mic_data$code_eindresultaat==30,"resp"]<- 1
    #mic_data[mic_data$code_eindresultaat!=30,"resp"]<- 0

    # Revise stratum class to factor
    mic_data$Stratum<- as.factor(mic_data$Stratum)

    # CV under adjustment
    indicator<- getRIndicator(resp ~ Stratum, mic_data, withPartials = FALSE)
    cv<- indicator$CV

    # export
    res<- data.frame(True=cv)
    row.names(res)<- paste0("w",i)

    return(res)}

  # total waves of new survey
  t<- ncol(n.dat[,-1])/2

  # CV in wave 0
  CV_w0<- data.frame(True=0)
  row.names(CV_w0)<- "w0"

  # list for adjusted cv in each wave
  CV_adj_w<- lapply(c(1:t),bias.adjust)
  # convert list to df
  CV_adj<- do.call(rbind,CV_adj_w)
  res<- rbind(CV_w0,CV_adj)

  return(res)
}



#' The indicators are calculated for the whole microdata before breaking it down to waves
#'
#' @param whole.dat original new survey data set without considering waves
#' @param s.w stratum proportion in population
#'
#' @return dataframe with overall RR, CV and adjusted CV
overall.true.theta<- function(whole.dat,s.w){
  data<- whole.dat
  # remove NA
  data<- data[!is.na(data$code_eindresultaat),]

  # revise stratum class
  data$Stratum<- as.factor(data$Stratum)

  # total strata number
  g<- nlevels(data$Stratum)

  # Derive respondent number and sample size relative to each stratum
  fun1<- function(i){
    # sample size
    n<- nrow(data[data$Stratum==i,])
    # respondents
    r<- nrow(data[data$Stratum==i&data$code_eindresultaat==30,])
    # true response propensity
    rho<- r/n

    res<- data.frame(stratum=i,rho=rho)
    return(res)}
  # do computation
  lst<- lapply(1:g,fun1)
  # convert lst to df
  df<- do.call(rbind,lst)
  df$stratum<- factor(df$stratum)

  # stratum proportion
  w<-  t(s.w$stratum.allocation)

  # True RR
  rr<- t(w%*%(as.matrix(df$rho)))

  # identity matrix
  I<- matrix(1,nrow = nrow(df),ncol = 1)

  # True CV
  cv0<-  (sqrt(w%*%(as.matrix(df$rho)-(I%*%rr))^2))/rr

  # adjust CV
  data$resp<- rep(2,nrow(data))
  data[data$code_eindresultaat==30,"resp"]<- 1
  data[data$code_eindresultaat!=30,"resp"]<- 0
  indicator<- getRIndicator(resp ~ Stratum, data, withPartials = FALSE)
  cv1<- indicator$CVUnadj


  result<- data.frame(RR=rr,CV_noadj=cv0,CV_unadj=cv1)
  return(result)
}

