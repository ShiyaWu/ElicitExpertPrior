#' True response propensities
#'
#' Based on observation realized stratum response propensities are derived
#' and used to calculate true overall response propensity.
#'
#' @param n.dat new survey data sets at wave levels
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return data frame with true response propensity from wave 0 to data collection end
rho.true.w<- function(n.dat,h.dat,s.score,a=1,b=1,svy.ref){
  p.para<- pprior(n.hat,h.dat,s.score,a,b,svy.ref)
  dat.w0<- data.frame(w0=(p.para$shape1-1)/(p.para$shape1+p.para$shape2-2))

  lst<- redat(n.dat)
  r<- lst$r
  n<- lst$n
  rho<- r/n
  colnames(rho)<- paste0("w",1:ncol(r))

  res<- cbind(dat.w0,rho)
  return(res)
}

#' True overall RR
#'
#' The weighted response rates as first quality indicator.
#'
#' @param n.dat n.dat new survey data sets at wave levels
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param s.w subgroup proportion
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return dataframe with wave-level true rr
true.RR<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref){
  rho_true<- rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref)
  w<- t(s.w$stratum.allocation)
  rho_true<- as.matrix(rho_true)
  rr_true<- as.data.frame(t(w%*%rho_true))
  colnames(rr_true)<- "True"
  return(rr_true)
}


#' True overall CV
#'
#' The Coefficeint of variation of response propensity as second quality indicator.
#'
#' @param n.dat  new survey data sets at wave levels
#' @param s.score historcial-level similarity scores dependent on feature-level weights
#' @param s.w subgroup proportion
#' @param a shape1 parameter with default value 1
#' @param b shape2 parameter with default value 1
#' @param svy.ref historical reference exists with 1 and 0 otherwise
#' @param h.dat historical survey data sets
#'
#' @return dataframe
true.CV<- function(n.dat,h.dat,s.score,s.w,a=1,b=1,svy.ref){
  rho_true<- as.matrix(rho.true.w(n.dat,h.dat,s.score,a,b,svy.ref))
  rr_true<- t(true.RR(n.dat,h.dat,s.score,s.w,a,b,svy.ref))

  t<- ncol(rho_true)
  w<- t(s.w$stratum.allocation)

  I<- matrix(1,nrow = nrow(rho_true),ncol = 1)
  res<- sqrt(w%*%(rho_true-(I%*%rr_true))^2)

  cv<- as.data.frame(t(res/rr_true))
  colnames(cv)<- "True"
  return(cv)
}



#' True overall CV with bias adjustment
#'
#' Need mircodata.
#'
#' @param n.dat new survey data sets at wave levels
#' @return dataframe from wave 0 to the end
adj.true.CV<- function(n.dat){
  bias.adjust<- function(i){
    path<- paste0("D:/R/Project1/Input/tp_",i,".csv")
    mic_data<- read.csv(path)[,-1]

    mic_data$resp<- rep(0,nrow(mic_data))
    mic_data[mic_data$code_eindresultaat==30,"resp"]<- 1

    mic_data$Stratum<- as.factor(mic_data$Stratum)

    indicator<- getRIndicator(resp ~ Stratum, mic_data, withPartials = FALSE)
    cv<- indicator$CV

    res<- data.frame(True=cv)
    row.names(res)<- paste0("w",i)
    return(res)
    }
  t<- ncol(n.dat[,-1])/2
  CV_w0<- data.frame(True=0)
  row.names(CV_w0)<- "w0"
  CV_adj<- do.call(rbind,lapply(c(1:t),bias.adjust))
  res<- rbind(CV_w0,CV_adj)

  return(res)
}



#' #' The indicators are calculated for the whole microdata before breaking it down to waves
#' #'
#' #' @param whole.dat original new survey data set without considering waves
#' #' @param s.w stratum proportion in population
#' #'
#' #' @return dataframe with overall RR, CV and adjusted CV
#' overall.true.theta<- function(whole.dat,s.w){
#'   data<- whole.dat
#'   # remove NA
#'   data<- data[!is.na(data$code_eindresultaat),]
#'
#'   # revise stratum class
#'   data$Stratum<- as.factor(data$Stratum)
#'
#'   # total strata number
#'   g<- nlevels(data$Stratum)
#'
#'   # Derive respondent number and sample size relative to each stratum
#'   fun1<- function(i){
#'     # sample size
#'     n<- nrow(data[data$Stratum==i,])
#'     # respondents
#'     r<- nrow(data[data$Stratum==i&data$code_eindresultaat==30,])
#'     # true response propensity
#'     rho<- r/n
#'
#'     res<- data.frame(stratum=i,rho=rho)
#'     return(res)}
#'   # do computation
#'   lst<- lapply(1:g,fun1)
#'   # convert lst to df
#'   df<- do.call(rbind,lst)
#'   df$stratum<- factor(df$stratum)
#'
#'   # stratum proportion
#'   w<-  t(s.w$stratum.allocation)
#'
#'   # True RR
#'   rr<- t(w%*%(as.matrix(df$rho)))
#'
#'   # identity matrix
#'   I<- matrix(1,nrow = nrow(df),ncol = 1)
#'
#'   # True CV
#'   cv0<-  (sqrt(w%*%(as.matrix(df$rho)-(I%*%rr))^2))/rr
#'
#'   # adjust CV
#'   data$resp<- rep(2,nrow(data))
#'   data[data$code_eindresultaat==30,"resp"]<- 1
#'   data[data$code_eindresultaat!=30,"resp"]<- 0
#'   indicator<- getRIndicator(resp ~ Stratum, data, withPartials = FALSE)
#'   cv1<- indicator$CVUnadj
#'
#'
#'   result<- data.frame(RR=rr,CV_noadj=cv0,CV_unadj=cv1)
#'   return(result)
#' }
#'
