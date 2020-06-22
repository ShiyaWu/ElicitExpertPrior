#' Reshape survey data sets relative to respondents and sample size separately,
#' where for new survey subdata are indexed by waves, i.e. iteration
#' and for historcal surveys by historical surveys name
#'
#' @param svy.dat survey data sets
#'
#' @return a list
redat<- function(svy.dat){
  # remove strata column
  dat<- svy.dat[,-1]
  # group respondent
  r<- dplyr::bind_cols(dat[seq(2,ncol(dat),2)])
  # group sample
  n<- dplyr::bind_cols(dat[seq(1,ncol(dat),2)])
  #
  res<- list(r=r,n=n)
  return(res)
}
