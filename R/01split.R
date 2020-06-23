#' Reshape survey data sets with respondents and sample size separately.
#'
#' For new survey data sets are indicayed by wave, the iteration of new sample released.
#' For historical survey data sets are indicated by the historical names.
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
