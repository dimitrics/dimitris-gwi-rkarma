#' Find transformation steps required for stationarity.   
#'
#' @param yt A univariate time-series vector; type <numeric> or <ts>.
#' @param autolog Logarithmic search flag. Indicates whether log-transformations on the input series will be part of the search.
#' @param autodiffs Differencing search flag. Indicates whether differencing on the input series will be part of the search.
#' @param autolags Flag T/F indicating whether or not to set lags automatically as a function of the length of the series.
#' @param r2_criterion Flag T/F incidating whether or not to use adjusted R-square as an ADF model selection criterion. When FALSE, the simplest possible stationarity transformation will be preferred.
#' @param stdout Option to print out all search diagnostics; <logical>.
#' @return Object of class "karma.transform".
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Apply ADF test and print out diagnostics:
#' stationarity.options <- karma.transform(WWWusage, stdout = F, autolog = F, autodiffs = 1)
#' print(stationarity.options$model)  # best fitted model type
#' print(stationarity.options$diffs)  # best degree of differencing
#' print(stationarity.options$log)    # best use of logarithmic transformation 

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


#library(rminer)

# ----- Auto stationarity -----------------------

karma.transform <- function(yt, stdout = T, autolog = T, autodiffs = 1, autolags = F, r2_criterion = T){
  
  # Stationarity check:
  adf_vec = list() #all karma.adf model objects
  ar2_vec = c()   #all karma.adf$adj_R2 values
  elig_vec = c()  #all cdf$stationarity values
  
  #Lag length:
  if( autolags == F ){
    k = 1     
  }else if( autolags == T ){
    k = trunc((length(yt) - 1)^(1/3)) #<-from edit(adf.test) #floor( ((7.4 * length(yt))^(1/3))/2 )   #set lag length heuristically according to series length (formula regressed from adf.test)
  }
  
  #Raw series:
  adf_vec[[1]] = karma.adf(yt, modeltype="none", lags=k, diffs=0, log=F, stdout=stdout)
  ar2_vec[1] = adf_vec[[1]]$adj_R2
  elig_vec[1] = adf_vec[[1]]$stationarity
  
  #1st difference series:
  adf_vec[[2]] = karma.adf(yt, modeltype="none", lags=k, diffs=autodiffs, log=F, stdout=stdout)
  ar2_vec[2] = adf_vec[[2]]$adj_R2
  elig_vec[2] = adf_vec[[2]]$stationarity
  
  #Log series:
  adf_vec[[3]] = karma.adf(yt, modeltype="none", lags=k, diffs=0, log=autolog, stdout=stdout)
  ar2_vec[3] = adf_vec[[3]]$adj_R2
  elig_vec[3] = adf_vec[[3]]$stationarity
  
  #Log-difference (1st diff) series:
  adf_vec[[4]] = karma.adf(yt, modeltype="none", lags=k, diffs=autodiffs, log=autolog, stdout=stdout)
  ar2_vec[4] = adf_vec[[4]]$adj_R2
  elig_vec[4] = adf_vec[[4]]$stationarity
  
  #2nd difference series:
  adf_vec[[5]] = karma.adf(yt, modeltype="none", lags=k, diffs=autodiffs*2, log=F, stdout=stdout)
  ar2_vec[5] = adf_vec[[5]]$adj_R2
  elig_vec[5] = adf_vec[[5]]$stationarity
  
  #Log-difference (2nd diff) series:
  adf_vec[[6]] = karma.adf(yt, modeltype="none", lags=k, diffs=autodiffs*2, log=autolog, stdout=stdout)
  ar2_vec[6] = adf_vec[[6]]$adj_R2
  elig_vec[6] = adf_vec[[6]]$stationarity
  
  
  if( r2_criterion == T ){
    #Get indexes of all karma.adf models ranked by adj_r2 value:
    rank_asc = order( -c(ar2_vec[1], ar2_vec[2], ar2_vec[3], ar2_vec[4], ar2_vec[5], ar2_vec[6]) )
  }else if ( r2_criterion == F ){
    #Get indexes from simplest to most complex model:
    rank_asc = c(1, 2, 4, 6, 3, 5)
  }
  
  #Get eligible (stationary) karma.adf model of maximum adj-R2:
  maxR2_ind = rank_asc[elig_vec[rank_asc]][1]
  
  #Print R2 and plot best model series:
  if( is.na(maxR2_ind) ){
    warning('karma.transform(): No default transformation could be found to achieve stationarity; setting higher order differencing (!)')
    model0 = adf_vec[[6]] 
    model0$diffs = 2
  }else{
    model0 = adf_vec[[maxR2_ind]]   #best starting model (as a 'karma.adf' object)
  }
  model0$y = yt
  model0$all_adf_models = adf_vec #all models tested for stationarity
  model0$all_model_rank = rank_asc #all tested models ranked by adj-R2
  model0$all_model_elig = elig_vec #all stationary models eligible for time-series analysis
  model0$all_model_ar2 = ar2_vec
  class(model0) = "karma.transform"
  
  return(model0)
}