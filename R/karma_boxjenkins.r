#' Box-Jenkins method for ARMA model selection on a stationary process.  
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param diffs Differencing step: Indicates whether the input series needs to be differenced for stationarity (and to what degree); {0,1,...,n}; type <int>.
#' @param log Logarithmic transformation flag. Indicates whether the input series needs to be log-transformed for stationarity; {T, F}; type <logical>.
#' @param fixed Fixed term flag. Indicates whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param box_test T/F flag. Indicates whether or not a Box-Pierce test for autocorrelation should be performed at every algorithm iteration.
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param N Maximum lag at which to calculate autocorrelation and partial autocorrelatin functions; see documentation for acf(), pacf(). 
#' @param max_ar Maximum AR term (value of p).
#' @param max_ma Maximum MA term (value of q).
#' @param max_conv Maximum number of iterations without improvement before the algorithm converges forcefully (stuck to a local optimum).
#' @param max_iter Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Find transformation steps required for stationarity:
#' stationarity.options <- karma.transform(ldeaths, stdout = F, autolog = F, autodiffs = 1)
#' # Apply Box-Jenkins method on the stationary series:
#' boxj.fit <- karma.boxjenkins(ldeaths, diffs = stationarity.options$diffs, log = stationarity.options$log)
#' # Apply cross-validation and calculate MAPE on out-of-sample (test) data:
#' karma.cv(boxj.fit)

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


#library(forecast)

# ----- Automatic model selection -----------------------

karma.boxjenkins <- function(y, diffs = 0, log = F, fixed = F, xreg = NULL, N = 100, box_test = F, max_ar = 20, max_ma = 20, max_conv = 2, max_iter = 200, plot = F, stdout = T){
  
  if( is.null(N) ){
    N = length(y)  #round(10*log10(length(y)))  #as per acf()
  }  
  
  #AC and PAC correlograms: [can get all values using print(acf(yt_req, lag.max=N)) ]
  acf_obj = acf(y, lag.max=N, plot = plot)   #AC object for MA terms
  acf_ma = acf_obj$acf[-1]  #keep only AC values minus the the 0th one.
  pacf_obj = pacf(y, lag.max=N, plot = plot)  #PAC object for AR termss
  pacf_ar = pacf_obj$acf  #keep only PAC values minus the the 0th one.
  
  #Get lag values k for significant AC: (for MA models) 
  ac_cutoff = abs(qnorm((1-0.95)/2)) /sqrt(length(acf_ma))  #AC cutoff (value of dashed lines in R) as per https://stats.stackexchange.com/questions/17760/dashed-lines-in-acf-plot-in-r
  ac_lags = which(abs(acf_ma) >= ac_cutoff )  #lags k with significant AC
  
  #Get lag values k for significant PAC: (for AR models)
  pac_cutoff = abs(qnorm((1-0.95)/2)) /sqrt(length(pacf_ar))   #PAC cutoff (value of dashed lines in R) 
  pac_lags = which(abs(pacf_ar) >= pac_cutoff)  #lags k with significant PAC
  
  
  
  #Apply classic BOX-JENKINS method:
  
  ac_lags_curr = ac_lags
  pac_lags_curr = pac_lags
  arma_model_curr = c()
  karma_terms = c(0,0,0)   #arma order (fixed=F)
  ma_terms = c()     #individual ma terms (fixed=T)
  ar_terms = c()     #individual ma terms (fixed=F)
  iter = 0  #selection algorithm iteration counter
  conv = 0 #counter of times the ARMA terms were forced to a lower value 
  warning_msg = c() #log key warning messages
  
  while(length(ac_lags_curr) != 0 | length(pac_lags_curr) != 0){
    
    iter = iter + 1 #increment interations counter
    
    if(stdout == T){
      cat('Karma Box-Jenkins iterations:', iter, '\n')
    }
    
    # Encode empty lag lists: (for the if statements)
    if(length(ac_lags_curr) == 0){
      ac_lags_curr = Inf
    }
    if(length(pac_lags_curr) == 0){
      pac_lags_curr = Inf
    }
    
    
    # Decide which ARMA terms to add:
    if(ac_lags_curr[1] != Inf & ac_lags_curr[1] < pac_lags_curr[1]){ #if AC lag is smaller than PAC lag
      if(karma_terms[3] != ac_lags_curr[1]){ #check if new MA term does not already exist
        karma_terms[3] = ac_lags_curr[1]  #get max MA term (fixed=F)
        ma_terms = c(ma_terms, ac_lags_curr[1]) #append new MA term (fixed=T)
      }else{
        if( length(ac_lags_curr)>1 ){
          karma_terms[3] = ac_lags_curr[2]
          ma_terms = c(ma_terms, ac_lags_curr[2])
        }
      }
    }else if(pac_lags_curr[1] != Inf & pac_lags_curr[1] < ac_lags_curr[1]){  #if PAC lag is smaller than AC lag
      if(karma_terms[1] != pac_lags_curr[1]){ #check if new AR term does not already exist
        karma_terms[1] = pac_lags_curr[1]     #get max AR term (fixed=F)
        ar_terms = c(ar_terms, pac_lags_curr[1]) #append new AR term (fixed=T)
      }else{
        if( length(pac_lags_curr)>1 ){
          karma_terms[1] = pac_lags_curr[2]
          ar_terms = c(ar_terms, pac_lags_curr[2])
        }
      }
    }else if(ac_lags_curr[1] != Inf & pac_lags_curr[1] != Inf & ac_lags_curr[1] == pac_lags_curr[1]){  #if AC lag is equal to PAC lag
      # Flip a coin
      if(runif(1, min=0, max=1) > 0.5){    # Heads -> add AC
        if(karma_terms[3] != ac_lags_curr[1]){ #check if new MA term does not already exist
          karma_terms[3] = ac_lags_curr[1]
          ma_terms = c(ma_terms, ac_lags_curr[1])
        }else{
          if( length(ac_lags_curr)>1 ){
            karma_terms[3] = ac_lags_curr[2]
            ma_terms = c(ma_terms, ac_lags_curr[2])
          }
        }
      }else{    # Tails -> add PAC
        if(karma_terms[1] != pac_lags_curr[1]){ #check if new AR term does not already exist
          karma_terms[1] = pac_lags_curr[1] 
          ar_terms = c(ar_terms, pac_lags_curr[1])
        }else{
          if( length(pac_lags_curr)>1 ){
            karma_terms[1] = pac_lags_curr[2]
            ar_terms = c(ar_terms, pac_lags_curr[2])
          }
        }      
      }
    }
    
    # Check max size of ARMA terms (if set):
    if(max_ar > 0 & karma_terms[1] > max_ar){
      karma_terms[1] = max_ar
      conv = conv + 1 #increment convergence counter
    }
    if(max_ma > 0 & karma_terms[3] > max_ma){
      karma_terms[3] = max_ma
      conv = conv + 1 #increment convergence counter
    }
    
    # Fit ARIMA model with current ARMA terms:
    if( is.na(diffs) | diffs <= 0 ){
      diffs = 0
      karma_terms[2] = diffs
    }    
    if( diffs > 0 ){
      karma_terms[2] = diffs
    }
    if( fixed == F ){
      if(stdout == T){
        cat('---> Fiting ARMA model of order:', karma_terms, '\n')
      }      
      carma0 = karma.fit(y, order=karma_terms, log = log, xreg = xreg, fixed=F)
    }else if ( fixed == T ){
      ar_terms = unique(ar_terms)[ unique(ar_terms) <= karma_terms[1] ] #adjust non-fixed terms to fixed terms (change that later by separating fixed term creation via hill climbing)
      ma_terms = unique(ma_terms)[ unique(ma_terms) <= karma_terms[3] ]      
      if(stdout == T){
        cat('---> Fiting ARMA model of fixed order:', 'ar={', ar_terms, '}', ' diff={', karma_terms[[2]], '} ma={', ma_terms, '}', '\n')
      }      
      carma0 = karma.fit(y, order=list(ar_terms, karma_terms[2], ma_terms), log = log, xreg = xreg, fixed=T)
    }
   
    # Get AC/PAC values for all lags:
    #AC and PAC correlograms: 
    acf_obj_curr = acf(carma0$residuals, lag.max=N, plot = plot)   #AC object for MA terms
    acf_ma_curr = acf_obj_curr$acf[-1]  #keep only AC values minus the the 0th one.
    pacf_obj_curr = pacf(carma0$residuals, lag.max=N, plot = plot)  #PAC object for AR termss
    pacf_ar_curr = pacf_obj_curr$acf  #keep only PAC values minus the the 0th one.
    
    #Get lag values k for significant AC: (for MA models) 
    ac_cutoff_curr = abs(qnorm((1-0.95)/2)) /sqrt(length(acf_ma_curr)) 
    ac_lags_curr = which(abs(acf_ma_curr) >= ac_cutoff_curr )  #lags k with significant AC
    #Get lag values k for significant PAC: (for AR models)
    pac_cutoff_curr = abs(qnorm((1-0.95)/2)) /sqrt(length(pacf_ar_curr))
    pac_lags_curr = which(abs(pacf_ar_curr) >= pac_cutoff_curr)  #lags k with significant PAC
    
    
    # Decode empty lag lists: (for the while loop)
    if(length(ac_lags_curr) == -1){
      ac_lags_curr = c()
    }
    if(length(pac_lags_curr) == -1){
      pac_lags_curr = c()
    }  
    
    #Exit while() if box_test option is set and the Box test does not reject H0 (i.e. no significant AC):
    if( box_test == T ){
      box_pval = karma.portmanteau(carma0$residuals, plot = F)$box_test_pval
      if( box_pval > 0.05 ){
        if( stdout == T ){
          print("Box-test null hypothesis was not rejected - no significant overall AC detected - selecting current model as optimal.")
        }
        break
      }
    }
    #Exit while() if forced convergence threshold has been reached:
    if( ((max_ar > 0) | (max_ma == 0)) & (conv == max_conv) ){ 
      msg1 = 'Karma Box-Jenkins algorithm did not converge naturally <- forced convergence applied <- ARMA term max threshold reached.'
      warning(msg1)
      warning_msg = c(warning_msg, msg1)
      break
    }
    #Force convergence after max_ter iterations:
    if( iter >= max_iter ){ 
      msg2 = 'Karma Box-Jenkins algorithm did not converge naturally <- forced convergence applied <- number of maximum iterations reached.'
      warning(msg2)
      warning_msg = c(warning_msg, msg2)      
      break
    }
    
  }
  
  if(!exists('carma0')){
    msg3 = 'Karma - WARNING: No adequate model could be fit - the series is white noise - exiting with NULL object.'
    warning(msg3)
    warning_msg = c(warning_msg, msg3)    
    carma0 = NULL
  }
  carma0$model_order = karma_terms
  #carma0$ar_terms = ar_terms
  #carma0$ma_terms = ma_terms
  carma0$iter = iter
  carma0$warning_msg = warning_msg
  #class(carma0) = c("karma", class(carma0))
  
  return(carma0)
  
}



