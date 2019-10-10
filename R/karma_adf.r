#' Augmented Dickey-Fuller test for stationarity.   
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param modeltype Set ADF model with intercept and trend, intercept only, or neither. Takes values: "trend", "drift", "none".
#' @param lags Length of lags for autoregression.
#' @param diffs Differencing step. Indicate whether the input series needs to be differenced for stationarity (and to what degree); {0,1,...,n}; type <int>.
#' @param log Logarithmic transformation flag. Indicate whether the input series needs to be log-transformed for stationarity; {T, F}; type <logical>.
#' @param stdout Option to print out all test diagnostics; <logical>.
#' @return Object of class "karma.adf".
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}, \code{\link{car}}
#' @export
#' @examples
#' # Apply ADF test and print out diagnostics:
#' adf.object <- karma.adf(WWWusage, modeltype = "drift", lags = 1, diffs = 0, log = F, stdout = T)
#' print(adf.object$stationarity)

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

#------------------------------------------------------------------------
#Author: Dimitris Tziotis (2016)
#Arguments:
#modeltype: equiv. to eviews 'include in test equation' - Can be an ADF model with intercept and trend, intercept only, or neither. Takes values: "trend", "drift", "none".
#lags: equiv. to eviews Lag length - User Specified - User sets the length of lags (integer in [0...inf]).
#diffs: equiv. to eviews 1st difference, 2nd difference.
#------------------------------------------------------------------------


# Augmented Dickey Fuller: (adding lags=1 and term dy_ar_k)
#library(car)

karma.adf <- function(y, modeltype="none", lags=1, diffs=0, log=F, dw=F, stdout=T){
  
  #Output object:
  adfout = list(stationarity=F, pval=0, tstat=0, adj_R2=0, model=modeltype, k=lags, diffs=diffs, log=log)

  #Original series:
  adfout$x = y  
    
  #--------------------------------
  #Take logarithm (to remove non-constant variance): 
  
  if(log == T){
    y = log(y)   
    adfout$log = T
  }  
  #--------------------------------   
  
  #--------------------------------
  #Take user-requested differences (to remove trend):
  if( diffs > 0 ){
    y = diff(y, lag=diffs)
  }
  adfout$diffs = diffs
  
  # #Take first difference:
  # if(diffs == 1){
  #   y = diff(y, lag=1)   
  #   adfout$diffs = 1
  # }
  # #Take second difference:
  # if(diffs == 2){
  #   y = diff(y, lag=2)   
  #   adfout$diffs = 2
  # }
  # #Take third difference:
  # if(diffs == 3){
  #   y = diff(y, lag=3)   
  #   adfout$diffs = 3
  # }  
  #--------------------------------
  
  
  #Take ADF difference: (needed for ADF's delta-y)
  dy = diff(y, lag=1)  #delta-y
  
  #DF terms: (basic model terms)---------------------------------
  n = length(dy)
  dy_ar1 = embed(dy, lags+1)[,1]
  y_ar1 = y[(lags+1):n]    #y_[t-1]
  k = lags+1
  #ADF term: (augmented DF term -> multiple lags)
  dy_ar_k = embed(dy, lags+1)[, 2:k]  #<- augmented term y_[t-k] - new term in comparison to the plain: y_[t-k], for k lags
  
  #----------------Apply ADF on selected model:-------------------------
  if(stdout==T){
    cat("\n----------- Augmented Dickey-Fuller Unit Root Test ----------------\n")
  }
  # Model without constant and trend: dyt = b*y_[t-1] + e_t
  if(modeltype=="none"){ 
    #constant term mu = 0 (no coefficient)
    adf = lm(dy_ar1 ~ 0 + y_ar1 + dy_ar_k )     #dyt = b*y_[t-1] + ... + b*y_[t-k]    #<- errors e_t are removed
    b_ind = 1 #index of beta coefficient (needed to take b's t_stat)
  }
  # Model with constant: dyt = mu + b*y_[t-1] + e_t
  if(modeltype=="drift"){ 
    #coefficient added: constant term mu
    adf = lm(dy_ar1 ~ 1 + y_ar1 + dy_ar_k )     #dyt = mu + b*y_[t-1] + ... + b*y_[t-k]    #<- errors e_t are removed
    b_ind = 2 #index of beta coefficient (needed to take b's t_stat)
  }    
  # Model with constant and trend: dyt = mu + b*y_[t-1] + gamma*t + e_t
  if(modeltype=="trend"){ 
    #coefficient added: constant term mu
    temps=(lags+1):n  #trend: gamma*t
    adf = lm(dy_ar1 ~ 1 + y_ar1 + dy_ar_k + temps)     #dyt = mu + b*y_[t-1] + ... + b*y_[t-k] + gamma*t   #<- errors e_t are removed
    b_ind = 2 #index of beta coefficient (needed to take b's t_stat)
  }  
  if(stdout==T){
    print(summary(adf))
  }
  #---------------------------------------------------
  
  #1. Get value for DW (must be close to 2 in order to not have residuals autocorrelation)
  if(dw==T & stdout==T){
    cat("\n-------------- Durbin Watson test ------------------------\n\n")
    print(durbinWatsonTest(adf), lags)  #<-the advantage of not using ur.df() is that we can pass the model object here and get the DW value.
    cat("\n----------------------------------------------------------\n")
  }
  #-----------------------------------------------------  
  
  #2. Get test statistic (test if slope coefficient is equal to 0, i.e. H0: b=0, H1: b<0
  t_stat = summary(adf)$coefficients[b_ind, 3]   #<- At diffs=1, this should return the same value as: #library(urca); ur.df(y,type="none",lags=lags)
  adfout$tstat = t_stat    #pt(t_stat,2)   #significance of b1 (is there a linear relationship in y~y_[t-k] or not)
  adfout$pval = approx( c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15), c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99), t_stat, rule = 2 )$y #<- from adf.test() 
  adfout$adj_R2 = summary(adf)$adj.r.squared
  
  if(stdout==T){
    cat("----------- ADF Unit Root Test results: ------------------\n\n")
    cat("-- ADF Statistic for test H0: b=0, H1: b<0:\n")
    print(t_stat)
    cat("\n-- Critical values:\n")
    print(qnorm(c(.01,.05,.1)/2))
    cat("\n-- Test p-value:\n")   
    print(adfout$pval) 
  }
  
  #library(urca); print(ur.df(y, type=modeltype, lags=lags)) #<-confirm result
  
  #3. Compare test statistic to critical value (t-test on model coefficients): #qnorm(c(.01,.05,.1)/2) #<- If the statistics exceeds those values, then the series is not stationary
  #NB: This output is valid only when DW is almost 2 (i.e. no residuals autocorrelation)
  if(stdout==T){
    cat("\n-- Stationarity diagnostics:")
    cat("\n---------------------------------------------------------------------------------------------\n")
  }
  if(is.numeric(t_stat) & !is.nan(t_stat)){
    #if( adfout$pval > 0.05 ){  
    if( sum(qnorm(c(.01,.05,.1)/2) > t_stat) == 0 ){  #should be compared to: library(tseries); adf.test(y,k=lags) or: summary(ur.df(y,type="none",lags=lags))
      if(stdout==T){
        cat("\nThe series is not stationary (there is no unit root).")
        cat("\nTry differencing to remove trend or log-transformation to achieve constant variance.\n")
      }
      adfout$stationarity = F
    }else{
      if(stdout==T){
        cat("\nThe series is stationary (there is a unit root). You may proceed to ARIMA modelling.\n")
      }
      adfout$stationarity = T
    }
  }else{ #if(!is.numeric(t_stat))
    if(stdout==T){
      cat("\nThe series seems to be missing some crucial data. Aborting.")
    }
    adfout$stationarity = F      
  }
  
  if(stdout==T){
    cat("\n---------------------------------------------------------------------------------------------\n")
  }
  adfout$fit = adf
  adfout$yt = y
  class(adfout) = "karma.adf"
  return(adfout)
  
}

  

