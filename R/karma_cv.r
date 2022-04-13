#' Fitted model validation and cross-validation (and related plots).   
#'
#' @param karma0 Object of class "karma.fit" or "ARIMA".
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param fixed Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30), positive integer for "window" or "percentage" test_type; "auto" to read from karma.fit object or generate; negative integer value to set window size to a multiple of the series' frequency.
#' @param test_type Train-test split type, i.e. percentage or fixed window; "auto": will try to read from karma.fit object or generate; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series; "auto" if input series is a ts() object, test_type is set to "window" and test_pct is set to twice the frequency of the series - if test_pct is given a negative factor, then test_pct (window size) will be set to the frequency of the series times the absolute value of that negative number.
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection.
#' @param cv Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation.
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Compatible with "forecast" package objects: Arima, auto.arima, nnetar, ets, bats, tbats, baggedETS, HoltWinters
#' # Using auto.karma():
#' karma.cv(auto.karma(y), cv="in")  # in-sample forecast
#' karma.cv(auto.karma(y), cv="out")    # out-of-sample forecast
#' # Using karma.fit()
#' karma.cv(karma.fit(y, order=c(2,1,4))); 
#' # Using karma.fit() with fixed terms:
#' karma.cv( karma.fit(y, order=list(c(1,2), 0, c(3,4)), fixed=T) ); 
#' # Using auto.arima():
#' karma.cv(auto.arima(y));
#' #' # Using nnetar():
#' karma.cv(nnetar(y));
#' #' # Using ets():
#' karma.cv(ets(y));


# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


karma.cv <- function(karma0, y = c(), test_pct = "auto", test_type = "auto", log = F, metric = "MAPE", xreg = NULL, fixed = F, plot = T, ensemble = F, cv="out", stdout = F){
  
  #--------- Train-test set cross-validation -----------------
  
  if( is.null(xreg) ){
    if( !is.null(karma0$xreg) ){
      if(stdout == T){
        cat("Found exogenous regressors in karma.fit object.\n")
      }
      xreg = karma0$xreg
    }
  }
  
  if(length(y)==0){
    if(class(karma0)[1]=="bats" | class(karma0)[1]=="tbats" | class(karma0)[1]=="baggedETS" ){
      y = karma0$y
    }else{
      y = karma0$x
    }
  }else{ #do nothing
    #y = as.numeric(y)
  }
  
  if( class(karma0)[1] == "karma.fit" ){
    fixed = karma0$fixed_flag
    log = karma0$log   #warning('karma.cv: "log" parameter is read from karma.fit object; ignoring function argument.')
    karma0$seasonal_terms = karma0$model_seasonal
  }else if (class(karma0)[1] == "forecast_ARIMA" | class(karma0)[1] == "ARIMA"){
    karma0 = karma.cast(karma0)
  }else if( class(karma0)[1] == "Arima" ){
    stop("karma.cv: Please use model of type 'karma.fit' or 'ARIMA'")
  }else if( class(karma0)[1] == "nnetar" | class(karma0)[1] == "ets" | class(karma0)[1] == "bats" | class(karma0)[1]=="tbats" | class(karma0)[1]=="baggedETS" ){
    fixed = F
    log = F
  }else if( class(karma0)[1] == "karma.ensemble" ){
    #Plot ensemble fit: y + forecast
    if(cv == 'in'){
      plot(as.numeric(karma0$x), type='l', col='red', xlab = 'time', ylab = 'y', main="In-sample forecast")
      points(c(karma0$agg_mean_fitted, karma0$agg_mean_forecast), type="l")#, pch=22, lty=2)
      return(karma0$agg_mape_in)
    }    
    else if(cv == 'out'){
      plot(as.numeric(karma0$x), type='l', col='red', xlab = 'time', ylab = 'y', main="Out-of-sample forecast")
      points((length(karma0$agg_mean_fitted) + 1):(length(karma0$x)), karma0$agg_mean_forecast, type="l")#, pch=22, lty=2)
      return(karma0$agg_mape_out)
    }
  }
  
  #Get test size from karma object:
  if( test_type == "auto" & test_pct == "auto" & !is.null(karma0$test_type) & !is.null(karma0$test_pct) ){
    test_type = karma0$test_type
    test_pct = karma0$test_pct
  }
  
  #Get test size from seasonal cycle:
  if( test_type == "auto" ){
    if( class(y) == "ts" ){
      test_type = "window"
      if( frequency(y) > 1 ){
        if( test_pct < 0 ){
          test_pct = abs(test_pct)*frequency(y)
        }else{
          if( test_pct == "auto"){
            test_pct = 2*frequency(y)
          }
        }
        # if( test_pct > 12 & test_pct > length(yt)/3 ){
        #   test_pct = 12
        # }
      }else{
        test_pct = 12        
      }
    }else{
      test_type = "percentage"
      test_pct = 20
    }
  }
  
  
  # Get CV metrics (in-sample):
  if(cv == "in"){
    
    # 'forecast':
    if(log == T){
      y_hat = exp( fitted( karma0 ) )
      y = exp(y)
    }else if (log == F){    
      y_hat = fitted( karma0 )  #'forecast' package
      if(class(karma0)[1] == "nnetar"){
        if( sum(is.na(y_hat)) > 0 ){
          start_ind = sum(is.na(y_hat)) + 1   #set starting index right after the NAs (needed for nnetar mostly)
          y_hat = y_hat[start_ind:length(y_hat)]
          y = y[start_ind:(length(y_hat)+start_ind-1)]
        }
      }
    }
    mape_in = karma.validate(y[1:length(y_hat)], y_hat, metric = metric) #mmetric(y_hat, y[1:length(y_hat)], "MAPE") 
    
    #Plots:
    if(plot == T){
      plot(as.numeric(y_hat), type='l', xlab = 'time', ylab = 'y', main="In-sample forecast")
      points(as.numeric(y), type='l', col='red')
      #legend(1, as.numeric(quantile(y)[4]), c('Predicted','Actual'), cex=0.8, fill=c("black", "red"))
      # #Add prediction intervals:
      # alpha = 0.05
      # moe = qnorm(1-alpha/2) * sd(y_hat) / sqrt(length(y_hat)) #margin of error
      # points(JohnsonJohnson+moe, col='red', type="l", pch=22, lty=2)
      # points(JohnsonJohnson-moe, col='red', type="l", pch=22, lty=2)
    }
    
    return(mape_in)
  }
  
  
  # Get CV metrics (out-of-sample): i.e. cross-validation
  if(cv == "out"){
    
    if(test_type == "percentage" ){
      test_size = round(length(y)*test_pct/100)  #size as percentage of the length of the series
    }else if(test_type == "window" ){ 
      test_size = test_pct  #fixed window size (e.g. in months)
    }
    train_size = length(y)-test_size
    if(class(y) == "numeric" | class(y) == "integer"){
      y_train = y[1:train_size]
    }else if(class(y) == "ts"){  #Arima( y = JohnsonJohnson, order = karma0$model_terms, seasonal = karma0$seasonal_terms)
      y_train = ts( y[1:train_size], start = time(y)[1], end = time(y)[train_size], frequency = frequency(y) )
    }
    
    #Create xreg train and test sets (if exists):
    if(!is.null(xreg)){
      xreg_train = xreg[1:train_size,]  
      xreg_test = xreg[(train_size+1):(train_size+test_size), ]
    }else{
      xreg_train = NULL
      xreg_test = NULL
    }
    
    #-----------------------------------------------------------------------------------
    #Train input model on train-set
    #-----------------------------------------------------------------------------------
    if(fixed == F){
      if(class(karma0)[1] == "nnetar"){
        if(!is.null(xreg)){
          fit_tmp = nnetar( y = y_train, p = karma0$p, P = karma0$P, size = karma0$size, xreg = xreg_train )
        }else{
          fit_tmp = nnetar( y = y_train, p = karma0$p, P = karma0$P, size = karma0$size, xreg = NULL )
        }
      }else if(class(karma0)[1] == "ets"){
        fit_tmp = ets( y = y_train )   #doesn't fit to test set(!): fit_tmp = ets( y = y_train, model = etsfit, use.initial.values = T )
        if(!is.null(xreg)){
          message("karma.cv - message: Exogenous regressors will not be considered in exponential smoothing.")
        }
      }else if(class(karma0)[1] == "bats"){
        fit_tmp = bats( y = y_train )  
        if(!is.null(xreg)){
          message("karma.cv - message: Exogenous regressors will not be considered in exponential smoothing.")
        }
      }else if(class(karma0)[1] == "tbats"){
        fit_tmp = tbats( y = y_train )  
        if(!is.null(xreg)){
          message("karma.cv - message: Exogenous regressors will not be considered in exponential smoothing.")
        }
      }else if(class(karma0)[1] == "baggedETS"){
        fit_tmp = baggedETS( y = y_train )  
        if(!is.null(xreg)){
          message("karma.cv - message: Exogenous regressors will not be considered in exponential smoothing.")
        }
      }else{
        fit_tmp = tryCatch({
          estimation_method = karma0$call$method   #read estimation method from fitted model
          if(is.null(estimation_method) | sum( estimation_method==c("CSS-ML", "ML", "CSS") )==0 ){    #if argument absent from call, then set to default method
            estimation_method = "CSS-ML"  
          }
          if( stdout == T ){
            cat("Attempting estimation method:", estimation_method, "\n")
          }
          
          #----------------Estimate training set model:----------------------------------------------
          if(!is.null(xreg)){
            fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = xreg_train, method = estimation_method )   #fit_tmp = nnetar(y_train) #<- autoregressive ANN
          }else{
            fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = NULL, method = estimation_method )            
          }
          #------------------------------------------------------------------------------------------
          
        }, error=function(e){
          if( stdout == T ){
            cat("WARNING:",conditionMessage(e), "\nFitting model with CSS.")
          }
          fit_tmp = tryCatch({
            if(!is.null(xreg)){
              fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = xreg_train, method = 'CSS' ) 
            }else{
              fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = NULL, method = 'CSS' ) 
            }
          }, error=function(e1){
            if( stdout == T ){
              cat("WARNING:",conditionMessage(e1), "\nFitting model with ML.")
            }            
            fit_tmp = tryCatch({
              if(!is.null(xreg)){
                fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = xreg_train, method = 'ML' ) 
              }else{
                fit_tmp = Arima( y = y_train, order = karma0$model_terms, seasonal = karma0$seasonal_terms, xreg = NULL, method = 'ML' ) 
              }
            },  error=function(e2){
              cat("WARNING:",conditionMessage(e2), "\n")
              msg1 = "ABORTING karma.cv: Arima() failed to fit model with ML-CSS, CSS, and ML methods."
              print(msg1)              
              return(0)
            })
          })
        })  
      }
    }else if(fixed == T){
      if( class(karma0$model_terms) == "list" ){
        carma_terms = karma0$model_terms
        ar_terms = carma_terms[[1]]
        ma_terms = carma_terms[[3]]
        if(length(ar_terms) == 0){
          ar_terms = 0
        }
        if(length(ma_terms) == 0){
          ma_terms = 0
        }
      }
      else{
        stop("For fixed = T, parameter arma order has to be a list of vectors of type 'numeric'.")
      }
      
      fixed_terms = karma.orderbin(ar_terms, ma_terms, diffs = carma_terms[[2]])
      
      fit_tmp = tryCatch({
        estimation_method = karma0$call$method   #read estimation method from fitted model
        if(is.null(estimation_method) | sum( estimation_method==c("CSS-ML", "ML", "CSS") )==0 ){    #if argument absent from call, then set to default method
          estimation_method = "CSS-ML"  
        }
        if( stdout == T ){
          cat("Attempting estimation method:", estimation_method, "\n")
        }
        if(!is.null(xreg)){
          fit_tmp = Arima(y = y_train, order = c(max(ar_terms), carma_terms[[2]], max(ma_terms)), xreg = xreg_train, fixed = fixed_terms, transform.pars=F, method = estimation_method)
        }else{
          fit_tmp = Arima(y = y_train, order = c(max(ar_terms), carma_terms[[2]], max(ma_terms)), xreg = NULL, fixed = fixed_terms, transform.pars=F, method = estimation_method)
        }
      }, error=function(e){
        if( stdout == T ){
          cat("WARNING:",conditionMessage(e), "\nFitting model with CSS.")
        }
        if(!is.null(xreg)){
          fit_tmp = Arima(y = y_train, order = c(max(ar_terms), carma_terms[[2]], max(ma_terms)), xreg = xreg_train, fixed = fixed_terms, transform.pars=F, method='CSS')
        }else{
          fit_tmp = Arima(y = y_train, order = c(max(ar_terms), carma_terms[[2]], max(ma_terms)), xreg = NULL, fixed = fixed_terms, transform.pars=F, method='CSS')
        }
      })          
    }
    
    if(log == T){
      if(!is.null(xreg)){
        y_hat2 = exp(forecast( fit_tmp, h=test_size, xreg = xreg_test )$mean)
      }else{
        y_hat2 = exp(forecast( fit_tmp, h=test_size, xreg = NULL )$mean)
      }
      y_train = exp(y_train)
      y = exp(y)
    }else if(log == F){
      if(!is.null(xreg)){
        y_hat2 = forecast( fit_tmp, h=test_size, xreg = xreg_test )$mean
      }else{
        y_hat2 = forecast( fit_tmp, h=test_size, xreg = NULL )$mean
      }
    }
    
    #Plots:
    if(plot == T){
      plot(as.numeric(c(y_train, y_hat2)), type='l', xlab = 'time', ylab = 'y', main="Out-of-sample forecast")   #plot( forecast( fit_tmp, h=test_size ) )
      points(as.numeric(y), type="l", col="red")  
      #legend(1, as.numeric(quantile(y)[4]), c('Predicted','Actual'), cex=0.8, fill=c("black", "red"))
      # #Add prediction intervals: p1 = 1.96*sd(fitted(kfit))/sqrt(length(y)); plot(JohnsonJohnson); points(JohnsonJohnson+p1, col='red', type="l", pch=22, lty=2); points(JohnsonJohnson-p1, col='red', type="l", pch=22, lty=2); 
      # alpha = 0.05
      # moe = qnorm(1-alpha/2) * sd(y_hat2) / sqrt(length(y_hat2))   #margin of error
      # points(JohnsonJohnson+moe, col='red', type="l", pch=22, lty=2)
      # points(JohnsonJohnson-moe, col='red', type="l", pch=22, lty=2)      
    }
    #MAPE:
    fit_tmp$mape_out = karma.validate(y[(train_size+1):(train_size+test_size)], y_hat2, metric = metric) #mmetric(y_hat2, y[train_size:(train_size+test_size-1)], "MAPE")  #mmetric(y_hat2, yt[109:120], "MAPE")
    
    #Validation data:
    fit_tmp$y_train = y_train
    fit_tmp$y_test = y[(train_size+1):(train_size+test_size)]
    fit_tmp$y_hat_test = y_hat2
    
    if(ensemble == F){
      return(fit_tmp$mape_out)  #return mape only
    }else if(ensemble == T){   #<- this flag is used in training by karma.ensemble()
      y_actual_oos = y[(train_size+1):(train_size+test_size)] 
      if(!is.null(xreg)){
        forecast_obj = forecast( fit_tmp, h=test_size, xreg = xreg_test )   #forecast() object with mean forecast and prediction intervals
      }else{
        forecast_obj = forecast( fit_tmp, h=test_size, xreg = NULL )
      }
      forecast_obj$predicted_residuals = y_actual_oos - y_hat2   #forecast on test set residuals (out-of-sample residuals)
      forecast_obj$y_actual_oos = y_actual_oos   #test set (out-of-sample y)
      forecast_obj$y_hat2 = y_hat2 #forecast on test set (out-f-sample y-hat)
      forecast_obj$mape_out = fit_tmp$mape_out               #test error (out-of-sample forecast error)
      forecast_obj$mape_in = karma.validate(fit_tmp$y_train, fit_tmp$fitted)    #training error (in-sample forecast error)
      forecast_obj$test_pct_valid = test_pct
      forecast_obj$test_type_valid = test_type
      
      return(forecast_obj)     #return out-of-sample forecast object
    }
  }
  
} #karma.cv(karmas[[1]], y)  or: karma.cv(magic.karma(y));  karma.cv(karma.fit(y, order=c(2,1,4))); karma.cv( karma.fit(y, order=list(c(1,2), 0, c(3,4)), fixed=T) ); karma.cv(auto.arima(y));  
#Usage:
# magic.k = magic.karma(y)
# karma.cv(magic.k, cv="out")
# karma.cv(magic.k, cv="in")

