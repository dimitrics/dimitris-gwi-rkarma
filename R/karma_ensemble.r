#' Ensemble learning for time-series: Train an ensemble of weak learners (trained models) that produce a strong learner via the aggregation of multiple out-of-sample forecasts.  
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param nsamples Maximum number of models in the ensemble (unique models will be usually less); type <numeric> integer.
#' @param family Family of ARMA models to choose from; "boxjenkins": non-seasonal ARIMA models with or without fixed terms; "sarima": standard seasonal and non-seasonal ARIMA models.  
#' @param method Box-jenkins model selection algorithm (applicable only when model="boxjenkins"); "greedy": a fully automated karma.boxjenkins in-sample search (default options make it similar to forward selection); "karma": A custom stochastic local search algorithm.
#' @param optimiser Option on the "neighbourhood function" of the optimisation algorithm (applicable only when model="boxjenkins"); "semi-stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate solutions are chosen deterministically; "stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate neighbour solutions are chosen stochastically.
#' @param fixed Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection (applicable only when model="boxjenkins"); {T, F}; type <logical>. 
#' @param box_test T/F flag. Indicates whether or not a Box-Pierce test for autocorrelation should be performed at every algorithm iteration (applicable only when model="boxjenkins").
#' @param autolog Logarithmic search flag. Indicates whether log-transformations on the input series will be part of the search (applicable only when model="boxjenkins").
#' @param autodiffs Differencing search flag. Indicates whether differencing on the input series will be part of the search (applicable only when model="boxjenkins").
#' @param autolags Flag T/F indicating whether or not to set lags automatically as a function of the length of the series (applicable only when model="boxjenkins").
#' @param r2_criterion Flag T/F incidating whether or not to use adjusted R-square as an ADF model selection criterion. When FALSE, the simplest possible stationarity transformation will be preferred (applicable only when model="boxjenkins").
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30), positive integer for "window" or "percentage" test_type; "auto" to read from karma.fit object or generate; negative integer value to set window size to a multiple of the series' frequency.
#' @param test_type Train-test split type, i.e. percentage or fixed window; "auto": will try to read from karma.fit object or generate; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series; "auto" if input series is a ts() object, test_type is set to "window" and test_pct is set to twice the frequency of the series - if test_pct is given a negative factor, then test_pct (window size) will be set to the frequency of the series times the absolute value of that negative number.
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection.
#' @param cv Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation.
#' @param ac_criterion Aucocorrelation / Partial autocorrelation test flag on/off (applicable only when model="boxjenkins"); An optional optimisation constraint which applies portmanteau test on every candidate solution and rejects solutions that do not improve AC/PAC.
#' @param mutations Optional neighbourhood operator (applicable only when model="boxjenkins"); Mutations flag {T, F}: whether or not to apply random "mutations" (term borrowed from evolutionary algorithms) on a candidate solution when the optimiser is about to converge (a way to escape local optima - works somewhat like an inverse simulated annealing).
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param N Maximum lag at which to calculate autocorrelation and partial autocorrelatin functions (applicable only when model="boxjenkins"); see documentation for acf(), pacf(). 
#' @param max_ar Maximum AR term (value of p).
#' @param max_ma Maximum MA term (value of q).
#' @param max_conv For karma.boxjenkins() (applicable only when model="boxjenkins"): Maximum number of iterations without improvement before the algorithm converges forcefully (stuck to a local optimum).
#' @param max_iter For karma.boxjenkins() (applicable only when model="boxjenkins"): Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).
#' @param max_rep For karma-search: Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).  
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @param std_smoothing Option to filter out predicted values if they exceed a user-defined number of standard deviations away from the mean of the ensemble at any time period. If std_smoothing is equal to or less than 0, no std filter is applied. If std_smoothing is set to a positive integer, then the filter threshold will be that input integer times the standard deviation of the ensemble at that period (added to and subtracted from the mean predicted value of the ensemble at that period).
#' @param ci_smoothing Addition (stricter) ensemble filter which leaves out all predicted values that exceed the aggregated confidence interval of the predicted value of the ensemble at that time period; type <boolean> T/F.
#' @return Object of class "karma.ensemble".
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Create ensemble of 10 Box-Jenkins models:
#' kensemble <- karma.ensemble(JohnsonJohnson, nsamples = 10, fixed = T)
#' # Apply cross-validation and calculate MAPE of the aggregated prediction on the out-of-sample data:
#' karma.cv(kensemble)
#' # Forecast and plot 12 periods into the future:
#' karma.forecast(kensemble, h = 12)
#' # All in one line:
#' karma.forecast( karma.ensemble(JohnsonJohnson) )
#' 
#' # Create ensemble of 10 SARIMA models:
#' karma.forecast( karma.ensemble(JohnsonJohnson, family = "sarima" ) )

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


karma.ensemble <- function(y, nsamples = 10, family = "sarima", method = "greedy", optimiser = "semi-stochastic", fixed = F, box_test = F, autolog = F, autodiffs = 1,  autolags = F, r2_criterion = T, test_type = "auto", test_pct = "auto", metric = "MAPE", cv = "out", ac_criterion = F, mutations = F, xreg = NULL, N = 100, max_ar = 20, max_ma = 20, max_conv = 2, max_rep = 1, max_iter = 200, max_sdiff = T, std_smoothing = 1, ci_smoothing = F, plot = T, stdout = T){
  
  if(test_type == "auto" & frequency(y) > 20){
    test_type = "percentage"
    test_pct = 20
    max_sdiff = T
    karma_transform = F
    max_iter = 1
  }else{
    #Detect seasonal cycle:
    if( test_type == "auto" ){
      if( class(y) == "ts" ){
        test_type = "window"
        if( frequency(y) > 1 ){
          if( test_pct < 0 ){
            test_pct = abs(test_pct)*frequency(y)
          }else{
            test_pct = 1*frequency(y)
          }
        }else{
          test_pct = 12
        }
      }else{
        test_type = "percentage"
        test_pct = 20
      }
    }else if( test_type == "percentage" ){
      if(stdout == T){
        cat("Forcing test-set to type 'window' and size 12.\n")
      }
      test_type = "window"
      test_pct = 12
    }
  }
  
  #Initialise ensemble object:
  kensemble = list()
  kensemble$x = y
  kensemble$unique_model_count = 1
  kensemble$model_list = list()
  kensemble$seasonal_list = list()
  kensemble$concat_list = list()
  kensemble$kfit_list = list()
  
  kfit_vec = list()
  mean_fitted_mat = matrix(0, nsamples, length(y)-test_pct) #matrix with average in-sample fitted values (y-hat) of all models
  mean_pred_mat = matrix(0, nsamples, test_pct)   #matrix with the average out-of-sample predictions (y-pred) of all models
  lower_ci_mat = matrix(0, nsamples, test_pct)    #matrix with the lower prediction interval of all models
  upper_ci_mat = matrix(0, nsamples, test_pct)    #matrix with the upper prediction interval of all models
  uniqueness_flag_vec = !logical(nsamples)        #vector of flags (T/F) indidating which models appear for the first time
  #weighted_model_vec = numeric(nsamples)         #weights of bagged models 
  mape_out_vec = c()                              #MAPE or MSE of bagged model
  #pred_residuals_mat = matrix(0, nsamples, test_pct) #matrix with the out-of-sample residuals of all models
  
  for(i in 1:nsamples){
    
    #Fit and validate model:
    if( family == "boxjenkins" ){
      #----------------- Fit Box-Jenkins model: --------------------------------------------------
      kfit_vec[[i]] = tryCatch({
        kfit_vec[[i]] = auto.boxjen(y, method = method, optimiser = optimiser, fixed = fixed, box_test = box_test, autolog = autolog, autodiffs = autodiffs, autolags = autolags, r2_criterion = r2_criterion, test_pct = test_pct, test_type = test_type, metric = metric, cv = cv, ac_criterion = ac_criterion, mutations = mutations, xreg = xreg, N = N, max_ar = max_ar, max_ma = max_ma, max_conv = max_conv, max_rep = max_rep, max_iter = max_iter, plot = F, stdout = F)    
      }, error=function(e2){
        cat("Skipping model - caught fitting error:", conditionMessage(e2), "\n")
        return(0)
      })
      if(as.character(kfit_vec[[i]]) == "0"){
        next
      }
      
      #Reject null model:
      if(sum(kfit_vec[[i]]$model_order) == 0){
        next
      }     
    }else if( family == "sarima" ){
      #----------------- Fit ARIMA/SARIMA model: --------------------------------------------------
      kfit_vec[[i]] = tryCatch({
        kfit_vec[[i]] = auto.karma( y, method = "random-first", xreg = xreg, metric = "AICc", cv = cv, max_iter = max_iter, test_pct = test_pct, test_type = test_type, plot = plot, max_sdiff = max_sdiff, stdout = stdout )
      }, error=function(e2){
        cat("Skipping model - caught fitting error:", conditionMessage(e2), "\n")
        return(0)
      })
      if(as.character(kfit_vec[[i]]) == "0"){
        next
      }
      
      #Reject null model:
      if(sum(kfit_vec[[i]]$model_terms) == 0 & sum(kfit_vec[[i]]$model_seasonal) == 0){
        next
      }     
    }
    
    
    #Validate model:
    kfit_vec[[i]]$cv_obj = karma.cv(kfit_vec[[i]], ensemble = T, test_pct = test_pct, test_type = test_type, plot = plot, metric = metric)
    
    #Store forecast data:
    mean_fitted_mat[i,] = fitted(kfit_vec[[i]]$cv_obj)
    mean_pred_mat[i,] = kfit_vec[[i]]$cv_obj$mean
    lower_ci_mat[i, ] = kfit_vec[[i]]$cv_obj$lower[,1] #[,1]: 80% c.i., [,2]: 95% c.i.
    upper_ci_mat[i, ] = kfit_vec[[i]]$cv_obj$upper[,1]
    #pred_residuals_mat[i, ] = kfit_vec[[i]]$cv_obj$predicted_residuals
    #mape_out_vec[i] = karma.validate(y = kfit_vec[[i]]$cv_obj$y_actual_oos, y_hat = kfit_vec[[i]]$cv_obj$mean, metric = metric )
    mape_out_vec[i] = kfit_vec[[i]]$cv_obj$mape_out
    kfit_vec[[i]]$mape_out = mape_out_vec[i]
    kfit_vec[[i]]$concat_terms = c(kfit_vec[[i]]$model_terms, kfit_vec[[i]]$model_seasonal)
    
    if( stdout == T ){
      if( family == "boxjenkins" ){
        cat(metric, ":", mape_out_vec[i], "    AICc:", kfit_vec[[i]]$aicc)
        cat("\nFitted Box-Jenkins model [", i, "]: \n" )
        print(kfit_vec[[i]]$model_terms)
        cat("\n----------------\n\n")
      }else if( family == "sarima" ){
        cat("\n[ Model ", i, "]   Order: ", kfit_vec[[i]]$model_terms, " ---  Seasonal: ", kfit_vec[[i]]$model_seasonal, " --- ", as.character(kfit_vec[[i]]), "\n")
        cat("\n----------------\n\n")
      }
    }
    
    if(i > 1){  #see if this model has been picked before
      if( family == "boxjenkins" ){
        if( !is_list_member( item_ = kfit_vec[[i]]$model_terms, list_ = kensemble$model_list ) ){
          kensemble$unique_model_count = kensemble$unique_model_count + 1
        }else{
          uniqueness_flag_vec[i] = F  #mark model as innactive
        }
      }else if( family == "sarima" ){
        if( !is_list_member( item_ = kfit_vec[[i]]$concat_terms, list_ = kensemble$concat_list ) ){
          kensemble$unique_model_count = kensemble$unique_model_count + 1
        }else{
          uniqueness_flag_vec[i] = F  #mark model as innactive
        }        
      }
    }
    #Store model data:
    kensemble$model_list[[i]] = kfit_vec[[i]]$model_terms
    kensemble$seasonal_list[[i]] = kfit_vec[[i]]$model_seasonal
    kensemble$concat_list[[i]] = kfit_vec[[i]]$concat_terms
    kensemble$kfit_list[[i]] = kfit_vec[[i]]  
    
  } #end-for
  
  
  #Remove predicted outliers from bad models:
  if( std_smoothing > 0 ){
    #Keep only predicted values that fall within the confidence interval of the mean of all predictions:
    for( j in 1:dim(mean_pred_mat)[2] ){
      # alpha = 0.01 #confidence level
      # moe = qnorm(1-alpha/2) * sd(mean_pred_mat[,j]) / sqrt(length(mean_pred_mat[,j])) #margin of error
      # bound_up = mean(mean_pred_mat[,j]) + moe
      # bound_do = mean(mean_pred_mat[,j]) - moe
      bound_up = mean(mean_pred_mat[,j]) + std_smoothing * sd(mean_pred_mat[,j])
      bound_do = mean(mean_pred_mat[,j]) - std_smoothing * sd(mean_pred_mat[,j])
      outlier_inds = ((mean_pred_mat[,j] > bound_up) | (mean_pred_mat[,j] < bound_do))
      
      if( sum(outlier_inds) > 0 ){
        if( stdout == T ){
          cat("\nkarma.ensemble: Detected", sum(outlier_inds), "predicted outliers <- neutralising.\n")
        }
        mean_pred_mat[ outlier_inds, j ] = mean(mean_pred_mat[!outlier_inds,j])  #neutralise outlisers by assining them to the mean of the non-outliers
      }else if( sum(outlier_inds) == length(outlier_inds) ){
        if( stdout == T ){
          cat("\nkarma.ensemble: Detected too many predicted outliers (!) <- leaving as is.\n")
        }
      }
    }
  }
  
  #Get aggregated forecast data:
  if(sum(uniqueness_flag_vec) > 1){
    kensemble$agg_mean_fitted = apply(mean_fitted_mat[uniqueness_flag_vec, ], 2, mean)
    kensemble$agg_mean_forecast = apply(mean_pred_mat[uniqueness_flag_vec, ], 2, mean)
    kensemble$agg_lower_ci = apply(lower_ci_mat[uniqueness_flag_vec, ], 2, mean)
    kensemble$agg_upper_ci = apply(upper_ci_mat[uniqueness_flag_vec, ], 2, mean)    
  }else if(sum(uniqueness_flag_vec) == 1){
    kensemble$agg_mean_fitted = mean_fitted_mat[uniqueness_flag_vec, ]
    kensemble$agg_mean_forecast = mean_pred_mat[uniqueness_flag_vec, ]
    kensemble$agg_lower_ci = lower_ci_mat[uniqueness_flag_vec, ]
    kensemble$agg_upper_ci = upper_ci_mat[uniqueness_flag_vec, ]
  }else{
    stop("karma.ensemble(): No acceptable models were fitted - exiting...")
  }
  
  kensemble$y_actual_oos = kfit_vec[[i]]$cv_obj$y_actual_oos    #out-of-sample actual (same for all models)
  kensemble$agg_mape_out = karma.validate(y = kensemble$y_actual_oos, y_hat = kensemble$agg_mean_forecast, metric = metric) #aggregated mean out-of-sample mape
  kensemble$agg_mape_in = karma.validate(y = kensemble$x, y_hat = c(kensemble$agg_mean_fitted, kensemble$agg_mean_forecast), metric = metric) #aggregated mean out-of-sample mape
  
  #Smoothen aggregated forecast to fall within aggregated confidence intervals: (heuristic approach, not statistically correct)
  if( ci_smoothing == T ){
    if( sum(kensemble$agg_mean_forecast < kensemble$agg_lower_ci) > 0 ){
      if( stdout == T ){
        cat("\nkarma.ensemble: Ensemble prediction too low <- using lower C.I. smoothening.")
      }
      kensemble$agg_mean_forecast[kensemble$agg_mean_forecast < kensemble$agg_lower_ci] = kensemble$agg_lower_ci[kensemble$agg_mean_forecast < kensemble$agg_lower_ci]
    }
    if( sum(kensemble$agg_mean_forecast > kensemble$agg_upper_ci) > 0 ){
      if( stdout == T ){
        cat("\nkarma.ensemble: Ensemble prediction too high <- using upper C.I. smoothening.")
      }
      kensemble$agg_mean_forecast[kensemble$agg_mean_forecast > kensemble$agg_upper_ci] = kensemble$agg_upper_ci[kensemble$agg_mean_forecast > kensemble$agg_upper_ci]
    }
  }  
  kensemble$agg_pred_residuals = kensemble$y_actual_oos - kensemble$agg_mean_forecast #apply(pred_residuals_mat, 2, mean)
  
  if(plot == T){
    plot(as.numeric(kensemble$x), type='l', col='red', xlab = 'time', ylab = 'y', main="Boosted out-of-sample forecast")
    points((length(kensemble$agg_mean_fitted) + 1):(length(kensemble$x)), kensemble$agg_mean_forecast, type="l")#, pch=22, lty=2)
  }  
  
  kensemble$uniqueness_flag_vec = uniqueness_flag_vec
  kensemble$test_type = test_type
  kensemble$test_pct = test_pct
  
  class(kensemble) = "karma.ensemble"
  return(kensemble)
  
}




#Return T/F whether item_ (list or vector) is member of list_ (list)

is_list_member <- function(item_, list_){
  
  if(class(item_) == "numeric"){
    for(i in 1:length(list_)){
      if(length(item_) == length(list_[[i]])){
        if( sum(item_ != list_[[i]]) == 0 ){   #if there's a match
          return(T)
        }
      }
    }
    return(F)
  }else if(class(item_) == "list"){
    for(i in 1:length(list_)){
      item_bin = karma.orderbin(ar_terms = item_[[1]], ma_terms = item_[[3]], diffs = item_[[2]], format = 'true-false') 
      list_i_bin = karma.orderbin(ar_terms = list_[[i]][[1]], ma_terms = list_[[i]][[3]], diffs = list_[[i]][[2]], format = 'true-false') 
      if(length(item_bin) == length(list_i_bin)){
        if( sum(item_bin == list_i_bin) == length(item_bin) ){   #if there's a match
          return(T)
        }
      }
    }
    return(F)    
  }else{
    stop("Function is_list_member() accepts only vectors and lists as input.")
  }
  
}
#e.g: is_list_member(c(1,2,3), list(c(1,2,3), c(4,5,6)))
#or:  is_list_member(list(c(1,2,3), c(2), c(1,2)), list( list(c(1,2,3,4), c(2), c(1,2)), list(c(1,2,3), c(2), c(1,2)) ) )


