#' @title Train ARIMA model for automated univariate time series forecasting
#' @description Improved ARIMA model selection using smart training algorithms and heuristic search options.
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>
#' @param method Model selection algorithm; "markov-selection": Trains MS-ARIMA model; "local-descent": Trains LD-ARIMA model; "random-walk": Trains RW-ARIMA model
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30), positive integer for "window" or "percentage" test_type; "auto" to read from karma.fit object or generate; negative integer value to set window size to a multiple of the series' frequency
#' @param test_type Train-test split type, i.e. percentage or fixed window; "auto": will try to read from karma.fit object or generate; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series; "auto" if input series is a ts() object, test_type is set to "window" and test_pct is set to twice the frequency of the series - if test_pct is given a negative factor, then test_pct (window size) will be set to the frequency of the series times the absolute value of that negative number
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection
#' @param cv Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'
#' @param max_iter For karma.boxjenkins(): Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum)
#' @param max_sdiff Option to set a threshold to the order of seasonal differencing. (T/F) 
#' @param violate_aic Heuristic option to escape local optima by preferring a model of lower AICc even if it's not comparable in theory; (T/F)
#' @param karma_transform Option to automatically detect and apply the right transformation for wide-sense stationarity (T/F)
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast')
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # One-line MS-ARIMA fit:
#' sfit = auto.karma(mdeaths)
#' # Validate trained model on test set:
#' karma.cv(sfit)
#' # One-line LD-ARIMA fit with in-sample validation and verbose search output:
#' sfit = auto.karma(y = mdeaths, method = "local-descent", cv = "in", stdout = T )
#' karma.cv(sfit)
#' # One-line RW-ARIMA using nonparametric training:
#' sfit = auto.karma(y = mdeaths, method = "random-walk", stdout = T )
#' karma.cv(sfit)

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


#ARIMA(p,q,d)(P,Q,D)_m
auto.karma <- function( y, method = "best-first", xreg = NULL, metric = "AICc", cv = "out", max_iter = 5, test_pct = "auto", test_type = "auto", max_sdiff = F, violate_aic = T, karma_transform = F, plot = F, stdout = F ){
  
  period_vec = c(12, 4, 3, 2, 1, 5:11)
  
  if(class(y) == "ts"){
    m = frequency(y)
  }else{
    m_tmp = period_vec[1]
    y = ts(y, frequency = m_tmp)  #ts(y, start=c(2007,5), frequency = 12)
  }
  
  #Adjust test size for higher frequency: (e.g. weekly)
  if(test_type == "auto" & frequency(y) > 20){
    test_type = "percentage"
    test_pct = 20
    max_sdiff = T
    karma_transform = F
    max_iter = 1
  }
  
  #Set initial SARMA parameters:
  sol_local = list()
  sol_local$order = c(0,0,0)
  sol_local$seasonal = c(0,0,0)
  sol_local$eval = 0
  
  #Apply transformation for wide-sense stationarity (if set):
  if( karma_transform == T ){
    sol_local$order[2] = karma.transform(as.numeric(y), stdout = F)$diffs
  }
  
  #Get initial fit:
  sol_local$sfit = karma.fit( y = y, order = sol_local$order, seasonal = sol_local$seasonal, xreg = xreg ) 
  #sol_local$eval = sol_local$sfit$aicc
  if( metric == "AICc" ){
    sol_local$eval = sol_local$sfit$aicc
  }else if( metric == "BIC" ){
    sol_local$eval = sol_local$sfit$bic
  }else if( metric == "AIC" ){
    sol_local$eval = sol_local$sfit$aic
  }else if( metric == "MAPE"){
    sol_local$eval = suppressWarnings( karma.cv(sol_local$sfit, test_pct = test_pct, test_type = test_type, cv=cv, plot = F, xreg = xreg) )
  }  
  
  #Set initial solution to local solution:
  #sol_tmp = sol_local
  
  
  # Greedy search:
  if( method == "local-descent" | method == "better-first" ){               # <----- LD-ARIMA (Local Descent ARIMA) -----
    
    # Going up the hill search:
    counter = 0
    for(i in 1:max_iter){
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "d", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "p", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "q", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "P", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "Q", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "D", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      counter = counter + sol_local$improved
      
      if(plot == T){
        suppressWarnings( karma.cv(sol_local$sfit, cv = cv, test_type = test_type, test_pct = test_pct, xreg = xreg) )
      }
      
      if(counter == 0){
        break
      }else{
        counter = 0
      }
    }
    
  }else if( method == "random-walk" | method == "random-first" ){             #<----- RW-ARIMA (Random Walk ARIMA) ------
    
    #Intialise random variable (throw a dice 6 times):
    rv = sample(6, replace = F)   #random variable
    
    #Shuffle operations:
    dir_vec = c("p", "d", "q", "P", "D", "Q")   #operations
    dir_vec = dir_vec[rv]
    
    #Apply operations:    
    counter = 0
    for(i in 1:max_iter){
      
      for( i in 1:length(dir_vec)){
        sol_local = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = dir_vec[i], direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )  
        if( sol_local$seasonal[2] > 1 ){
          sol_local$seasonal[2] = 1   #compulsory seasonal difference threshold 
        }
        counter = counter + sol_local$improved
      }
      
      if(plot == T){
        suppressWarnings( karma.cv(sol_local$sfit, cv = cv, test_type = test_type, test_pct = test_pct, xreg = xreg) )
      }
      
      if(counter == 0){
        break
      }else{
        counter = 0
      }      
    }
    
  }else if( method == "markov-selection" | method == "best-first" ){       #<----- MS-ARIMA (Markov Selection ARIMA) ------
    
    #Basic greedy best-first search:
    counter = 0
    for( i in 1:max_iter){
      sol_tmp_list = list()
      sol_aicc_vec = c()
      
      sol_tmp_list[[1]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "p", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[1] = sol_tmp_list[[1]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[1]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[1]]$order[2]){
          sol_aicc_vec[1] = Inf
        }
      }      
      
      sol_tmp_list[[2]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "q", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[2] = sol_tmp_list[[2]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[2]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[2]]$order[2]){
          sol_aicc_vec[2] = Inf
        }
      }
      
      sol_tmp_list[[3]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "d", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[3] = sol_tmp_list[[3]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[3]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[3]]$order[2]){
          sol_aicc_vec[3] = Inf
        }
      }
      
      sol_tmp_list[[4]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "P", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[4] = sol_tmp_list[[4]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[4]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[4]]$order[2]){
          sol_aicc_vec[4] = Inf
        }
      }
      
      sol_tmp_list[[5]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "Q", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[5] = sol_tmp_list[[5]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[5]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[5]]$order[2]){
          sol_aicc_vec[5] = Inf
        }
      }
      
      sol_tmp_list[[6]] = GetBestNeighbourSol(y = y, sol_local_ = sol_local, dimension = "D", direction = "+", step_size = 1, metric = metric, test_pct_ = test_pct, test_type_ = test_type, max_sdiff_ = max_sdiff, cv_ = cv, xreg_ = xreg, stdout = stdout )
      sol_aicc_vec[6] = sol_tmp_list[[6]]$eval
      if(violate_aic == F){
        if(sol_local$seasonal[2] != sol_tmp_list[[6]]$seasonal[2] | sol_local$order[2] != sol_tmp_list[[6]]$order[2]){
          sol_aicc_vec[6] = Inf
        }
      }
      
      #Get solution of lowest AICc (or whatever evaluation metric):
      sol_local = sol_tmp_list[[ order(sol_aicc_vec)[1] ]]
      counter = counter + sol_local$improved
      
      if(plot == T){
        suppressWarnings( karma.cv(sol_local$sfit, cv = cv, test_type = test_type, test_pct = test_pct, xreg = xreg) )
      }
      
      if(counter == 0){
        break
      }else{
        counter = 0
      }
    }
    
  }
  
  sol_local$sfit$cv = cv
  sol_local$sfit$test_pct = test_pct
  sol_local$sfit$test_type = test_type
  sol_local$sfit$method = method
  sol_local$sfit$metric = metric
  return(sol_local$sfit)
}





# Find and return best neighbour solution in skarma search:
GetBestNeighbourSol <- function( y, sol_local_, dimension = "p", direction = "+", step_size = 1, metric = "AICc", test_pct_ = 20, test_type_ = "percentage", cv_ = "out", xreg_ = NULL, max_sdiff_ = F, stdout = F ){   #dimension: {p,q,d,P,Q,D,(m)}; direction: {-1, 1}
  
  sol_tmp_ = sol_local_
  
  #Get new neighbour:
  if(dimension == "p"){
    if(direction == "+"){
      sol_tmp_$order[1] = sol_tmp_$order[1] + step_size   #Take +MA step:
    }else if(direction =="-" & sol_tmp_$order[1] > 0){
      sol_tmp_$order[1] = sol_tmp_$order[1] - step_size   #Take -MA step:
    }
  }else if(dimension == "d"){
    if(direction == "+"){
      sol_tmp_$order[2] = sol_tmp_$order[2] + step_size   #Take +d step:
    }else if(direction =="-" & sol_tmp_$order[2] > 0){
      sol_tmp_$order[2] = sol_tmp_$order[2] - step_size   #Take -d step:
    }
  }else if(dimension == "q"){
    if(direction == "+"){
      sol_tmp_$order[3] = sol_tmp_$order[3] + step_size   #Take +AR step:
    }else if(direction =="-" & sol_tmp_$order[3] > 0){
      sol_tmp_$order[3] = sol_tmp_$order[3] - step_size   #Take -AR step:
    }
  }else if(dimension == "P"){
    if(direction == "+"){
      sol_tmp_$seasonal[1] = sol_tmp_$seasonal[1] + step_size   #Take +SMA step:
    }else if(direction =="-" & sol_tmp_$seasonal[1] > 0){
      sol_tmp_$seasonal[1] = sol_tmp_$seasonal[1] - step_size   #Take -SMA step:
    }
  }else if(dimension == "D"){
    if(direction == "+"){
      sol_tmp_$seasonal[2] = sol_tmp_$seasonal[2] + step_size   #Take +D step:
    }else if(direction =="-" & sol_tmp_$seasonal[2] > 0){
      sol_tmp_$seasonal[2] = sol_tmp_$seasonal[2] - step_size   #Take -D step:
    }
  }else if(dimension == "Q"){
    if(direction == "+"){
      sol_tmp_$seasonal[3] = sol_tmp_$seasonal[3] + step_size   #Take +SAR step:
    }else if(direction =="-" & sol_tmp_$seasonal[3] > 0){
      sol_tmp_$seasonal[3] = sol_tmp_$seasonal[3] - step_size   #Take -SAR step:
    }
  }
  
  #Apply seasonal differencing treshold:
  if( max_sdiff_ == T ){
    if( sol_tmp_$seasonal[2] > 1 ){
      sol_tmp_$seasonal[2] = 1
    }
  }
  
  #Get new fit:
  sol_tmp_$sfit = karma.fit( y = y, order = sol_tmp_$order, seasonal = sol_tmp_$seasonal, stdout = stdout, xreg = xreg_ ) 
  
  #Check goodness of fit and increase d term for stationarity: [could alternatively use karma.transform here]
  if(sol_tmp_$sfit$call$method == "CSS"){
    kt = karma.transform(y, stdout = F)
    sol_tmp_$order[2] = kt$diffs
    sol_tmp_$sfit = karma.fit( y = y, order = sol_tmp_$order, seasonal = sol_tmp_$seasonal, log = kt$log, xreg = xreg_ )  #apply logarithmic transformation                    
    #sol_tmp_$sfit = karma.fit( y = y, order = sol_tmp_$order, seasonal = sol_tmp_$seasonal, log = T )  #apply logarithmic transformation                    
    if(sol_tmp_$sfit$call$method == "CSS"){
      if( stdout == T ){
        cat("WARNING -", "ARIMA(", sol_tmp_$order, ")(", sol_tmp_$seasonal, ")[", frequency(sol_tmp_), "]:", "At order 3 of nonseasonal differencing (d=3) and logarithmic transformation, the series does not appear to be stationary (!)\n")
      }
    }
  } #DEBUG: print(sol_tmp_$sfit$call$method)
  
  if( metric == "AICc" ){
    sol_tmp_$eval = sol_tmp_$sfit$aicc
  }else if( metric == "BIC" ){
    sol_tmp_$eval = sol_tmp_$sfit$bic
  }else if( metric == "AIC" ){
    sol_tmp_$eval = sol_tmp_$sfit$aic
  }else if( metric == "MAPE"){
    sol_tmp_$eval = suppressWarnings( karma.cv(sol_tmp_$sfit, test_pct = test_pct_, test_type = test_type_, cv=cv_, plot = F, xreg = xreg_) )
  }
  
  #Evaluate new fit:
  check_nan = suppressWarnings( is.nan(sum(sqrt(diag(sol_tmp_$sfit$var.coef)))) ) #check if there are NaNs in coefficient S.E.
  if( cv_ == "none" ){
    check_cv = T
  }else{
    check_cv = suppressWarnings( karma.cv(sol_tmp_$sfit, test_pct = test_pct_, test_type = test_type_, cv=cv_, plot = F, xreg = xreg_) < karma.cv(sol_local_$sfit, test_pct = test_pct_, test_type = test_type_, cv=cv_, plot = F, xreg = xreg_) )
  }
  if( is.null(check_cv) | is.nan(check_cv) ){
    check_cv = F
  }
  if( is.na(sol_tmp_$eval) | is.na(sol_local_$eval) ){
    check_cv = F
  }
  
  if( (sol_tmp_$eval < sol_local_$eval) & (check_nan == F) & (check_cv == T) ){
    sol_tmp_$alter_eval = sol_local_$eval
    sol_local_ = sol_tmp_
    sol_local_$improved = T
    if( stdout == T ){
      cat("Solution improved - updated current optimum: ", sol_local_$alter_eval, "-> ", sol_local_$eval, " --- ", as.character(sol_local_$sfit), "\n")
    }
  }else{
    sol_local_$alter_eval = sol_tmp_$eval
    sol_local_$improved = F
    # if( stdout == T ){
    #   cat("Solution NOT improved - kept current optimum: ", sol_local_$eval, "\n")
    # }
  }
  
  return(sol_local_)
  
}


