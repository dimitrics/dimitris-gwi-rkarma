#' Automated Box-Jenkins method - Improved model selection with flexible optimisation criteria.  
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param method Generic selection algorithm; "greedy": a fully automated karma.boxjenkins in-sample search (default options make it similar to forward selection); "karma": A custom stochastic local search algorithm.
#' @param optimiser Option on the "neighbourhood function" of the optimisation algorithm; "semi-stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate solutions are chosen deterministically; "stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate neighbour solutions are chosen stochastically.
#' @param fixed Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param box_test T/F flag. Indicates whether or not a Box-Pierce test for autocorrelation should be performed at every algorithm iteration.
#' @param autolog Logarithmic search flag. Indicates whether or not log-transformation on the series will be decided algorithmically. Value "force": will force log-transformation.
#' @param autodiffs Differencing search flag. T/F - Indicates whether or not the order of differencing will be algorithmically determined. Values "skip" or F of 0: will force to skip differencing. Negative values will force their absolute value to be set as the order of differencing.
#' @param autolags Flag T/F indicating whether or not to set lags automatically as a function of the length of the series.
#' @param r2_criterion Flag T/F incidating whether or not to use adjusted R-square as an ADF model selection criterion. When FALSE, the simplest possible stationarity transformation will be preferred.
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30).
#' @param test_type Train-test split type, i.e. percentage or fixed window; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series.
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection.
#' @param cv Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation.
#' @param ac_criterion Aucocorrelation / Partial autocorrelation test flag on/off; An optional optimisation constraint which applies portmanteau test on every candidate solution and rejects solutions that do not improve AC/PAC.
#' @param mutations Optional neighbourhood operator; Mutations flag {T, F}: whether or not to apply random "mutations" (term borrowed from evolutionary algorithms) on a candidate solution when the optimiser is about to converge (a way to escape local optima - works somewhat like an inverse simulated annealing).
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param N Maximum lag at which to calculate autocorrelation and partial autocorrelatin functions; see documentation for acf(), pacf(). 
#' @param max_ar Maximum AR term (value of p).
#' @param max_ma Maximum MA term (value of q).
#' @param max_conv For karma.boxjenkins(): Maximum number of iterations without improvement before the algorithm converges forcefully (stuck to a local optimum).
#' @param max_iter For karma.boxjenkins(): Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).
#' @param max_rep For karma-search: Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).  
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Automatic fit: (default: method = "greedy"; box-jenkins on insample CV)
#' magic.fit <- auto.boxjen(ldeaths)   
#' # Apply cross-validation and calculate MAPE on out-of-sample (test) data:
#' karma.cv(magic.fit)

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


auto.boxjen <- function(y, method = "greedy", optimiser = "semi-stochastic", fixed = F, box_test = F, autolog = F, autodiffs = 1, autolags = F, r2_criterion = T, test_pct = 20, test_type = "percentage", metric = "MAPE", cv = "out", ac_criterion = F, mutations = F, xreg = NULL, N = 100, max_ar = 20, max_ma = 20, max_conv = 2, max_rep = 1, max_iter = 200, plot = T, stdout = F){
  
  if( autolog == T ){
    if( sum(y<=0) > 0 ){
      autolog = F
      print("There are zero or negative values in the series - 'autolog' is set to FALSE.")
    }
  }
  if( fixed == T ){
    max_ar = 30
    max_ma = 30
  }
  if( autodiffs == "skip" ){
    autodiffs = F   #force integrated term to 0
  }  
  
  #Find transformation for weak stationarity:  
  karma.t = karma.transform(y, autolog = autolog, autodiffs = abs(autodiffs), autolags = autolags, r2_criterion = r2_criterion, stdout=F)
  
  #If forced log option is set:
  if( autolog == "force" ){
    karma.t$log = T   #force log parameter to T
  }
  #If autodiffs option is forced to skip differencing:
  if( autodiffs == F ){
    karma.t$diffs = 0   #force integrated term to 0
  }
  #If autodiffs is negative, then force the integrated term to the absolute value of the negative number: 
  if( autodiffs < 0 ){
    karma.t$diffs = abs(autodiffs)
  }
  
  # ----- KARMA model selection (combinatorial optimisation on max F_[mape_out]) -----
  
  if(method == "karma" | method == "greedy-karma"){
    
    # Local optimisation ----- 
    #Create initial solution vectors:
    ar_bool = logical(max_ar)
    ma_bool = logical(max_ma)    
    
    if(method == "karma"){
      
      #Initialisation parameters:
      diffs = karma.t$diffs    #diffs = karma_tmp$model_terms[[2]]
      
    }else if(method == "greedy-karma"){
      
      #Apply greedy Box-Jenkins method as initial solution:
      karma_tmp = karma.boxjenkins( y = y, log = karma.t$log, diffs = karma.t$diffs, xreg = xreg, N=N, box_test = box_test, fixed=fixed, max_ar=max_ar, max_ma=max_ma, max_conv=max_conv, max_iter=max_iter, plot = F, stdout = stdout)  #karma_tmp = auto.boxjen(y, fixed=F)    #<- method = "box-jenkins-karma"
      
      #Initialisation from greedy solution:
      diffs = karma_tmp$model_terms[[2]]
      
      # Local optimisation ----- 
      #Create initial solution vectors:
      ar_bool[ 1:karma_tmp$model_terms[[1]] ] = T      #<- method = "greedy-karma"
      ma_bool[ 1:karma_tmp$model_terms[[3]] ] = T
      
    }
    
    #Create solution vector:
    sol = list(ar = ar_bool, ma = ma_bool) #sol = karma.orderbin(ar_terms = max_ar, ma_terms = max_ma, diffs = diffs, format = 'true-false') #list(ar_terms=c(1,2), ma_terms=c(3,4)) 
    
    #Evaluate initial solution:
    karma_tmp = ModelEvaluationFunction(sol = sol, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = T, test_pct = test_pct, test_type = test_type, plot = plot ) 
    eval = karma_tmp$eval
    
    #Minimise objective function:
    L = karma.hillclimbing( sol = sol, y = y, ac_criterion = ac_criterion, xreg = xreg, metric = metric, cv = cv, diffs = diffs, max_rep = max_rep, optimiser = optimiser, fixed = T, test_pct = test_pct, test_type = test_type, mutations = mutations, plot = plot )
    
    maxeval = as.numeric(L[1])
    solopt = L[[2]]    
    
    karma.opt = ModelEvaluationFunction(sol = solopt, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = T, test_pct = test_pct, test_type = test_type, plot = plot );
    #print( karma.opt$eval ) #<- must be equal to maxeval
    
    karma.opt$solopt = solopt
    return(karma.opt)
  }
  
  # ------
  
  # ----- Box-Jenkins model selection -----------------------  
  if( method == "greedy" ){
    #Apply Box-Jenkins method:
    karma.bj = karma.boxjenkins(y = y, log = karma.t$log, diffs = karma.t$diffs, xreg = xreg, N=N, box_test = box_test, fixed=fixed, max_ar=max_ar, max_ma=max_ma, max_conv=max_conv, max_iter=max_iter, plot = F, stdout = stdout)    
    
    return(karma.bj)
  }
  
}
#Usage: 
#magic.k = auto.boxjen(y)
#karma.cv(magic.k)
#plot(forecast(magic.k))
#or:
#auto.boxjen(y, method = 'karma')

#NB: Add modes: box-jenkins-karma, stochastic neighbourhood function in hill climbing (find random neighbours)


#' Automated ARIMA model selection according to flexible optimisation criteria.  
#'
#' @param sol ARIMA order in the form of a candidate solution to the optimisation problem.
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param ar_terms Autoregressive terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param ma_terms Move average terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param diffs_ Differencing step: Indicates whether the candidate solution needs to be differenced for stationarity (and to what degree); {0,1,...,n}; type <int>.
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection.
#' @param fixed_ Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param cv_ Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation.
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30).
#' @param test_type Train-test split type, i.e. percentage or fixed window; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series.
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' sol = list(ar_terms=c(1,2), ma_terms=c(3,4)) 
#' ModelEvaluationFunction(y, sol, 0)

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

ModelEvaluationFunction <- function(sol, y, diffs_ = 0, metric = "MAPE", fixed_ = T, xreg = NULL, cv_ = 'out', test_pct = 20, test_type = "percentage", plot = T){
  
  #Break down current solution:
  ar_bool = sol$ar
  ma_bool = sol$ma
  
  #Get ARMA terms of current solution:
  ar_terms = which(ar_bool==T)
  ma_terms = which(ma_bool==T)
  
  if(length(ar_terms) == 0){
    ar_terms = 0
  }
  if(length(ma_terms) == 0){
    ma_terms = 0
  }
  
  #Fit fixed model:
  karma0 = karma.fit(y = y, order = list(ar_terms, diffs_, ma_terms), xreg = xreg, fixed = fixed_)  
  
  if( metric == "AIC" ){
    karma0$eval = karma0$aic
  }else if( metric == "BIC" ){
    karma0$eval = karma0$bic
  }else if( metric == "AICc" ){
    karma0$eval = karma0$aicc
  }else{
    karma0$eval = karma.cv(karma0 = karma0, xreg = xreg, metric = metric, cv = cv_, test_pct = test_pct, test_type = test_type, plot = plot)
  }
  
  if( is.na(karma0$eval) ){  #if for any reason no value is assigned to the solution 
    karma0$eval = Inf   #mark solution for rejection
  }
  
  return(karma0)
  
}   

#delme-test: 
#sol = list(ar_terms=c(1,2), ma_terms=c(3,4)) 
#ModelEvaluationFunction(y, sol, 0)



#' Combinatorial optimisation local search algorithm (modified hill-climbing) for ARIMA model selection - part of the model optimisation API.
#'
#' @param sol ARIMA order in the form of a candidate solution to the optimisation problem.
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param ar_terms Autoregressive terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param ma_terms Move average terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param ac_criterion Aucocorrelation / Partial autocorrelation test flag on/off; An optional optimisation constraint which applies portmanteau test on every candidate solution and rejects solutions that do not improve AC/PAC.
#' @param diffs Differencing step: Indicates whether the candidate solution needs to be differenced for stationarity (and to what degree); {0,1,...,n}; type <int>.
#' @param metric Choose a model validation metric that will be used as the main optimisation criterion during model selection.
#' @param fixed Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param cv Choose cross-validation dataset to be used during model selection; "out": Performance of out-of-sample forecast (classic train/test split) will be used for model validation; "in": Performance of in-sample forecast (classic parametric regression type of validation) will be used for model validation.
#' @param test_pct Percentage of train-test split in cross-validation (e.g. 70-30).
#' @param test_type Train-test split type, i.e. percentage or fixed window; "percentage": test_pct = 12 will be read as the 12 of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series.
#' @param optimiser Option on the "neighbourhood function" of the optimisation algorithm; "semi-stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate solutions are chosen deterministically; "stochastic": Once a neighbourhood region (of either AR and MA terms) has been selected randomly, the candidate neighbour solutions are chosen stochastically.
#' @param max_conv For karma.boxjenkins(): Maximum number of iterations without improvement before the algorithm converges forcefully (stuck to a local optimum).
#' @param max_iter For karma.boxjenkins(): Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).
#' @param max_rep For karma-search: Maximum number of iterations without improvement before the algorithm converges naturally (reached a global or local optimum).  
#' @param mutations Optional neighbourhood operator; Mutations flag {T, F}: whether or not to apply random "mutations" (term borrowed from evolutionary algorithms) on a candidate solution when the optimiser is about to converge (a way to escape local optima - works somewhat like an inverse simulated annealing).
#' @return List with optimisation results.
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' kfit = auto.boxjen(y, method = "karma", optimiser = "stochastic", mutations = T, max_rep = 15, max_ar = 4, max_ma = 4)
#' tmp = karma.hillclimbing( sol = kfit$solopt, y = y, fixed=T, diffs = kfit$model_terms[[2]], optimiser = "stochastic")
#' tmp = karma.hillclimbing( sol = tmp[[2]], y = y, fixed=T, diffs = kfit$model_terms[[2]],  optimiser = "stochastic")

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

karma.hillclimbing <- function(sol, y, ac_criterion = F, xreg = NULL, metric = "MAPE", cv = "out", diffs = 0, max_rep = 1, optimiser = "semi-stochastic", fixed = F, test_pct = 20, test_type = "percentage", mutations = F, plot = T ){
  
  sol_local = sol;
  
  i=2;
  c=1;
  mutations_num = 1  #default value (no mutations)
  
  while( c <= max_rep ){    #(i <= length(sol$ar)){
    
    if( mutations == T ){ #escape local optima
      mutations_num = round( (c+1)/2 )   #start applying mutations after local neighbours have been exhausted; at polynomial rate: plot(round( (1:100+1)/2 ), type='l')
      #print(mutations_num)
    }
    
    karma_tmp = ModelEvaluationFunction(sol = sol_local, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = fixed, test_pct = test_pct, test_type = test_type, plot = plot );    
    maxeval_local = karma_tmp$eval
    karma_local = karma_tmp
    
    #print("TEST1")    
    #print(sol_tmp)
    
    #Check neighbourhood of AR:
    if( runif(1) <= 0.5 ){
      
      sol_tmp = sol_local;    
      #sol_tmp$ar[i] = !sol_local$ar[i];
      if(optimiser == "semi-stochastic"){
        sol_tmp$ar[i] = !sol_local$ar[i];  #search local neighbourhood
      }else if(optimiser == "stochastic"){
        if( mutations_num > length(sol_local$ar) ){
          mutations_num_tmp = round(length(sol_local$ar)/2)
        }else{
          mutations_num_tmp = mutations_num
        }
        i_tmp = sample(length(sol_local$ar), mutations_num_tmp)
        sol_tmp$ar[i_tmp] = !sol_local$ar[ i_tmp ];   #search random neighbourhood
      }
      karma_tmp = ModelEvaluationFunction(sol = sol_tmp, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = fixed, test_pct = test_pct, test_type = test_type, plot = plot );
      maxeval_tmp = karma_tmp$eval   
      
      if (maxeval_tmp < maxeval_local){
        if( ac_criterion == T & optimiser == "stochastic" ){  #check if number of autocorrelated lags has decreased
          autocorr_tmp = karma.portmanteau(karma_tmp$residuals)
          autocorr_local = karma.portmanteau(karma_local$residuals)
          if( length(autocorr_tmp$ac_lags) + length(autocorr_tmp$pac_lags) >= length(autocorr_local$ac_lags) + length(autocorr_local$pac_lags) ){
            print("Rejecting neighbour solution due to higher autocorrelation.")
            next;
          }
        }
        #print("TEST2")
        maxeval_local = maxeval_tmp;
        sol_local$ar = sol_tmp$ar;
        print("AR Solution improved: ^")
        cat("AR terms:", which(sol_local$ar), ";", "MA terms:", which(sol_local$ma), "\n") #print(sol_local)
        cat(metric, ": ", maxeval_local, "\n\n") #print(maxeval_local)
        #fprintf('Solution improved ^...(%.4f, %.4f, %.4f, %.4f)\n',sol_local(1),sol_local(2),sol_local(3),sol_local(4));
        #fprintf('Objective value: %f\n', maxeval_local);             
        i=1;
        next;
      }else{
        #sol_tmp = sol_local;    
        #sol_tmp$ar[i] = !sol_local$ar[i];
        if(optimiser == "semi-stochastic"){
          sol_tmp$ar[i] = !sol_local$ar[i];
        }else if(optimiser == "stochastic"){
          #i_tmp = sample(length(sol_local$ar), c)
          sol_tmp$ar[i_tmp] = !sol_local$ar[ i_tmp ];
        }
        karma_tmp = ModelEvaluationFunction(sol = sol_tmp, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = fixed, test_pct = test_pct, test_type = test_type, plot = plot );
        maxeval_tmp = karma_tmp$eval
        i=i+1;
      }
    }
    
    #Check neighbourhood of MA:
    else if( runif(1) > 0.5 ){
      
      sol_tmp = sol_local;    
      if(optimiser == "semi-stochastic"){
        sol_tmp$ma[i] = !sol_local$ma[i];  #search local neighbourhood
      }else if(optimiser == "stochastic"){
        if( mutations_num > length(sol_local$ma) ){
          mutations_num_tmp = round(length(sol_local$ma)/2)
        }else{
          mutations_num_tmp = mutations_num
        }        
        i_tmp = sample(length(sol_local$ma), mutations_num_tmp)
        sol_tmp$ma[i_tmp] = !sol_local$ma[ i_tmp ];   #search random neighbourhood
      }
      karma_tmp = ModelEvaluationFunction(sol = sol_tmp, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = fixed, test_pct = test_pct, test_type = test_type, plot = plot );
      maxeval_tmp = karma_tmp$eval         
      
      if (maxeval_tmp < maxeval_local){
        if( ac_criterion == T & optimiser == "stochastic" ){  #check if number of autocorrelated lags has decreased
          autocorr_tmp = karma.portmanteau(karma_tmp$residuals)
          autocorr_local = karma.portmanteau(karma_local$residuals)
          if( length(autocorr_tmp$ac_lags) + length(autocorr_tmp$pac_lags) >= length(autocorr_local$ac_lags) + length(autocorr_local$pac_lags) ){
            print("Rejecting neighbour solution due to higher autocorrelation.")
            next;
          }
        }        
        #print("TEST2")
        maxeval_local = maxeval_tmp;
        sol_local$ma = sol_tmp$ma;
        print("MA Solution improved: ^")
        cat("AR terms:", which(sol_local$ar), ";", "MA terms:", which(sol_local$ma), "\n") #print(sol_local)
        cat(metric, ": ", maxeval_local, "\n\n") #print(maxeval_local)
        #fprintf('Solution improved ^...(%.4f, %.4f, %.4f, %.4f)\n',sol_local(1),sol_local(2),sol_local(3),sol_local(4));
        #fprintf('Objective value: %f\n', maxeval_local);             
        i=1;
        next;
      }else{
        #sol_tmp = sol_local;    
        #sol_tmp$ma[i] = !sol_local$ma[i];   
        if(optimiser == "semi-stochastic"){
          sol_tmp$ma[i] = !sol_local$ma[i];
        }else if(optimiser == "stochastic"){
          sol_tmp$ma[i_tmp] = !sol_local$ma[ i_tmp ];
        }        
        karma_tmp = ModelEvaluationFunction(sol = sol_tmp, y = y, metric = metric, cv = cv, xreg = xreg, diffs_ = diffs, fixed_ = fixed, test_pct = test_pct, test_type = test_type, plot = plot );
        maxeval_tmp = karma_tmp$eval
        i=i+1;
      }
    }
    
    #i=i+1;
    
    
    if( i >= length(sol$ar) ){
      c = c+1;
    }
  }  
  return(list(maxeval_local, sol_local))
}



#' Convert fixed order terms to 0/NA mask (needed for function Arima()). 
#'
#' @param ar_terms Autoregressive terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param ma_terms Move average terms; can be in the form of a vector (fixed_ = FALSE) or a list of vectors (fixed_ = TRUE).
#' @param diffs Differencing step: Indicates whether the candidate solution needs to be differenced for stationarity (and to what degree); {0,1,...,n}; type <int>.
#' @param format Type of transformation; "na-zero": converts numeric fixed order to 0/NA mask; "boolean": converts numeric fixed order to T/F.
#' @return Converted vector.
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' karma.orderbin(ar_terms = 15, ma_terms = 15, diffs = 1, format = 'true-false')

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

karma.orderbin <- function(ar_terms, ma_terms, diffs = 0, format = 'na-zero' ){
  
  if( length(ar_terms) == 0 ){
    ar_terms = 0
  }
  if( length(ma_terms) == 0 ){
    ma_terms = 0
  }
  
  fixed_terms = c( numeric(max(ar_terms)), numeric(max(ma_terms)) )  #initialise zero array
  
  #Add fixed ARMA terms:
  fixed_terms[c(ar_terms, max(ar_terms) + ma_terms)] = NA   #fix constant terms (set to NA)
  
  #Add constant term only if diffs are 0:
  if(is.na(diffs)){
    diffs = 0
  }
  if(diffs <= 0){
    fixed_terms = c( fixed_terms, NA ) #add constant term
  }
  
  #Output in 2 forms:
  if( format == 'na-zero' ){
    return(fixed_terms)
  }else{
    fixed_terms[is.na(fixed_terms)] = T
    fixed_terms[fixed_terms == 0] = F
    return(fixed_terms)
  }
}




