#' Fit ARIMA model to univariate time-series. Mostly a wrapper to 'Arima' function in the 'forecast' package; the difference is that it facilitates fitting fixed term models. 
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param order Specification of the non-seasonal part of the ARIMA model (see documentation for Arima(), argument 'order').
#' @param log Logarithmic transformation flag; Indicates whether the input series needs to be log-transformed for stationarity; {T, F}; type <logical>.
#' @param fixed Fixed term flag. Indicate whether the fixed term option in Arima() needs to be switched on during model selection; {T, F}; type <logical>. 
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "ARIMA" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Fit and summarise model:
#' summary(karma.fit(WWWusage, order=list(c(1,2), 1, c(3,4)), fixed=T))
#' # Using forecast() from package 'forecast':
#' plot(forecast( karma.fit(WWWusage, order=c(2,1,3), fixed=F) )) 

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

#----- Fit ARMA model --------

karma.fit <- function(y, order = c(0, 0, 0), seasonal = c(0, 0, 0), log = F, xreg = NULL, fixed = F, stdout = F){
  
  #Init:
  karma_terms = order
  fixed_terms = NA
  ar_terms = 0
  ma_terms = 0
  if( log == T ){
    y = log(y)
    #lambda = !log
  }
  
  #Fit:
  carma0 = tryCatch({
    if( fixed == F ){
      if(is.numeric(karma_terms) & length(karma_terms) == 3){
        carma0 = Arima(y, order = karma_terms, seasonal = seasonal, xreg = xreg, method = 'CSS-ML')
      }else{
        stop("For fixed = F, parameter karma_terms has to be a vector of type 'numeric' of length equal to 3.")
      }
    }else if ( fixed == T ){
      if( class(karma_terms) == "list" ){
        ar_terms = karma_terms[[1]]
        ma_terms = karma_terms[[3]]
        if(length(ar_terms) == 0){
          ar_terms = 0
        }
        if(length(ma_terms) == 0){
          ma_terms = 0
        }
      }
      else{
        stop("For fixed = T, parameter karma_terms has to be a list of vectors of type 'numeric'.")
      }
      
      #Create fixed-order vector:
      fixed_terms = karma.orderbin( ar_terms, ma_terms, diffs = karma_terms[[2]] )
      
      #Fit fixed model:
      carma0 = Arima(y, order = c(max(ar_terms), karma_terms[[2]], max(ma_terms)), seasonal = seasonal, xreg = xreg, fixed = fixed_terms, transform.pars = F, method = 'CSS-ML')
    }
  }, error=function(e){
    if( stdout == T ){
      cat("WARNING:",conditionMessage(e), "\n")
      msg0 = "WARNING: Arima() with method='CSS-ML' failed; trying CSS estimation."
      print(msg0)
    }
    #warning_msg = c(warning_msg, msg0)
    if( fixed == F ){
      if(is.numeric(karma_terms) & length(karma_terms) == 3){
        carma0 = tryCatch({
          carma0 = Arima(y, order = karma_terms, seasonal = seasonal, xreg = xreg, method = "CSS")  #Try again fitting using the more robust conditional sum of squares likelihood
        }, error=function(e1){
          if( stdout == T ){
            cat("WARNING:",conditionMessage(e1), "\n")
            msg1 = "WARNING: Arima() with method='CSS' failed; trying ML estimation."
            print(msg1)
          }
          carma0 = tryCatch({
            carma0 = Arima(y, order = karma_terms, seasonal = seasonal, xreg = xreg, method = "ML")
          }, error=function(e2){
            cat("WARNING:",conditionMessage(e2), "\n")
            cat("\n-------- ABORTING karma.fit: Arima() failed to fit model with ML-CSS, CSS, and ML methods.-------\n")
            cat("Order:", karma_terms, "Seasonal:", seasonal, "\n\n")
            return(0)
          })
        })
      }
      else{
        stop("For fixed = F, parameter karma_terms has to be a vector of type 'numeric' of length equal to 3.")
      }
    }else if ( fixed == T){
      if( class(karma_terms) == "list" ){
        ar_terms = karma_terms[[1]]
        ma_terms = karma_terms[[3]]
        if(length(ar_terms) == 0){
          ar_terms = 0
        }
        if(length(ma_terms) == 0){
          ma_terms = 0
        }        
      }
      else{
        stop("For fixed = T, parameter karma_terms has to be a list of vectors of type 'numeric'.")
      }      
      
      #Create fixed-order vector:
      fixed_terms = karma.orderbin( ar_terms, ma_terms, diffs = karma_terms[[2]] )
      
      #Fit fixed model:
      carma0 = Arima(y, order = c(max(ar_terms), karma_terms[[2]], max(ma_terms)), seasonal = seasonal, xreg = xreg, method = 'CSS', fixed = fixed_terms, transform.pars=F)
    }
    
    # carma0$fixed_flag = fixed
    # carma0$fixed_terms = fixed_terms
    # carma0$model_terms = karma_terms
    # class(carma0) = c("karma", class(carma0))
    # return(carma0)
  })
  
  carma0$R2_approx = suppressWarnings( cor(fitted(carma0), carma0$x)^2 )
  carma0$pval_coef = suppressWarnings( round( (1-pnorm(abs(carma0$coef)/sqrt(diag(carma0$var.coef))))*2, 3 ) )
  carma0$xreg = xreg
  carma0$log = log
  carma0$fixed_flag = fixed
  carma0$fixed_terms = fixed_terms
  carma0$model_terms = karma_terms
  carma0$model_seasonal = seasonal
  carma0$model_order = c(max(ar_terms), karma_terms[[2]], max(ma_terms))
  carma0$ar_terms = ar_terms
  carma0$ma_terms = ma_terms
  class(carma0) = c("karma.fit", class(carma0))
  return(carma0)  
  
} #test: plot(forecast( karma.fit(y, karma_terms=c(2,1,3), fixed=F) ))  and: summary(karma.fit(y, karma_terms=list(c(1,2), 1, c(3,4)), fixed=T))



