#' Model validation metrics.  
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param y_hat Predicted values of univariate time series.
#' @param metric A model validation metric (MAPE, MAE, MSE, RMSE, R2).
#' @return Output value of selected metric; <numeric>.
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' kfit = magic.karma(WWWusage); karma.validate(fitted(kfit), WWWusage, metric = "MAPE")

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------

karma.validate <- function(y, y_hat, metric = "MAPE"){
  y = as.numeric(y)
  y_hat = as.numeric(y_hat)
  
  if(metric == "MAPE"){
    eval = sum(abs((y - y_hat)/y))*100 / length(y_hat) 
  }else if(metric == "MAE"){
    eval = sum(abs(y - y_hat)) / length(y_hat) 
  }else if(metric == "MSE"){
    eval = sum((y_hat - y)^2) / length(y_hat)
  }else if(metric == "RMSE"){
    eval = sqrt( sum((y_hat - y)^2) / length(y_hat) )
  }else if(metric == "R2"){
    eval = cor(y_hat, y)^2
  }else{
    stop('Unknown validation metric.')
  }
  # if(metric == "MASE"){
  #   T = length(y); 
  #   eval = sum(abs(residuals(kfit))) / ( (T/(T-1))*sum(abs(diff(y)))  )
  # }
  
  return(eval)
} #usage: kfit = magic.karma(y); karma.validate(fitted(kfit), y, metric = "MAPE")  #<- compare to: mmetric(fitted(kfit), y, "MAPE")

