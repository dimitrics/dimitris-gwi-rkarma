#' #Forecast and plot future periods: Combine multiple weak learners (fitted models) from a pre-trained "karma.ensemble" object to produce future time periods.
#'
#' @param kensemble Object of class "karma.ensemble"; ensemble of unique trained models. <karma.ensemble>
#' @param h Number of periods to forecast ahead. <numeric> integer
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.forecast".
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' # Create and forecast ensemble of 10 weak learners:
#' karma.forecast( karma.ensemble(JohnsonJohnson) )

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


karma.forecast <- function(kensemble, h = 12, plot = T, xreg = NULL, stdout = F){
  
  kforecast = list()
  #kforecast$objects = list()
  
  y = kensemble$x
  nsamples = length(kensemble$kfit_list)
  mean_pred_mat = matrix(0, nsamples, h)   #matrix with the average out-of-sample predictions (y-pred) of all models
  lower_ci_mat = matrix(0, nsamples, h)    #matrix with the lower prediction interval of all models
  upper_ci_mat = matrix(0, nsamples, h)
  
  for(i in 1:nsamples){
    kforecast[[i]] = forecast( kensemble$kfit_list[[i]], h, xreg = xreg ) 
    mean_pred_mat[i,] = kforecast[[i]]$mean
    lower_ci_mat[i, ] = kforecast[[i]]$lower[,1]   #[,1]: 80% c.i., [,2]: 95% c.i.
    upper_ci_mat[i, ] = kforecast[[i]]$upper[,1]      
  }
  
  if(sum(kensemble$uniqueness_flag_vec) > 1){
    kforecast$agg_unseen_forecast = apply(mean_pred_mat[kensemble$uniqueness_flag_vec, ], 2, mean)
    kforecast$agg_lower_ci = apply(lower_ci_mat[kensemble$uniqueness_flag_vec, ], 2, mean)
    kforecast$agg_upper_ci = apply(upper_ci_mat[kensemble$uniqueness_flag_vec, ], 2, mean)
  }else if(sum(kensemble$uniqueness_flag_vec) == 1){
    kforecast$agg_unseen_forecast = mean_pred_mat[kensemble$uniqueness_flag_vec, ]
    kforecast$agg_lower_ci = lower_ci_mat[kensemble$uniqueness_flag_vec, ]
    kforecast$agg_upper_ci = upper_ci_mat[kensemble$uniqueness_flag_vec, ]
  }else{
    stop("karma.forecast(): No acceptable models were fitted - exiting...")
  }
  
  if(plot == T){
    #Plot forecast with C.I.:
    #plot(kforecast$agg_unseen_forecast, type="l", ylim=c(min(kforecast$agg_lower_ci), max(kforecast$agg_upper_ci)));
    #points(kforecast$agg_lower_ci, col='red', type="l", pch=22, lty=2); 
    #points(kforecast$agg_upper_ci, col='red', type="l", pch=22, lty=2); 
    
    title_ = paste("Boosted forecast","\n(", h, "periods ahead )")
    plot(c(y, kforecast$agg_unseen_forecast), type="l", xlab = 'time', ylab = 'y', main = title_);
    points((length(y) + 1):(length(y) + h), kforecast$agg_unseen_forecast, col = "blue", type = "l", pch=22, lwd=2);
    points((length(y) + 1):(length(y) + h), kforecast$agg_lower_ci, col='green', type="l", pch=22, lty=2);
    points((length(y) + 1):(length(y) + h), kforecast$agg_upper_ci, col='green', type="l", pch=22, lty=2);
  }
  
  class(kforecast) = "karma.forecast"
  return(kforecast)
}

