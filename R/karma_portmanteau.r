#' Portmanteau test to detect statistically significant autocorrelation and partial autocorrelation lags.   
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param N Maximum lag at which to calculate autocorrelation and partial autocorrelatin functions; see documentation for acf(), pacf(). 
#' @return Object of class "karma.portmanteau".
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}, \code{\link{car}}
#' @export
#' @examples
#' autocorrelation.lags = karma.portmanteau( magic.karma(WWWusage)$residuals )

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


karma.portmanteau <- function(y, plot = F, N = 100){
  
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
  
  ac_test = list()
  ac_test$ac_lags = ac_lags
  ac_test$pac_lags = pac_lags
  ac_test$box_test_pval = Box.test(x = y, type = "Box-Pierce", lag = 2)$p.value   #Ljung-Box test of overall autocorrelation
  class(ac_test) = "karma.portmanteau"
  
  return(ac_test)
} #usage: karma.portmanteau(y, plot = T)

