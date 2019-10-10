#' Convert "ARIMA" object to "karma.fit" object.
#'
#' @param arma0 "ARIMA" object.
#' @return Object of class "karma.fit"; (extends class "ARIMA" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#' karma.cv( karma.cast(auto.arima(WWWusage)) )

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


#Typecast 'forecast' object to 'karma':
karma.cast <- function(arma0){

  #Get AR/MA terms:
  terms = c()
  terms$ar = c()
  terms$ma = c()
  terms$sar = c()
  terms$sma = c()  
  
  for(i in 1:length(arma0$coef)){
    term_str = names(arma0$coef)[i]
    
    if( str_sub(term_str, 1, 1) != 's' ){
      term_name = str_sub(term_str, 1, 2)
      term_value = str_sub(term_str, 3, -1)
      if(term_name=='ar'){
        terms$ar = c(terms$ar, as.numeric(term_value))
      }else if(term_name=='ma'){
        terms$ma = c(terms$ma, as.numeric(term_value))
      }
    }else if( str_sub(term_str, 1, 1) == 's' ){
      term_name = str_sub(term_str, 1, 3)
      term_value = str_sub(term_str, 4, -1)      
      if(term_name=='sar'){
        terms$sar = c(terms$sar, as.numeric(term_value))
      }else if(term_name=='sma'){
        terms$sma = c(terms$sma, as.numeric(term_value))
      }
    }
  }
  
  if( length(terms$ar) == 0){
    terms$ar = 0
  }
  if( length(terms$ma) == 0){
    terms$ma = 0
  }
  if( length(terms$sar) == 0){
    terms$sar = 0
  }
  if( length(terms$sma) == 0){
    terms$sma = 0
  }    
  # if( length(arma0$mask) > (max(terms$ar)+max(terms$ma)+max(terms$sar)+max(terms$sma)) ){ #check if there is a constant term in the model
  #   diffs = 0 #constant term means no differencing
  # }else{
  #   diffs = 1
  # }
  
  arma0$model_terms = c(max(terms$ar), arma0$arma[length(arma0$arma)-1], max(terms$ma))
  arma0$seasonal_terms = c(max(terms$sar), arma0$arma[length(arma0$arma)], max(terms$sma))
  return(arma0)
}
#Usage: karma.cv( karma.cast(auto.arima(y)) )  #karma.cv( magic.karma(y) )
#karma.cast( Arima( y = WWWusage, order = c(1,1,0)) )$model_terms

