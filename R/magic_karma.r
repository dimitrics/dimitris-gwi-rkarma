#' Train multitude of models on a univariate time series 
#'
#' @param y A univariate time-series vector; type <numeric> or <ts>.
#' @param model_list List of indexes of training models/algorithms in that order: ms-arima, auto-arima, nms-arima, ld-arima, rw-arima, sbj-arima, smbj-arima, nnetar, ets, bats, tbats; <list>
#' @param stacking Whether to use ensemble learning algorithms or not; <T/F>
#' @param test_pct_train Percentage of train-test split in model training (e.g. 70-30), after the model is trained.
#' @param test_type_train Train-test split type for training, i.e. percentage or fixed window; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series.
#' @param test_pct_valid Percentage of train-test split in model validation (e.g. 70-30), after the model is trained.
#' @param test_type_valid Train-test split type for validation, i.e. percentage or fixed window; "percentage": test_pct = 12 will be read as the 12 percent of the length of the series; "window": test_pct = 12 will be read as the 12 last time points (e.g. months) of the series.
#' @param xreg Optional vector or matrix of exogenous regressors; see documentation for Arima(), package 'forecast'.
#' @param plot Option to depict plots during local search; if TRUE (default), AC and PAC plots are active. <logical>
#' @param stdout Option to output optimisation diagnostics during local search; <logical>
#' @return Object of class "karma.fit"; (extends class "Arima" from package 'forecast').
#' @seealso \code{\link{tseries}}, \code{\link{forecast}}
#' @export
#' @examples
#'  kmodels = magic.karma(JohnsonJohnson)
#'  kmodels[[1]]$fit_obj$aicc
#'  kmodels[[1]]$cv_obj$mape_in
#'  kmodels[[1]]$cv_obj$mape_out  

# Author: Dimitris Tziotis (2016)
# email: dimitrios.tziotis@gmail.com
#----------------------------------------


magic.karma <- function(y, model_list = 1:11, stacking = F, test_pct_train = "auto", test_type_train = "auto",  test_pct_valid = "auto", test_type_valid = "auto", xreg = NULL, plot = F, stdout = F){
  
  model_number = length(model_list)
  
  if( stacking == F ){
    model_number = model_number
  }else{
    model_number = model_number + 1
  }
  
  #------------------------------------------------------------------------------
  #------------- Train models ----------------------------------  
  #------------------------------------------------------------------------------
  
  model_vec = rep(list(NULL), model_number)
  c = 1
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 1){
      cat("\nFitting Markov Selection model (MS-ARIMA) with auto.karma...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "ms_arima"
      
      model_vec[[c]]$fit_obj = tryCatch({
        auto.karma(y, method = "markov-selection", test_pct = test_pct_train, test_type = test_type_train, xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){
          #model_vec[[c]]$model_specs = paste("Unfitted", model_vec[[c]]$function_name , sep = "_")
        })
      }
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 2){
      cat("\nFitting Hyndman's ARIMA model with auto.arima...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "auto_arima"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = auto.arima(y, xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 3){
      cat("\nFitting Nonparametric Markov Selection model (NMS-ARIMA) with auto.karma...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "nms_arima"
      
      model_vec[[c]]$fit_obj = tryCatch({
        auto.karma(y, method = "markov-selection", test_pct = test_pct_train, test_type = test_type_train, metric = "MAPE", xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){
          #model_vec[[c]]$model_specs = paste("Unfitted", model_vec[[c]]$function_name , sep = "_")
        })
      }
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 4){
      cat("\nFitting Local Descent model (LD-ARIMA) with auto.karma...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "ld_arima"
      
      model_vec[[c]]$fit_obj = tryCatch({
        auto.karma(y, method = "local-descent", test_pct = test_pct_train, test_type = test_type_train, xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){
          #model_vec[[c]]$model_specs = paste("Unfitted", model_vec[[c]]$function_name , sep = "_")
        })
      }
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 5){
      cat("\nFitting Random Walk model (RW-ARIMA) with auto.karma...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "rw_arima"
      
      model_vec[[c]]$fit_obj = tryCatch({
        auto.karma(y, method = "random-walk", test_pct = test_pct_train, test_type = test_type_train, xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){
          #model_vec[[c]]$model_specs = paste("Unfitted", model_vec[[c]]$function_name , sep = "_")
        })
      }
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 6){
      cat("\nFitting Stochastic Box Jenkins model (SBJ-ARIMA) with auto.boxjen...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "auto_boxjen"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = auto.boxjen(y, test_pct = test_pct_train, test_type = test_type_train, xreg = xreg)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 7){
      cat("\nFitting Stochastic Multiseasonal Box Jenkins model (SMBJ-ARIMA) with auto.boxjen...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "auto_boxjen_multiseasonal"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = auto.boxjen(y, test_pct = test_pct_train, test_type = test_type_train, xreg = xreg, fixed = T)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1  
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 8){
      cat("\nFitting nnetar...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "nnetar"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = nnetar(y)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1    
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 9){
      cat("\nFitting ets...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "ets"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = ets(y)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1  
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 10){
      cat("\nFitting bats...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "bats"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = bats(y)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1  
    }
  }
  
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 11){
      cat("\nFitting tbats...\n")
      model_vec[[c]]$fit_obj = NA
      model_vec[[c]]$cv_obj = NA
      model_vec[[c]]$error_msg_fit = NA
      model_vec[[c]]$error_msg_cv = NA
      model_vec[[c]]$function_name = "tbats"
      
      model_vec[[c]]$fit_obj = tryCatch({
        model_vec[[c]]$fit_obj = tbats(y)
      }, error = function(e2) {
        model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      if(!is.logical(model_vec[[c]]$fit_obj)){
        model_vec[[c]]$model_specs = tryCatch({
          model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
        }, error = function(e2){})
      }  
      c = c+1  
    }
  }
  
  
  if( stacking==T ){
    
    cat("\nTraining Random Walk learners with Stacked Stochastic Forecast (SSF)...\n")
    model_vec[[c]]$fit_obj = NA
    model_vec[[c]]$cv_obj = NA
    model_vec[[c]]$error_msg_fit = NA
    model_vec[[c]]$error_msg_cv = NA
    model_vec[[c]]$function_name = "ms_arima"
    
    model_vec[[c]]$fit_obj = tryCatch({
      karma.ensemble(y, family = "sarima", test_pct = test_pct_train, test_type = test_type_train, xreg = xreg)
    }, error = function(e2) {
      model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
      handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
    })
    if(!is.logical(model_vec[[c]]$fit_obj)){
      model_vec[[c]]$model_specs = tryCatch({
        model_vec[[c]]$model_specs = karma.specs(model_vec[[c]]$fit_obj)
      }, error = function(e2){
        #model_vec[[c]]$model_specs = paste("Unfitted", model_vec[[c]]$function_name , sep = "_")
      })
    }
    c = c+1
    
  }
  
  
  # cat("\nFitting mlp...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "nnfor::mlp"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj = nnfor::mlp(y)
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nFitting tsml.rf...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "tsml.rf"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj = tsml.rf(y, ar_terms = frequency(y))
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nFitting prophet.fit...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "prophet.fit"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj = prophet.fit(y)
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nFitting HoltWinters multiplicative...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "HoltWinters-mult"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj =  HoltWinters(y, seasonal = "mult")
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nFitting HoltWinters additive...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "HoltWinters-add"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj =  HoltWinters(y, seasonal = "add")
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nFitting bsts...\n")
  # model_vec[[c]]$fit_obj = NA
  # model_vec[[c]]$cv_obj = NA
  # model_vec[[c]]$error_msg_fit = NA
  # model_vec[[c]]$error_msg_cv = NA
  # model_vec[[c]]$function_name = "bsts"
  # 
  # model_vec[[c]]$fit_obj = tryCatch({
  #   model_vec[[c]]$fit_obj = bsts.fit(y)
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_fit <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  
  #------------------------------------------------------------------------------------------------------------------
  #------------- Validate models ----------------------------------  
  #------------------------------------------------------------------------------------------------------------------
  
  c = 1
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 1){
      cat("\nValidating Markov Selection model (MS-ARIMA) from auto.karma...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 2){  
      cat("\nValidating Hyndman's ARIMA model from auto.arima...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1  
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 3){  
      cat("\nValidating Nonparmetric Markov Selection model (NMS-ARIMA) from auto.karma...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 4){
      cat("\nValidating Local Descent model (LD-ARIMA) from auto.karma...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 5){
      cat("\nValidating Random Walk model (RW-ARIMA) from auto.karma...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 6){
      cat("\nValidating Stochastic Box Jenkins model (SBJ-ARIMA) from auto.boxjen...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 7){
      cat("\nValidating Stochastic Multiseasonal Box Jenkins model (SMBJ-ARIMA) from auto.boxjen...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, plot = plot, test_pct = test_pct_valid, test_type = test_type_valid, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1  
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 8){
      cat("\nValidating nnetar...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, test_pct = test_pct_valid, test_type = test_type_valid, plot = plot, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1    
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 9){
      cat("\nValidating ets...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, test_pct = test_pct_valid, test_type = test_type_valid, plot = plot, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1  
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 10){
      cat("\nValidating bats...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, test_pct = test_pct_valid, test_type = test_type_valid, plot = plot, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1  
    }
  }
  
  if(!is.na(model_list[c])){
    if(model_list[c] == 11){   
      cat("\nValidating tbats...\n")
      model_vec[[c]]$cv_obj = tryCatch({
        model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, test_pct = test_pct_valid, test_type = test_type_valid, plot = plot, ensemble = T);
      }, error = function(e2) {
        model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
        handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
      })
      c = c+1  
    }
  }
  
  
  if( stacking == T ){
    
    cat("\nValidating Stacked Stochastic Forecast (SSF) from karma.ensemble...\n")
    model_vec[[c]]$cv_obj = tryCatch({
      model_vec[[c]]$cv_obj = karma.cv(model_vec[[c]]$fit_obj, test_pct = test_pct_valid, test_type = test_type_valid, plot = plot);
    }, error = function(e2) {
      model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
      handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
    })
    c = c+1  
    
  }
  
  # cat("\nValidating mlp...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj = mlp.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nValidating tsml.rf...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj = tsml.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nValidating prophet.fit...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj = prophet.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nValidating HoltWinters multiplicative...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj =  hw.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nValidating HoltWinters additive...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj =  hw.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  # 
  # cat("\nValidating bsts...\n")
  # model_vec[[c]]$cv_obj = tryCatch({
  #   model_vec[[c]]$cv_obj = bsts.cv(model_vec[[c]]$fit_obj, test_pct = test_pct, test_type = test_type, plot = plot);
  # }, error = function(e2) {
  #   model_vec[[c]]$error_msg_cv <<- conditionMessage(e2)
  #   handle_error(func_name = model_vec[[c]]$function_name, e2 = e2, stdout = stdout)
  # })
  # c = c+1  
  
  
  class(model_vec) = "magic.karma"
  return(model_vec)
} #usage: models = magic.karma(JohnsonJohnson); models[[1]]$fit_obj$aicc; models[[1]]$cv_obj$mape_in; models[[1]]$cv_obj$mape_out



#---- error handler ---------------------

handle_error <- function(func_name, e2, stdout = F){
  if( stdout == T ){
    cat("\n-------- WARNING:", func_name, "exited with error: -------\n")
    cat(conditionMessage(e2), "\n")
  }
  return(NA)
}


#--- Return model specifics from fit object ---------------

karma.specs <- function(kfit){
  if( class(kfit)[1] == "nnetar" ){
    layers = kfit$model[[1]][[1]]
    model_desc = paste("Shallow autoregressive neural network with ", layers[2], " hidden nodes", sep = "")
  }else{
    model_desc = as.character(kfit)
  }
  if(length(model_desc) > 1){
    model_desc = "Unfitted_BoxJenkins"
  }
  return(model_desc)
}

#--------------------------------------------------------------------------------------------------------------


#models = magic.karma(JohnsonJohnson)
#names(models[[1]])

#TODO:
# make karma.cv compatible with all *.cv wrappers
# make karma.forecast compatible with forecast() and all *.forecast() wrappers
# make all non-forecast cv functions in karma_extras return cv object
# make size_cv() usable in karma_extras

# y = c(-1, 0, 40, 0, 0, 0, 40, 40, 0, 40, 0, 0, 0, 0, 40, 0, 0, 40, 0, 40, 0, 80, 0, 40, 0, 0, 40, 40, 0, 0, 40, 0, 0, 0, 22, 0, 40, 0, 0, 0, 0, 40, 0, 19, 40, -3, 40, 40, 0, 0, 80, 0, 0, 40, 0, 40, 0, 0, 40, 40, 0, 0, 40, 0, 0, 0, 0, 0, 40, 0, 0, 40, 40, 0, 80, 0, 40, 0, 80, -3, 0, 40, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 80, 0, 0, 0, 0, 0, 0, 40, 40, 0, 0, 0, 20, 0, 40, 40)
# y = ts(y, start = c(2008,1), end = c(2012,7), frequency = 26)
# models = magic.karma(y)
# 
# names(models[[1]])



# 
# magic.karma2 <- function(y, stdout = F, model_number = 20, algo_train_func = NA){
#   
#   if(is.NA(algo_train_func)){
#     algo_train_func = c(auto.karma, auto.arima, auto.boxjen, nnetar, ets, bats, tbats, nnfor::mlp, tsml.rf, prophet.fit, HoltWinters, bsts.fit )
#   }
#   #algo_cv_func = c(auto.karma, auto.arima, auto.boxjen, nnetar, ets, bats, tbats, nnfor::mlp, tsml.rf, prophet.fit, HoltWinters, bsts.fit )
#   
#   
#   model_number = 20#length(algo_train_func)
#   
#   #------------- Train ----------------------------------  
#   
#   model_vec = rep(list(NA), model_number)
#   
#   for(i in 1:length(algo_train_func)){
#     
#     cat("\nFitting SARIMA(X) model with auto.karma...\n")
#     model_vec[[i]]$fit_obj = tryCatch({
#       model_vec[[i]]$fit_obj = algo_train_func[[i]](y)
#     }, error = function(e2) {
#       model_vec[[i]]$error_msg_fit = conditionMessage(e2)
#       if( stdout == T ){
#         cat("\n-------- WARNING: auto.karma exited with error: -------\n")
#         cat(conditionMessage(e2), "\n")
#       }
#     })
#     model_vec[[i]]$function_name = algo_train_func[[i]]
#     
#   }
#   
# }
# magic.karma2(JohnsonJohnson)



