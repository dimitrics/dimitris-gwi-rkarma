# karma
<i> A library for automated time series prediction</i> 

<b>Statistical forecasting with automated model selection and hyperparameter optimization</b>
- One-line ARIMA selection and fitting with both stochastic and deterministic training options (meant as an alternative to auto.arima)
- One-line Box-Jenkins model constuction with training options for approximation metaheuristic

<b>Manual time series analysis</b>
- Backwards compatibility with package 'forecast' (for Hyndman's ARIMA, neural network autoregression, and exponential smoothing)
- In-sample and out-of-sample model validation (diagnose overfitting/underfitting via train-test splitting)

<b>AutoML for univariate series</b>
- Train multiple weak learners to make stronger predictions (ensemble learning / stacking)
- Train all models in 'karma' and 'forecast' on a single time series dataset
- Train all models in 'karma' and 'forecast' on multiple time series datasets and benchmark their results


<br/>


*__Model selection algorithms__*

- MS-ARIMA (Markov Selection ARIMA)
- RW-ARIMA (Random Walk ARIMA)
- LD-ARIMA (Local Descent ARIMA)
- Nonparametric Markov Selection ARIMA (NMS-ARIMA)
- Stochastic Box-Jenkins (SBJ-ARIMA)
- Stochastic Multiseasonal Box-Jenkins (SMBJ-ARIMA)
- NP-Hard ARIMA (NPH-ARIMA) 
- Stacked Stochastic Forecast (SSF)

<br />

*__Instructions__*

<br />

__You first install these packages:__

install.packages(c('devtools', 'forecast', 'stringr'))

<br />

__Then you load devtools:__

library(devtools)

<br />

__Then you install this package:__

install_github("snarf-snarf/karma")

<br />

__Then you load all packages:__
```R
library(karma)

library(forecast)

library(stringr)
```

<br />


__Examples:__


Search for optimal SARIMA model using auto.karma with cross-validation:

```R
sfit <- auto.karma( mdeaths )   # train SARIMA/ARIMA model
```

Retrain selected model (training set) and validate prediction accuracy on held out data (test set):

```R
karma.cv( sfit )   # will plot predicted vs. actual test-set data and print out test MAPE
```

Forecast future periods:

```R
plot(forecast( sfit ))    # plot projection
```

Use karma API with 'forecast' package:

```R
afit <- auto.arima(mdeaths)   # auto.arima
karma.cv( afit )
plot(forecast( afit ))    
nfit <- nnetar(mdeaths)       # shallow neural network autoregression
karma.cv( nfit )
plot(forecast( nfit ))   
efit <- ets(mdeaths)          # exponential smoothing state space model
karma.cv( efit )
plot(forecast( efit ))    
                              # ...etc...
```


Train multiple models on a single series:

```R
kmodels <- magic.karma(mdeaths)    # train all models
karma.cv(kmodels[[1]]$fit_obj)     # validate 1st model
karma.cv(kmodels[[2]]$fit_obj)     # validate 2nd model
                                   # ...etc...
```


<br/>

__Or you can follow the examples on this script:__

https://github.com/snarf-snarf/karma/blob/master/karma_demo.r

Use help() on any function to see a more detailed documentation.
