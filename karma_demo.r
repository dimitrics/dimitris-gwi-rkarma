#---------------------------------------
# High-level functions:
#---------------------------------------
# auto.karma()        #<- Improved ARIMA and seasonal ARIMA selection algorithms (an alternative to auto.arima)
# auto.boxjen()       #<- Box-Jenkins model optimisation - improved automatic model selection (options for seasonal and non-seasonal data)
# karma.boxjenkins()  #<- The classic Box-Jenkins algorithm
# karma.cv()          #<- in-sample and out-of-sample cross-validation metrics and plots (compatible with arima, auto.arima, and nnetar)
# karma.ensemble()    #<- Train ensemble model (trains multiple weak learners) 
# karma.forecast()    #<- Forecasting unseen data using an ensemble model (aggregate multiple learners)
#
#---------------------------------------------------------------------------------
# Low-level functions:
#---------------------------------------
# karma.portmanteau() #<- Portmanteau test for autocorrelation
# karma.adf()         #<- Augmented Dickey-Fuller test for stationarity
# karma.transform()   #<- Detect transformation for wide-sense stationarity
# karma.fit()         #<- Fit ARIMA model to univariate time-series (basically a wrapper of Arima() - requires 'forecast' package)
# karma.cast()        #<- convert "ARIMA" object to "karma.fit" object
# karma.orderbin      #<- Convert fixed order terms to 0/NA mask
# karma.undiff()      #<- Reverse diff() - NOT IMPLEMENTED
# karma.validate()    #<- Options for model validation
# kforecastIn()       #<- in-sample forecast - UNUSED but useful for arma() from package 'tseries'
# kforecastOut()      #<- out-of-sample forecast - UNUSED
#
#---------------------------------------
# Optimisation API functions:
#---------------------------------------
# ModelEvaluationFunction()  #<- Objective function to model-fitting problem (MAPE minimisation)
# karma.hillclimbing()       #<- A configurable local-search hill-climbing type of algorithm with options to escape local optima
#
#---------------------------------------------------------------------------------



#library(devtools)
#install_github("snarf-snarf/karma")
library(forecast)
library(stringr)
library(karma)


#---------------------------------------------------
# Quick usage:
#---------------------------------------------------

kfit <- auto.karma(JohnsonJohnson)  #trains ARIMA/SARIMA model (selection and estimation)
karma.cv(kfit) #returns test-set MAPE to standard output
plot(forecast(kfit))  #make prediction


#---------------------------------------------------
# Basic functions:
#---------------------------------------------------

# Load series:
y = ldeaths   


# karma.cv: [works on both 'karma' and 'Arima' objects - #returns out-of-sample MAPE (default 80-20 validation) and plots]

karma.cv( auto.arima(y) )   #used with auto.arima() one-line model fit
karma.cv( nnetar(y) )         #used with nnetar() one-line model fit
karma.cv( auto.boxjen(y) )    #used with auto.boxjen() one-line Box-Jenkins model fit
karma.cv( auto.boxjen(y), cv = 'in' )   #in-sample option
plot(forecast(auto.boxjen(y))) #one-line forecast



# karma.fit():  [wraps Arima()]

karma.cv( karma.fit(y, order = c(10,1,5)) )
kfit = karma.fit(y, order = c(10,1,5), fixed = F)    #regular ARMA model
karma.cv(kfit, cv='in')  #in-sample MAPE and plot
karma.cv(kfit)    #out-of-sample MAPE and plot
plot(forecast(kfit)) #non-fixed model forecast

kfit = karma.fit(y, order = list(c(1,2,4,6), 1, c(1,3,6,7)), fixed = T)   #ARMA model with fixed terms
karma.cv(kfit)
plot(forecast(kfit)) #fixed model forecast


# karma.boxjenkins():   [automates classic Box-Jenkins approach - fixed terms or otherwise]

kfit = karma.boxjenkins(y, diffs = 1)   #one-line Box-Jenkins model
karma.cv(kfit)  #<- raw series (no stationarity transformation)
plot(forecast(kfit)) #forecast


# karma.transform(): [finds transformation steps for stationarity]

kfit = karma.boxjenkins(y, diffs = 1, log = T, fixed=T, max_iter = 100 )  #Box-Jenkins model with arbitrary selection of stationarity parameters
karma.cv(kfit)  #<- bad fit
plot(forecast(kfit)) #forecast

model1 = karma.transform(y, stdout = F, autolog = F, autodiffs = 1)  #Determine transfromation steps for stationarity
kfit = karma.boxjenkins(y, diffs = model1$diffs, log = model1$log, max_ar=15, max_ma=15, max_conv=2, max_iter=30, plot = F) #Box-Jenkins model with automatic selection of stationarity parameters
karma.cv(kfit)  #<- better fit 
plot(forecast(kfit)) #forecast



# auto.boxjen():  [automated Box-Jenkins model selection using forward selection, backwards elimination, or custom heuristic search options; NB: Configuration for best performance varies with input series ]

kfit = auto.boxjen(y)   #default: method = "greedy" (Box-jenkins mode | in-sample validation)
karma.cv(kfit)    
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "karma", metric = "AICc", optimiser = "semi-stochastic")   #method = "karma": combinatorial ARMA model selection (forward selection-type)
karma.cv(kfit)    
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "greedy-karma", optimiser = "semi-stochastic", metric = "AICc")   #backwards elimination starting with box-jenkins model as greedy solution (slower)
karma.cv(kfit)      
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "karma", optimiser = "stochastic", metric = "MAPE", cv = "in", max_rep = 1)   #Stochastic search - uses out-of-sample MAPE to validate every potential model
karma.cv(kfit)
plot(forecast(kfit)) #forecast



# auto.karma()   [alternative to auto.arima - automated ARIMA/SARIMA model selection]

#Good overall: (slower convergence)
sfit = auto.karma(y = y, method = "local-descent", cv = "in", stdout = T ); karma.cv(sfit) #sometimes better but generally slower convergence
sfit = auto.karma(y = y, method = "local-descent", cv = "out", stdout = T); karma.cv(sfit)

sfit = auto.karma(y = y, method = "random-walk", cv = "in", stdout = T, max_iter = 10); karma.cv(sfit)
sfit = auto.karma(y = y, method = "random-walk", cv = "out", stdout = T, max_iter = 10 ); karma.cv(sfit)

#Best overall: (faster convergence)
sfit = auto.karma(y = y, method = "markov-selection", cv = "in", stdout = F ); karma.cv(sfit)
sfit = auto.karma(y = y, method = "markov-selection", cv = "out", stdout = F ); karma.cv(sfit)

#More examples:
library(fpp2)
karma.cv(auto.karma(a10, test_type="auto", test_pct=-2), test_type="auto", test_pct=-2)
karma.cv(auto.karma(ausbeer, test_type="auto", test_pct=-5), test_type="auto", test_pct=-5)
karma.cv(auto.karma(beer, test_type="auto", test_pct=-1), test_type="auto", test_pct=-1)
karma.cv(auto.karma(auscafe, test_type="auto", test_pct=-3), test_type="auto", test_pct=-3)
karma.cv(auto.karma(austourists, test_type="auto", test_pct=-2), test_type="auto", test_pct=-2)
karma.cv(auto.karma(debitcards, test_type="auto", test_pct=-2), test_type="auto", test_pct=-2)
karma.cv(auto.karma(elecequip, test_type="auto", test_pct=-2), test_type="auto", test_pct=-2)
karma.cv(auto.karma(h02, test_type="auto", test_pct=-2), test_type="auto", test_pct=-2)


# karma.ensemble():  [creates ensemble of weak learners and aggregates their prediction on unseen data]

#Ex1:
kensemble = karma.ensemble(y = JohnsonJohnson, nsamples = 10, std_smoothing = 1)   #create ensemble of Box-jenkins models...
print(kensemble$unique_model_count)   #number of unique models in the ensemble
karma.cv(kensemble)   #model validation and out-of-sample MAPE
#karma.cv(kensemble$kfit_list[[1]])    #mape of the first model in the ensemble
kforecast = karma.forecast(kensemble) #mape of the aggregated forecast

#Ex2:
kensemble1 = karma.ensemble(y = JohnsonJohnson, nsamples = 10, fixed = T, metric = "MAPE", std_smoothing = 0) #ensemble of fixed terms ARMA models
print(kensemble1$unique_model_count)
karma.cv(kensemble1)
kforecast = karma.forecast(kensemble1) 

#Ex3:
kensemble1 = karma.ensemble(y = JohnsonJohnson, family = "sarima") #ensemble of fixed terms ARMA models
print(kensemble1$unique_model_count)
karma.cv(kensemble1)
kforecast = karma.forecast(kensemble1)


#---------------------------------------------------
# Box-Jenkins models:
#---------------------------------------------------

y = WWWusage

kfit = auto.boxjen(y, method = "karma", metric = "AICc")   #Optimise using in-sample MAPE
karma.cv(kfit, cv = 'in')  #<- will certainly overfit
karma.cv(kfit)   #<- proof of overfitting
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "karma", metric = "MSE", cv = "out", ac_criterion = T, autolog = F, optimiser = "semi-stochastic")   #Use portmonteau test as an extra optimisation constraint (ac_criterion=T) and MSE as the optimisation criterion
karma.cv(kfit)   #<- works good with WWWusage, underfits with JohnsonJohnson
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "greedy-karma", metric = "MAPE", cv = "out", ac_criterion = T, autolog = F, optimiser = "stochastic")   #Some further options (not suited for WWWusage)
karma.cv(kfit)  
plot(forecast(kfit)) #forecast

kfit = auto.boxjen(y, method = "karma", optimiser = "stochastic", mutations = T, max_rep = 15)   #Apply random mutations (random mutations are added to neighbours as local solutions are exhausted; it increases with c and stops at max_rep)
karma.cv(kfit) #Maybe the best setting for WWWusage
plot(forecast(kfit)) #forecast

# Manually further improve an solution: [using the lower level functions as an optimisation API]
kfit = auto.boxjen(y, method = "karma", optimiser = "stochastic", mutations = T, max_rep = 15, max_ar = 4, max_ma = 4)  #Just an initial solution...
tmp = karma.hillclimbing( sol = kfit$solopt, y = y, fixed=T, diffs = kfit$model_terms[[2]], optimiser = "semi-stochastic")  #Continue search using the initial solution as a staring point...
tmp = karma.hillclimbing( sol = tmp[[2]], y = y, fixed=T, diffs = kfit$model_terms[[2]],  optimiser = "semi-stochastic")  #Continue search from last solution using different settings...
#etc...



#---------------------------------------------------
#Algorithm comparison on multiple datasets:
#----------------------------------------------------

# One-line examples with various datasets: [WWWusage, austres, nhtemp, mdeaths, Nile, nottem, precip, presidents, Seatbelts, sunspot.month, sunspot.year, sunspots, treering, lynx, ldeaths, LakeHuron, EuStockMarkets, discoveries, BJsales, BJsales.lead]

library(fpp2)


#-- Lynx ----------------

#auto.arima
afit = auto.arima(lynx)
karma.cv(afit)
plot(forecast(afit))

#Greedy boxjenkins
bfit = auto.boxjen(lynx)
karma.cv(bfit)
plot(forecast(bfit))

#karma search
kfit = auto.boxjen(lynx, method="karma", cv="in", metric = "AICc")
karma.cv(kfit)
plot(forecast(kfit))


#-- nottem ----------------

#auto.arima
afit = auto.arima(nottem)
karma.cv(afit)
plot(forecast(afit))

#Greedy boxjenkins
bfit = auto.boxjen(nottem)
karma.cv(bfit)
plot(forecast(bfit))

#karma search
kfit = auto.boxjen(nottem, method="karma", cv="in", metric = "AICc")
karma.cv(kfit)
plot(forecast(kfit))


#-- h02 ----------------

#auto.arima
afit = auto.arima(h02)
karma.cv(afit)
plot(forecast(afit))

#Greedy boxjenkins
bfit = auto.boxjen(h02)
karma.cv(bfit)
plot(forecast(bfit))

#karma search
kfit = auto.boxjen(h02, method="karma", cv="in", metric = "AICc")
karma.cv(kfit)
plot(forecast(kfit))


#-- beer ----------------

#auto.arima
afit = auto.arima(beer)
karma.cv(afit)
plot(forecast(afit))

#Greedy boxjenkins
bfit = auto.boxjen(beer)
karma.cv(bfit)
plot(forecast(bfit))

#karma search
kfit = auto.boxjen(beer, method="karma", cv="in", metric = "AICc")
karma.cv(kfit)
plot(forecast(kfit))


#-- Other datasets: ----------------

karma.cv( auto.arima(austres) )
karma.cv( auto.karma(austres) )
karma.cv( auto.boxjen(y = austres, method = 'greedy') )

karma.cv( auto.arima(mdeaths) )
karma.cv( auto.karma(mdeaths) )
karma.cv( auto.boxjen(y = mdeaths, method = 'greedy') )

plot(forecast( auto.arima(Seatbelts[,1]) ))
plot(forecast( auto.karma(Seatbelts[,1]) ))
plot(forecast( auto.boxjen(y = Seatbelts[,1]) ))

plot(forecast( auto.arima(sunspot.year) ))
plot(forecast( auto.boxjen(sunspot.year) ))

plot(forecast( auto.arima(BJsales) ))
plot(forecast( auto.boxjen(BJsales) ))

#[from fpp2: goog, hyndman, pigs, oil]

#----------------------------------------------------




