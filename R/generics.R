## basic constructor
setGeneric("pomp",function(data,...)standardGeneric("pomp"))

setGeneric("print",function(x,...)standardGeneric("print"))
setGeneric("plot",function(x,y,...)standardGeneric("plot"))
setGeneric("summary",function(object,...)standardGeneric("summary"))
setGeneric("window",function(x,...)standardGeneric("window"))

## constituent components of a 'pomp' object
setGeneric("dmeasure",function(object,...)standardGeneric("dmeasure"))
setGeneric("rmeasure",function(object,...)standardGeneric("rmeasure"))
setGeneric("dprocess",function(object,...)standardGeneric("dprocess"))
setGeneric("rprocess",function(object,...)standardGeneric("rprocess"))
setGeneric("dprior",function(object,...)standardGeneric("dprior"))
setGeneric("rprior",function(object,...)standardGeneric("rprior"))
setGeneric("init.state",function(object,...)standardGeneric("init.state"))
setGeneric("skeleton",function(object,...)standardGeneric("skeleton"))

## functions to extract or call the components of a "pomp" object
setGeneric("obs",function(object,...)standardGeneric("obs"))
setGeneric("data.array",function(object,...)standardGeneric("data.array"))
setGeneric("time",function(x,...)standardGeneric("time"))
setGeneric("time<-",function(object,...,value)standardGeneric("time<-"))  
setGeneric("coef",function(object,...)standardGeneric("coef"))
setGeneric("coef<-",function(object,...,value)standardGeneric("coef<-"))
setGeneric("states",function(object,...)standardGeneric("states"))
setGeneric("timezero",function(object,...)standardGeneric("timezero"))
setGeneric("timezero<-",function(object,...,value)standardGeneric("timezero<-"))
setGeneric("partrans",function(object,params,dir=c("forward","inverse"),...)standardGeneric("partrans"))
setGeneric("logLik",function(object,...)standardGeneric("logLik"))

## internals
setGeneric("pomp.fun",function(f,...)standardGeneric("pomp.fun"))
setGeneric("plugin.handler",function(object,...)standardGeneric("plugin.handler"))

## prediction mean
setGeneric("pred.mean",function(object,...)standardGeneric("pred.mean"))
## prediction variance
setGeneric("pred.var",function(object,...)standardGeneric("pred.var"))
## filter mean
setGeneric("filter.mean",function(object,...)standardGeneric("filter.mean"))
## conditional log likelihood
setGeneric("cond.logLik",function(object,...)standardGeneric("cond.logLik"))
## effective sample size
setGeneric("eff.sample.size",function(object,...)standardGeneric("eff.sample.size"))
## convergence record
setGeneric("conv.rec",function(object,...)standardGeneric("conv.rec"))
## stochastic simulation
setGeneric("simulate",function(object,nsim=1,seed=NULL,...)standardGeneric("simulate"))

## deterministic trajectory computation
setGeneric("trajectory",function(object,...)standardGeneric("trajectory"))
## trajectory matching
setGeneric("traj.match.objfun",function(object,...)standardGeneric("traj.match.objfun"))
setGeneric("traj.match",function(object,...)standardGeneric("traj.match"))

## ABC algorithm functions
setGeneric('abc',function(object,...)standardGeneric("abc"))

## Bayesian SMC (Liu & West)
setGeneric("bsmc",function(object,...)standardGeneric("bsmc"))

## basic SMC (particle filter)
setGeneric("pfilter",function(object,...)standardGeneric("pfilter"))

## particle Markov chain Monte Carlo (PMCMC)
setGeneric('pmcmc',function(object,...)standardGeneric("pmcmc"))

## nonlinear forecasting
setGeneric('nlf',function(object,...)standardGeneric("nlf"))

## iterated filtering
setGeneric('mif',function(object,...)standardGeneric("mif"))
## generate new particles
setGeneric('particles',function(object,...)standardGeneric("particles"))

## synthetic likelihood
setGeneric("probe",function(object,probes,...)standardGeneric("probe"))
## probe matching
setGeneric("probe.match.objfun",function(object,...)standardGeneric("probe.match.objfun"))
setGeneric("probe.match",function(object,...)standardGeneric("probe.match"))

## power spectrum
setGeneric("spect",function(object,...)standardGeneric("spect"))

## continue an iteration
setGeneric("continue",function(object,...)standardGeneric("continue"))
