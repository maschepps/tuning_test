## Optimal enrollment design using upper caps in countries
## forecasting at study start-up 
## set-up
## min Cost
## given PoS >= Pplan

## Vlad Anisimov
## Sep 25 2021

## Using a normal approximation
## oriented to not so small number of countries >= 15
## to have a good approximation of PoS
## using a proportional representation of the main functions on vecN
## to accelerate time of computation


## PATH: C:/Users/vanisimo/Documents/R/Amgen21/OptCost-large-restrict-normal-approx-V8-1

# Input data

dataStudy <- read.csv("C:/Users/Admin/Desktop/Github/amg/OptCost-large-restrict-normal-approx-V8-1/OptCost-large-restrict-normal-approx-V8-1/Study-input-example-Caps-large.csv")

#dataStudy <- read.csv(file='C:/Users/vanisimo/Documents/R/Amgen21/OptCost-large-restrict-normal-approx-V8-1/Study-input-example.csv')

# Functions 
setwd('C:/Users/Admin/Desktop/Github/amg/OptCost-large-restrict-normal-approx-V8-1/OptCost-large-restrict-normal-approx-V8-1/')

source("functionsV14-meanvar.R")

source("functions-vecmeantimes-opt.R")

source("Fdistributionsitestimes.R")

source("functions-optim-normal-approx.R")


### CREATE input data for tool
# vecu;vecm;vectt;vecsi2

source("module-data-preparation-1.R")


### COSTS

## voluntary way 

AvcostPts <- 13000*seq(0.8,1.19,0.4/S)

## c(13600.76, 16821.52, 13683.21, 11936.29)

## SITES AND COUNTRY COSTS

AvcostSites <- rep(5000, length(AvcostPts))

AvcostCountries <- rep(20000, length(AvcostPts)) ## rep(100000, length(AvcostPts)) #rep(0, length(AvcostPts))

##
source("module-data-preparation-optim.R")


## FINAL FUNCTIONS TO USE IN OPTIMIZATION

#############################

## Use functions below to run optimization
## these functions can be used for any type of optimization algorithm


## the vectors vb,vs2,vAl,vprob,vecCap
## are precalculated and do not depend on the steps in optimization 

## New functions in a simplified form
## Global cost
## GC(vx)
## vx = vecNrun 

GC <- function(vx) GlobCostMeanPtsPG(vx,va,vb,vAl,vprob,vecCap) + FcostSiteCountry(vx,cs,cc)

## PoS in the interval [0,TT]
## FPoS(vx)  -- PoSNormal(vecN,vb,vs2,vAl,vprob,vecCap)

FPoS <- function(vx) PoSNormal(vx,vb,vs2,vAl,vprob,vecCap)  

library(metaheuristicOpt)
costFunctionC = function(param){
  # vecG = vecb + vecS2*vecmeantimes.22
  #constraints
  # sitesLow2 = sitesLow
  # sitesUpp2 = sitesUpp
  
  constraints = c(1:3)
  penalty = 0
  constraints[1] = any(param > vecUpp)
  constraints[2] = any(param < vecLow)
  # constraints[3] = nn - sum(vecb*cur_sites) + zq*sqrt(sum(cur_sites*vecG)) > 0
  constraints[3] = FPoS(param) < .9
  # constraints[3] = RunningPoS(debug_list$vecN, debug_list$vecCountries, debug_list$vecm, debug_list$vecsi2, debug_list$vecu, debug_list$vectt, debug_list$vecCap, debug_list$TT, debug_list$Mattimes, debug_list$vecSiteRate, debug_list$vecalfReg, debug_list$pts, unlist(debug_list$listuuN0))
  for(i in 1:3){
    penalty = sum(penalty, constraints[i])
  }
  obj1 = GC(round(param))
  return(obj1 + penalty * 100000000000000000000000000)
}

rangeV = matrix(t(data.frame(a = vecLow, b = vecUpp)), nrow = 2)
algos = c('GWO', 'MFO')

a = metaOpt(costFunctionC, optimType = 'MIN', algorithm = algos, numVar = length(vecUpp), rangeVar = rangeV, 
            control = list(numPopulation = 40, maxIter = 1000))


# a = metaOpt(costFunctionC, optimType = 'MIN', algorithm = 'PSO', numVar = length(vecUpp), rangeVar = rangeV, 
#             control = list(numPopulation = 40, maxIter = 1000))
library(doParallel)
registerDoParallel(20)
## usage
# vx = vecN
# GC(vx)
# FPoS(vx)


## Calculation of PoS for allocation - vecN 
vx <- vecN

before <- Sys.time()
PoSnorm <- FPoS(vx)
after <- Sys.time()
calc.time <- after-before
print(calc.time)
## Time difference of 0.01296496 secs

print(PoSnorm)


####################
####################

## The functions below can be used only for verification of PoS for particular allocations vecNrun
## using PG distributional approach 

source("functionsV14-prob.R")

## set of new functions
source("function-PoS-Caps-PG.R")

#vecNrun <- vecN

## PG distribution using proportional representation 

before <- Sys.time()
PoSPGprop <- FProbOfSuccessRestrictPG(pts,vecN,vAl,vprob,vecCap)
after <- Sys.time()
calc.time <- after-before
print(calc.time)

## Time difference of 2.446537 secs

## PoS using PG distribution
# print(ProbOfSuccessPG)

## running time is much longer, 188 times more, (2.446537 secs) compared to FPoS(vx) (0.01296496 secs)

## PoS - normal approximation
## using direct calculation of mean and Var of restricted PG processes without proportional representation
## using initial vecu created using vecN

PoSDirectNormal <- FPoSDirectNormal(pts,vecCap,vecCountries,vecu,vecm,vectt,vecsi2,veccountry,TT)
#print(PoSDirectNormal)

before <- Sys.time()
PoSDirectPG <- FPoSDirectPG(pts,vecCap,vecCountries,vecu,vecm,vectt,vecsi2,veccountry,TT)
after <- Sys.time()
calc.time <- after-before
print(calc.time)

print("PoS - norm. approx; PG proport; norm. not proport; PG direct - not proport")
print(c(PoSnorm,PoSPGprop,PoSDirectNormal,PoSDirectPG))
##


library(irace)

target_runner <- function(experiment, scenario){
  start = proc.time()
  configuration = experiment$configuration
  # par = runif(16, min = 0, max = 20)
  # a = PSO(costFunctionA, optimType = 'MIN', numVar = nrow(dataStudy), rangeVar = rangeV)
  a = genoud(fn = costFunctionC,
             starting.values = SitesUpp,
             nvars = length(SitesUpp),
             print.level = 0,
             MemoryMatrix = TRUE,
             Domains = cbind(as.numeric(SitesLow), as.numeric(SitesUpp)),
             boundary.enforcement = 2,
             P1 = configuration[['x']],
             P2 = configuration[['b']],
             P3 = configuration[['c']],
             P4 = configuration[['d']],
             P5 = configuration[['e']])
  time1 = proc.time() - start
  return(list(cost = a$value,
              time = time1[1]))
}
scenario <- list(targetRunner = target_runner,
                 instances = NULL,
                 maxExperiments = 200,
                 # Do not create a logFile
                 logFile = "")
parameters_table <- '
 x "" i (1, 500)
 b "" i (1, 500)
 c "" i (1, 500)
 d "" i (1, 500)
 e "" i (1, 500)
 '

## We use the irace function readParameters to read this table:
parameters <- readParameters(text = parameters_table)
setwd('C:/Users/Admin/Desktop/Github/tuning')
checkIraceScenario(scenario, parameters = parameters)
ss = irace(scenario = scenario, parameters = parameters)

library(tictoc)
tic()
a = genoud(fn = costFunctionC,
           starting.values = SitesUpp,
           nvars = length(SitesUpp),
           print.level = 0,
           MemoryMatrix = TRUE,
           boundary.enforcement = 2,
           Domains = cbind(as.numeric(SitesLow), as.numeric(SitesUpp)),
           P1 = 10,
           P2 = 10,
           P3 = 500,
           P4 = 100,
           P5 = 100)
toc()


