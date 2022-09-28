?metaOpt
library(metaheuristicOpt)
algos = 'GA'
rangeV = rbind(c(-10, -10, -10, -10, -10),
               c(10, 10, 10, 10, 10))
ackley <- function(xx, a=20, b=0.2, c=2*pi)
{
  ##########################################################################
  #
  # ACKLEY FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # a = constant (optional), with default value 20
  # b = constant (optional), with default value 0.2
  # c = constant (optional), with default value 2*pi
  #
  ##########################################################################
  
  d <- length(xx)
  
  sum1 <- sum(xx^2)
  sum2 <- sum(cos(c*xx))
  
  term1 <- -a * exp(-b*sqrt(sum1/d))
  term2 <- -exp(sum2/d)
  
  y <- term1 + term2 + a + exp(1)
  return(y)
}

inner_ga = function(x){
  a = metaOpt(FUN = ackley,
          optimType = 'MIN',
          algorithm = 'GA',
          numVar    = 5,
          rangeVar  = rangeV,
          control = list(maxIter = 80, 
                         Pm = x[1],
                         Pc = x[2]))
  return(a$optimumValue)
}

rangeV_GA = rbind(c(0, 0),
               c(1, 1))
metaOpt(FUN = inner_ga,
        optimType = 'MIN',
        algorithm = 'GA',
        numVar = 2, 
        rangeVar = rangeV_GA,
        control = list(maxIter = 10))


metaOpt(FUN = ackley,
        optimType = 'MIN',
        algorithm = 'GA',
        numVar    = 5,
        rangeVar  = rangeV,
        control = list(maxIter = 5, 
                       Pm = 0.978,
                       Pc = 0.826))

library(irace)

target_runner <- function(experiment, scenario){
  start = proc.time()
  configuration = experiment$configuration
  # par = runif(16, min = 0, max = 20)
  # a = PSO(costFunctionA, optimType = 'MIN', numVar = nrow(dataStudy), rangeVar = rangeV)
  # a = genoud(fn = costFunctionC,
  #            starting.values = SitesUpp,
  #            nvars = length(SitesUpp),
  #            print.level = 0,
  #            MemoryMatrix = TRUE,
  #            Domains = cbind(as.numeric(SitesLow), as.numeric(SitesUpp)),
  #            boundary.enforcement = 2,
  #            P1 = configuration[['x']],
  #            P2 = configuration[['b']],
  #            P3 = configuration[['c']],
  #            P4 = configuration[['d']],
  #            P5 = configuration[['e']])
  # a = metaOpt(FUN = ackley,
  #             optimType = 'MIN',
  #             algorithm = 'GA',
  #             numVar    = 5,
  #             rangeVar  = rangeV,
  #             control = list(maxIter = 80, 
  #                            Pm = configuration[['m']],
  #                            Pc = configuration[['c']]))
  # 
  a = metaOpt(costFunctionC, 
              optimType = 'MIN', 
              algorithm = "PSO", 
              numVar = length(vecUpp), 
              rangeVar = rangeV, 
              control = list(numPopulation = 40,
                             ci = configuration[['i']],
                             cg = configuration[['g']],
                             w = configuration[['w']],
                             maxIter = 1000))
  time1 = proc.time() - start
  return(list(#cost = a$value,
              cost = a$optimumValue,
              time = time1[1]))
}

targetrunn
library(tictoc)
tic()
scenario <- list(targetRunner = target_runner,
                 instances = NULL,
                 maxExperiments = 1000,
                 parallel = 2,
                 # Do not create a logFile
                 logFile = "amgen_pso_race.Rdata"
                 )
parameters_table <- '
 i "" r (0, 100)
 g "" r (0, 100)
 w "" r (0, 100)
 '

## We use the irace function readParameters to read this table:
parameters <- readParameters(text = parameters_table)
setwd('C:/Users/Admin/Desktop/Github/tuning')
checkIraceScenario(scenario, parameters = parameters)
ss = irace(scenario = scenario, parameters = parameters)
toc()

load("iwantthis2.Rdata")
View(iraceResults$allConfigurations)
plot(iraceResults$allConfigurations$m, iraceResults$allConfigurations$c)
View(iraceResults$experiments)

res = cbind(t(iraceResults$experiments), iraceResults$allConfigurations)
res$means = rowMeans(res[,1:10], na.rm = T)
plot(res$m, res$means)
library(plotly)

plot_ly( x = res$m, y = res$c, z = res$means)


library(doParallel)
registerDoParallel(20)
analysis = foreach(i = 1:20, .combine = rbind, .packages = 'metaheuristicOpt') %dopar% 
  metaOpt(FUN = ackley,
              optimType = 'MIN',
              algorithm = 'GA',
              numVar    = 5,
              rangeVar  = rangeV,
              control = list(maxIter = 80, 
                             Pm = .1595418,
                             Pc = .1860242))
mean(unlist(analysis)[101:120])
sd(unlist(analysis)[101:120])

analysis2 = foreach(i = 1:20, .combine = rbind, .packages = 'metaheuristicOpt') %dopar% 
  metaOpt(FUN = ackley,
          optimType = 'MIN',
          algorithm = 'GA',
          numVar    = 5,
          rangeVar  = rangeV,
          control = list(maxIter = 80, 
                         Pm = .2175,
                         Pc = .3566))
mean(unlist(analysis2)[101:120])
sd(unlist(analysis2)[101:120])

analysis3 = foreach(i = 1:20, 
                    .combine = rbind, 
                    .packages = 'metaheuristicOpt') %dopar% 
  metaOpt(FUN = ackley,
          optimType = 'MIN',
          algorithm = 'GA',
          numVar    = 5,
          rangeVar  = rangeV,
          control = list(maxIter = 80))
mean(unlist(analysis3)[101:120])
sd(unlist(analysis3)[101:120])

analysis4 = foreach(i = 1:20, 
                    .combine = rbind, 
                    .packages = 'metaheuristicOpt') %dopar% 
  metaOpt(FUN = ackley,
          optimType = 'MIN',
          algorithm = 'GA',
          numVar    = 5,
          rangeVar  = rangeV,
          control = list(maxIter = 80,
                         Pm = runif(1, min = 0, max = 1),
                         Pc = runif(1, min = 0, max = 1)))
mean(unlist(analysis4)[101:120])
sd(unlist(analysis4)[101:120])

# install.packages('rBayesianOptimization')
library(rBayesianOptimization)
search_bound <- list(Pm = c(0,1), Pc = c(0,1))
search_grid = data.frame(Pm = runif(20, 0, 1),
                         Pc = runif(20, 0, 1))
BayesianOptimization(FUN = inner_ga,
                     bounds = search_bound,
                     init_grid_dt = search_grid,
                     n_iter = 10)


inner_ga_bayes = function(Pm){
  x = c(Pm)
  a = metaOpt(FUN = ackley,
              optimType = 'MIN',
              algorithm = 'GA',
              numVar    = 5,
              rangeVar  = rangeV,
              control = list(maxIter = 80, 
                             Pm = x[1],
                             Pc = .35))
  return(-a$optimumValue)
}
library(ParBayesianOptimization)
FUN = function(Pm, Pc) list(Score = inner_ga_bayes(Pm, Pc))  

optObjSimp = bayesOpt(FUN = FUN,
                      bounds = search_bound,
                      initGrid = search_grid,
                      # initPoints = 10,
                      iters.n = 1000,
                      iters.k = 1000,
                      errorHandling = "continue")

optObj <- list()
class(optObj) <- "bayesOpt"
optObj$FUN <- FUN
optObj$bounds <- search_bound
optObj$iters <- 0
iters.n = 10
iters.k = 1
optObj$initPars <- list()
optObj$optPars <- list()
optObj$GauProList <- list()
checkParameters(
  bounds
  , iters.n
  , iters.k
  , otherHalting
  , acq
  , acqThresh
  , errorHandling
  , plotProgress
  , parallel
  , verbose
)
# Formatting

otherHalting = list(timeLimit = Inf,minUtility = 0)
acq = "ucb"
kappa = 2.576
eps = 0.0
parallel = FALSE
gsPoints = pmax(100,length(bounds)^3)
convThresh = 1e8
acqThresh = 1.000
errorHandling = "stop"
plotProgress = FALSE
verbose = 1
initGrid = search_grid
boundsDT <- boundsToDT(bounds)
otherHalting <- formatOtherHalting(otherHalting)
setDT(initGrid)
inBounds <- checkBounds(initGrid,bounds)
inBounds <- as.logical(apply(inBounds,1,prod))
if (any(!inBounds)) stop("initGrid not within bounds.")
optObj$initPars$initialSample <- "User Provided Grid"
initPoints <- nrow(initGrid)



library(mlrMBO)
obj.fun = makeSingleObjectiveFunction(
  name = 'ackley',
  fn = function(x) inner_ga_bayes(Pm = x[1]
                                  ),
  par.set = makeParamSet(
    makeNumericVectorParam("x", len = 1L, lower = 0.0001, upper = 1)
),
minimize = F
)

des = generateDesign(n = 5, par.set = getParamSet(obj.fun), fun = lhs::randomLHS)
des$y = apply(des, 1, obj.fun)

surr.km = makeLearner("regr.km", predict.type = "se", covtype = "matern3_2", control = list(trace = FALSE))

control = makeMBOControl()
control = setMBOControlTermination(control, iters = 50)
# control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())

inner_pso_bayes = function(Pb, Pg){
  x = c(Pb, Pg)
  a = metaOpt(FUN = ackley,
              optimType = 'MIN',
              algorithm = 'PSO',
              numVar    = 5,
              rangeVar  = rangeV,
              control = list(maxIter = 5, 
                             ci = x[1],
                             cg = x[2]))
  return(-a$optimumValue)
}

inner_pso_bayes = function(Pb, Pg){
  x = c(Pb, Pg)
  number_of_parameters = length(x)
  parameter_bounds = cbind(rep(0, number_of_parameters),rep(1, number_of_parameters))
  a = optim_pso(objective_function = ackley, 
                number_of_parameters = 2, 
                number_of_particles = 40, 
                max_number_of_iterations = 10,
                max_number_function_calls= 1000, 
                parameter_bounds = parameter_bounds, 
                initial_estimates=NULL, 
                Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
                C1 = Pb, C2 = Pg,
                lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)
                
  return(a$value)
}

# result6_cov_1grp<-optim_pso(objective_function = Dopt_lith_1grp_SS3, number_of_parameters = number_of_parameters, number_of_particles = 40, 
#                             max_number_of_iterations = 1000,max_number_function_calls= 1000, 
#                             parameter_bounds = parameter_bounds, 
#                             initial_estimates=NULL, Vmax = (parameter_bounds[, 2] - parameter_bounds[, 1])/3, 
#                             lhc_init=FALSE, do_plot = NULL, wait_for_keystroke = FALSE, logfile = NULL, projectfile =NULL)


run = mbo(obj.fun, design = des, learner = surr.km, control = control, show.info = TRUE)
run = exampleRun(obj.fun, design = des, learner = surr.km, control = control, show.info = TRUE)
run = exampleRun(obj.fun, points.per.dim = 10, control = control, show.info = T)

plotExampleRun(run)

plot(run)
