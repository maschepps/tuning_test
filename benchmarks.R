#ackley - 0, 0 = 0
#bohachevsky - 0, 0 = 0
#colville = 1, 1, 1, 1 = 0
##freudenstein_roth - 5, 4 = 0
#griewank_func = 0, 0 = 0
#miele_cantrell_func - 0, 1, 1, 1 = 0
#rastrigin_func = 0, 0 = 0
# rosenbrock_func - 1, 1 = 0
#schwefel_func - 420.9687, 420.9687 = -418.9829
#styblinski_tang_func - -2.903534 = -39.165

list_of_objectives = c(ackley_func,
  bohachevsky_func,
  colville_func,
  freudenstein_roth_func,
  griewank_func,
  miele_cantrell_func,
  rastrigin_func,
  rosenbrock_func,
  schwefel_func,
  styblinski_tang_func)

ackley_bounds = cbind(c(-32.76, 32.76), c(-32.76, 32.76)) 
bohachevsky_bounds = cbind(c(-100, 100), c(-100, 100))
colville_bounds =  cbind(c(-10, 10), c(-10, 10), c(-10, 10), c(-10, 10))
freudenstein_roth_bounds =  cbind(c(-10, 10), c(-10, 10))
griewank_bounds = cbind(c(-600, 600), c(-600, 600))
miele_cantrell_bounds =  cbind(c(-2, 2), c(-2, 2), c(-2, 2), c(-2, 2))
rastrigin_bounds =  cbind(c(-5.12, 5.12), c(-5.12, 5.12))
rosenbrock_bounds =  cbind(c(-5, 10), c(-5, 10))
schwefel_bounds = cbind(c(-500, 500), c(-500, 500))
styblinski_tang_bounds = cbind(c(-5, 5), c(-5, 5))

list_of_constraints =  list(ackley_bounds,
                         bohachevsky_bounds,
                         colville_bounds,
                         freudenstein_roth_bounds,
                         griewank_bounds,
                         miele_cantrell_bounds,
                         rastrigin_bounds,
                         rosenbrock_bounds,
                         schwefel_bounds,
                         styblinski_tang_bounds)



for(i in 1:length(list_of_objectives)){
  print(i)
  a = metaOpt(FUN = list_of_objectives[[i]],
              optimType = 'MIN',
              algorithm = 'GA',
              numVar    = ncol(list_of_constraints[[i]]),
              rangeVar  = list_of_constraints[[i]])
  print(a$result)
  
  p1 <- parameter("x1", list_of_constraints[[i]][1], list_of_constraints[[i]][2], FALSE)
  p2 <- parameter("x2", list_of_constraints[[i]][3], list_of_constraints[[i]][4], FALSE)
  p3 <- parameter("x3", list_of_constraints[[i]][5], list_of_constraints[[i]][6], FALSE)
  p4 <- parameter("x4", list_of_constraints[[i]][7], list_of_constraints[[i]][8], FALSE)
  
  conf_algo <- config_bat(iterations = 100, population_size = 100)
  params = c(list(p1, p2),
             list(p3, p4))
  if(ncol(list_of_constraints[[i]]) == 2){
    params = params[1:2]
  } else{
    params = params[1:4]
  }
  
  results <- minimize(algorithm_id = list("BAT"), 
                      obj_func = list_of_objectives[[i]], 
                      parameters = params,
                      config = conf_algo,
                      save_pop_history = T)
  print(results)
  
}

tic()
inner_ga_bayes = function(Pm){
  x = c(Pm)
  a = metaOpt(FUN = ackley,
              optimType = 'MIN',
              algorithm = 'GA',
              numVar    = 5,
              rangeVar  = rangeV,
              control = list(maxIter = 100, 
                             Pm = x[1],
                             Pc = x[2]))
  return(-a$optimumValue)
}

tryCatch({
  metaOpt(FUN =D1grp,
          optimType = "MIN",
          algorithm = 'GA',
          numVar = 6,
          rangeVar = t(parameter_bounds))
})
toc()

p1 <- parameter("x1", 0, 8, FALSE)
p2 <- parameter("x2", 0, 8, FALSE)
p3 <- parameter("x3", 0, 8, FALSE)
p4 <- parameter("x4", 0, 8, FALSE)
p5 <- parameter("x5", 0, 8, FALSE)
p6 <- parameter("x6", 0, 8, FALSE)
conf_algo <- config_bat(iterations = 100, population_size = 100)
results <- minimize(algorithm_id = list("BAT"), 
                    obj_func = D1grp, 
                    parameters = list(p1, p2, p3, p4, p5, p6),
                    config = conf_algo,
                    save_pop_history = T)
print(results)


