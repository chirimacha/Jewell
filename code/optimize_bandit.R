'
Functions for optimizing parameters
- currently implements toy optimization problem
- future: implement bandit

'

#https://cran.r-project.org/web/packages/neldermead/neldermead.pdf
library(neldermead)


objective_stochastic_toy <- function(x) {
#computes the objective.  the objectic is a parabola with a cusp giving the maximum, plus stochastic part which has a mean of zero
  true_max <- params[["true_max"]]
  ret <- (x[1]-true_max[1])^2 + (x[2]-true_max[2])^2 + sapply(c(1), runif) - 0.5
  return(ret)
}

objective_estimated <- function(x) {
#finds average of multiple replications of the stochastic objective.  the average is then used by the optimization
#for example, the return of a bandit would be stochastic, with anywhere between 0-20 house found out of 1000 searched at prevalence of 2%
#x as the form of a vector like c(2,3) expressing the parameter values  
  num_replications <- seq(params[["num_replications"]])
  returns <- c()
  for (r in seq(num_replications)) {
    returns <- c(returns, objective_stochastic_toy(x))
  }
  return(mean(returns))
}


z_test_evaluation <- function() {
  obj <- objective_estimated(x=c(3, 4))
  print(obj)
}

z_optize_grid <- function() {
  initial_guess <- c(0,0)
  
  number_pts_per_dimension <- 20  #ideally larger
  lower_bound <- c(-5,-5)
  upper_bound <- c(5,5)
  
  #for some odd reason, this function actually maximizes the objective, instead of Minimizing.  
  iteration_data <-fmin.gridsearch(fun = objective_estimated, 
                  x0 = initial_guess, 
                  xmin = lower_bound,
                  xmax = upper_bound, 
                  npts = number_pts_per_dimension)
  
  final_iteration <- tail(iteration_data, 1)
  xopt <-            final_iteration[,seq(length(initial_guess))]
  opt_is_feasible <- final_iteration$feasible
  
  cat("Best solution:", as.numeric(xopt), "\n")
  cat("   Feasible:", (opt_is_feasible>0), "\n")
}

z_optimize_neldermead <- function() {
#more advanced algorithm
  initial_guess <- c(0.5,0.5)
  
  lower_bound <- c(-5,-5)
  upper_bound <- c(5,5)
  
  iteration_data <-fminbnd(fun = objective_estimated, 
                                   x0 = initial_guess, 
                                   xmin = lower_bound,
                                   xmax = upper_bound) 
                                   #,options = list(MaxtIter=100))
  print(iteration_data$output)
  
  fopt <- iteration_data$optbase$fopt
  fx0 <- iteration_data$optbase$fx0
  xopt <- iteration_data$optbase$xopt
  #opt_is_feasible <- iteration_data$exitflag
  
  cat("Best solution: x=", as.numeric(xopt), "\n")
  cat("  f0:", fx0, " F_optimized: ", fopt, "\n")
  #cat("   Feasible:", (opt_is_feasible>0))
}

params <- list(num_replications=100, true_max=c(2,4))

z_optize_grid()    #good to try first
z_optimize_neldermead() #might never find good solutions, but if it does, could really tune it
z_optimize_neldermead() #restart a few times
z_optimize_neldermead() #restart a few times