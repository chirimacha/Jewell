'
Functions for optimizing parameters
- currently implements toy optimization problem
- future: implement bandit by writing a function similar to objective_stochastic_toy() below
(see zz_test_bandit_rc in bandit.R)

'

#https://cran.r-project.org/web/packages/neldermead/neldermead.pdf
library(neldermead)

#TODO: set this to point to your code, or create an environment variable SPATIAL_UNCERTAINTY

if(Sys.getenv("SPATIAL_UNCERTAINTY") != "") {
  setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep=""))  
} else if(Sys.getenv("USER") == "sgutfraind") {
  setwd("/Users/sgutfraind/academic/chagas/bandit/code")
} else if(Sys.getenv("USER") == "sasha") {
  setwd("/home/sasha/Documents/projects/chagas/shared/Bandit/code")
} else if(Sys.getenv("USER") == "mlevy") {
  #TODO: check your USER string (your username on your computer) and add appropriate directory here with 
  setwd("PATH_TO_BANDIT/code")
} else if(Sys.getenv("USER") == "rcastillo") {
  setwd("PATH_TO_BANDIT/code")
} else if(Sys.getenv("USER") == "snutman") {
  setwd("PATH_TO_BANDIT/code")
} 

#Call bandit simulations
source("bandit_simulations_for_optimization.R")

#timing
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  toc - tic
}

# objective_stochastic_toy <- function(x) {
# #computes the objective.  the objectic is a parabola with a cusp giving the maximum, plus stochastic part which has a mean of zero
#   true_max <- params[["true_max"]]
#   ret <- (x[1]-true_max[1])^2 + (x[2]-true_max[2])^2 + sapply(c(1), runif) - 0.5
#   return(ret)
# }

objective_estimated <- function(x) {
#finds average of multiple replications of the stochastic objective.  the average is then used by the optimization
#for example, the return of a bandit would be stochastic, with anywhere between 0-20 house found out of 1000 searched at prevalence of 2%
#x as the form of a vector like c(2,3) expressing the parameter values 
#x-c(learning_rate,discount_rate)  
  #assumptions: 
      'looking @500 houses (1/7 of tiabaya)
      10/houses per search
  '
  num_replications <- seq(params[["num_replications"]])
  returns <- c()
  for (r in seq(num_replications)) {
    returns <- c(returns, ZBanditGridProblemBlockedResults(test.time=50,params=params_grid_sim,
                                  params.arm=params.arm.1,new.infestations = "Yes",block.size=10,params.bandit=x))
  }
  return(mean(returns))
}


z_test_evaluation <- function(x) {
  obj <- objective_estimated(x)
  print(obj)
}


#grid optimization solver
'
@param guess:       vector x=c(a,b); a=learning rate, b=discount rate
@param points:      integer>=3: points per dimension. the larger the better. 
@param lower.bound: vector x=c(a,b) lower bound for parameters
@param lower.bound: vector x=c(a,b) upper bound for parameters
'

z_optize_grid <- function(guess,points,lower.bound,upper.bound) {

  
  initial_guess <- guess
  
  number_pts_per_dimension <- points  #ideally larger
  lower_bound <- lower.bound
  upper_bound <- upper.bound
  
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

#number of replications
params <- list(num_replications=20)

'
@param learning: bandit learning rate
@param discount: bandit discount rate

Bandit is updated with
F=preferences
r=reward
R=mean reward
A=discount factor
B=learning rate

F(t+1)=F(t)+B(r(t)-R(t)) <- update first
R(t+1)=A(r(t))+(1-A)(R(t))
'


#Infestation and search paramaters

'
@param max_cost:          limit for number of searches by RingSearch
@param ncol:              ncol in infestation
@param nrow:              nrow in infestation
@param num_days:          number of days of infestation if not using prevalence 
@param num_starts:        number of foci that begin each infestation [FIXME: Set to number of potential reinfested foci]
@param p_hop:             probability of hop
@param p_skip:            probability of skip
@param p_jump:            probability of jump
@param pMovePerDay:       probability of movement per site
@param pClearancePerday: spontaneous clearance per day
'
params_grid_sim <- list(
  max_cost=4000,
  ncol=25,
  nrow=45,
  num_days=365,
  num_starts=2,
  p_hop=0.8,
  p_skip=0.19,
  p_jump=0.01,
  pMovePerDay=0.01, #per site 
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  ring=1,
  verbose=FALSE
)

#params for arms 
params.arm.1 <- list (
  prev.arm1=.02,
  prev.arm2=.02,
  prev.arm3=.02
)

z_test_evaluation(c(.375,.25))

z_optize_grid (guess=c(0.1,0.1),points=10,lower.bound=c(0,0),upper.bound=c(1,1))
# Best solution: 0.375 0.25 


z_optize_grid()    #good to try first
z_optimize_neldermead() #might never find good solutions, but if it does, could really tune it
z_optimize_neldermead() #restart a few times
z_optimize_neldermead() #restart a few times
