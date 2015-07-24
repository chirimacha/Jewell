library(Rsymphony)

## Simple mixed integer linear program.
## maximize:    3 x_1 + 1 x_2 + 3 x_3
## subject to: -1 x_1 + 2 x_2 +   x_3 <= 4
##                      4 x_2 - 3 x_3 <= 2
##                x_1 - 3 x_2 + 2 x_3 <= 3
##                x_1, x_3 are non-negative integers
##                x_2 is a non-negative real number

#obj <- c(3, 1, 3)
#mat <- matrix(c(-1, 0, 1, 2, 4, -3, 1, -3, 2), nrow = 3)
#dirsign <- c("<=", "<=", "<=")
#rhs <- c(4, 2, 3)
#max <- TRUE
#types <- c("I", "C", "I")
#
#Rsymphony_solve_LP(obj, mat, dirsign, rhs, types = types, max = max)
#
### Same as before but with bounds replaced by
### -Inf <  x_1 <= 4
###    0 <= x_2 <= 100
###    2 <= x_3 <  Inf
#
#bounds <- list(lower = list(ind = c(1L, 3L), val = c(-Inf, 2)),
#               upper = list(ind = c(1L, 2L), val = c(4, 100)))
#
#Rsymphony_solve_LP(obj, mat, dirsign, rhs, types = types, max = max,
#                   bounds = bounds)
#

optimize_hunt <- function(strategies, params) {
#based on the stragies and parameters, formulate the problem as a MIP in matrix form, and solve it
  #order of the variables
  #H1, s_1^1, ...,  s_n^1, H2, s_1^2, ..., s_n^2, 
  num_strategies <- length(strategies)
  num_variables  <- (params$Tmax)*(num_strategies+1)
  types          <- c() #rep("C", num_variables)
  cat("Optimizing.  Num strategies: ", num_strategies, " Num weeks: ", params$Tmax, "\n")
  Sit_index      <- matrix(0, nrow=num_strategies, ncol=params$Tmax)
  Ht_index       <- rep(0,params$Tmax)  #Ht = number of HY sites available for t+1.
  next_idx <- 1
  for(t in 1:params$Tmax) {
    Ht_index[t] <- next_idx
    next_idx <- next_idx + 1
    types <- append(types, "C")
    for(i in 1:num_strategies) {
      Sit_index[i,t] <- next_idx
      next_idx <- next_idx + 1
    }
    types <- c(types, as.character(strategies$vartype))
  }
  #print(types)
  #print(Sit_index)
  #print(Ht_index)
  num_weeks      <- params$Tmax
  obj <- rep(0, num_variables)
  for(i in 1:num_strategies) {
    si_gain <- (params$Yh)*(strategies[i,]$HYvisits) + (params$Yl)*(strategies[i,]$LYvisits)
    for(t in 1:params$Tmax) {
      obj[Sit_index[i,t]] <- si_gain
    }
  }
  #print(obj)
  
  num_constraints <- 0  
  constraints <- matrix(0, nrow=(params$Tmax+params$Tmax), ncol=num_variables)
  rhs <- c()
  dirsign <- c()
  #weekly budgets
  for(t in 1:params$Tmax) {
    constr_line <- rep(0, num_variables)
    for(i in 1:num_strategies) {
      constr_line[Sit_index[i,t]] <- strategies[i,]$cost
    }    
    dirsign <- append(dirsign, "<=")
    rhs <- append(rhs, params$weeklybudget)
    num_constraints <- num_constraints + 1
    constraints[num_constraints,] <- constr_line
    cat(constr_line, dirsign[[num_constraints]], rhs[[num_constraints]], "\n")
  }
  #HY bags
  ## t=1
  constr_line <- rep(0, num_variables)
  for(i in 1:num_strategies) {
    constr_line[Sit_index[i,1]] <- -1*strategies[i,]$HYvisits
  }
  constr_line[Ht_index[1]] <- -1
  dirsign <- append(dirsign, "==")
  rhs <- append(rhs, -1*params$H0)
  num_constraints <- num_constraints + 1
  constraints[num_constraints,] <- constr_line
  cat(constr_line, dirsign[[num_constraints]], rhs[[num_constraints]], "\n")
  ## t>=2
  for(t in 2:params$Tmax) {
    constr_line <- rep(0, num_variables)
    for(i in 1:num_strategies) {
      constr_line[Sit_index[i,t-1]] <- strategies[i,]$FoundHY
      constr_line[Sit_index[i,t]] <- -1*strategies[i,]$HYvisits
    }
    constr_line[Ht_index[t-1]] <- 1
    constr_line[Ht_index[t]]   <- -1
    dirsign <- append(dirsign, "==")
    rhs <- append(rhs, 0)
    num_constraints <- num_constraints + 1
    constraints[num_constraints,] <- constr_line
    cat(constr_line, dirsign[[num_constraints]], rhs[[num_constraints]], "\n")
  }
  mat <- matrix(constraints, nrow=num_constraints)
  
  #bounds
  lower_bounds <- rep(0, num_variables)
  upper_bounds <- rep(0, num_variables)
  for(t in 1:params$Tmax) {
    lower_bounds[Ht_index[t]] <- 0
    upper_bounds[Ht_index[t]] <- 1000 #FIXME
    for(i in 1:num_strategies) {
      lower_bounds[Sit_index[i,t]] <- 0
      upper_bounds[Sit_index[i,t]] <- strategies[i,]$MaxUnits
    }    
  }  
  bounds <- list(lower = list(ind=seq(num_variables), val=lower_bounds), upper = list(ind=seq(num_variables), val=upper_bounds))
  print(lower_bounds)
  print(upper_bounds)
  
  soln <- Rsymphony_solve_LP(obj, mat, dirsign, rhs, types = types, max = TRUE, bounds = bounds)
  return (list(soln=soln, Sit_index=Sit_index, Ht_index=Ht_index))
}

##########
print_soln <- function(soln, strategies, params) {
  num_strategies <- length(strategies)
  num_variables  <- (params$Tmax)*(num_strategies+1)
  num_weekly_vars <- (num_strategies+1)
  
  cat("Strategy names: \n", rownames(strategies), "\n")
  if(soln$soln$status == 0) {
    print("Solution found")
  } else {
    print("Solution NOT found")
  }
  
  var_vals <- soln$soln$solution
  soln_mat <- matrix(NA, ncol=num_weekly_vars, nrow=params$Tmax)
  for(t in 1:params$Tmax) {
    soln_mat[t,1] <- var_vals[(soln$Ht_index[t])]
    for(i in 1:num_strategies) {
      soln_mat[t,i+1] <- var_vals[(soln$Sit_index[i,t])]
    }
  }
  pretty_soln <-as.data.frame(soln_mat)
  names(pretty_soln) <- c("H", rownames(strategies))
  print(pretty_soln)
  
  #TODO: report the number of finds
  cat("Expected sites found: ", soln$soln$obj)
  
  #print(soln)
}
params <- c()
params$weeklybudget <- 5
params$Tmax         <- 4
params$Yh           <- 0.2
params$Yl           <- 0.02
params$H0           <- 3 #params$weeklybudget / 3  #sensible choice
params$max_house_visits <- 100
radio_ad              <- data.frame(cost=10, HYvisits=0,  LYvisits=0,   FoundHY=20, MaxUnits=2, vartype="I")
grid_survey           <- data.frame(cost=1, HYvisits=0.1, LYvisits=0.9, FoundHY=0.1, MaxUnits=params$max_house_visits, vartype="C")
high_only             <- data.frame(cost=1, HYvisits=1,   LYvisits=0,   FoundHY=2, MaxUnits=params$max_house_visits, vartype="C")
two_hops_from_infest  <- data.frame(cost=1, HYvisits=0.5, LYvisits=0.5, FoundHY=1, MaxUnits=params$max_house_visits, vartype="C")
coerce_missed         <- data.frame(cost=2, HYvisits=0.9, LYvisits=0,   FoundHY=0, MaxUnits=params$max_house_visits, vartype="C")
follow_passive_reports<- data.frame(cost=15, HYvisits=10, LYvisits=5,   FoundHY=20, MaxUnits=1, vartype="I")

strategies <- rbind(radio_ad=radio_ad, grid_survey=grid_survey, high_only=high_only, two_hops_from_infest=two_hops_from_infest, coerce_missed=coerce_missed, follow_passive_reports=follow_passive_reports)

solnB5 <- optimize_hunt(strategies=strategies, params=params)
print_soln(solnB5, strategies, params)
