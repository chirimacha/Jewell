#this code is for some new algorithms to simulate stochastic search
library(reshape2)
library(ggplot2)
library(pROC)

#TODO: set this to point to your code, or create an environment variable SPATIAL_UNCERTAINTY
#on Mac, run something like
#launchctl setenv SPATIAL_UNCERTAINTY "/Users/sgutfraind/academic_research/chagas/bandits"
#Note: The path above cannot have any spaces
#setwd("/home/sasha/Dropbox/chagas_models_aim3_spatial_uncertainty/code")
setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep=""))
source("bandit.R")

Cerro_Colorado = '
* very low prevalence
* spraying will
* perhaps reward based on finding clusters, rather than individual houses? (e.g. reward for any house within > 10m of existing house)
* using maps before spray, kids reporting as priors to bias the search
'



params_block <- list(
  sim_prob_shape1=0.8,
  sim_prob_shape2=0.8,
  denuncias_path="../data/surveillance/vigilancia_denuncias.csv",
  bases_tiabaya_casa_path="../data/bases_tiabaya/Tiabaya_Points_Casa-blocks.csv",
  byManz_fullEID_path="../data/corentinEID/byManz_fullEID.csv",
  byHouse_fullEID_path="../data/corentinEID/byHouse_fullEID.csv"
)

params_grid_sim <- list(
  max_cost=100,
  ncol=20,
  nrow=20,
  num_days=365,
  num_starts=1,
  p_hop=0.6,
  p_skip=0.3,
  p_jump=0.1,
  pMovePerDay=0.01, #per site 
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  ring=1,
  verbose=TRUE
)

params_manzana_sim <- list(
  ncol=10,
  nrow=3  
)

generate_infestation <- function(params) {
  #matrix with 1=infested, 0=not infested
  #1. some houses are initially infested.  
  #2. each turn there is a probability of starting a new infestation
  #assumptions: growth is linear in the infestation size (future: use Hussl and Riker growth models)
  #for (p in names(alt_params)) {
  #  stopifnot(p %in% names(params_def_liberia))
  #  params_liberia[p] <- alt_params[p]
  #}
  print("generating infestation ...")
  make_move <- function(name, house, nrow, ncol) {
    if(name == "hop") {
      stepsize <- 1
    } else if(name == "skip") {
      stepsize <- 2
    } else {
      stopifnot(name == "jump")
      stepsize <- 3
    }
    new_house <- house + round(runif(2)*stepsize)*sign(runif(2)-c(0.5, 0.5))
    ##wraparound - toroidal grid
    #new_house[1] <- 1 + ((new_house[1]-1)%%nrow)
    #new_house[2] <- 1 + ((new_house[2]-1)%%ncol)
    new_house$lat <- min(nrow, max(1, new_house$lat))
    new_house$lon <- min(ncol, max(1, new_house$lon))
    return(new_house)
  }
  
  infested_coords <- data.frame()
  for(trial in 1:params$num_starts) {
    new_house <- c(sample(1:params$nrow, 1), sample(1:params$ncol, 1))
    infested_coords <- rbind(infested_coords, new_house)
  }
  names(infested_coords) <- c("lat", "lon")
  
  infested_coords_prime <- infested_coords
  for(t in 1:params$num_days) {
    #clearing
    num_infested <- dim(infested_coords)[1]
    num_cleared <- rbinom(1, num_infested, params$pClearancePerDay)
    infested_coords_prime <- infested_coords[sample.int(num_infested, num_infested-num_cleared),]

    #could even add breeding ...
    #now move
    num_infested <- dim(infested_coords_prime)[1]
    num_spreading <- rbinom(1, num_infested, params$pMovePerDay)
    for(house_num in sample.int(num_infested, num_spreading)) {
      house <- infested_coords_prime[house_num,]
      r <- runif(1)
      if (r<params$p_hop) {
        new_house <- make_move("hop", house, params$nrow, params$ncol)
      } else if (r < params$p_hop + params$p_skip) {
        new_house <- make_move("skip", house, params$nrow, params$ncol)
      } else {
        new_house <- make_move("jump", house, params$nrow, params$ncol)
      }
      #print(new_house)
      infested_coords_prime <- rbind(infested_coords_prime, new_house)
    }
    infested_coords <- infested_coords_prime
  }
  #print(infested_coords)
  num_infested <- dim(infested_coords)[1]
  m <- matrix(0, nrow=params$nrow, ncol=params$ncol)
  for(house_num in 1:num_infested) {
    house <- unlist(infested_coords[house_num,])
    #weighted
    m[house[1], house[2]] <- m[house[1], house[2]] + 1
    
    #0,1
    #m[house[1], house[2]] <- 1
  }
  return(m)    
}

ranked_search <- function(infestation, st=NULL, max_actions=Inf, params=NULL) {
  #visited sites based on their ranks
  if(is.null(st)) {
    st <- list()
    st$total_squares <- length(infestation)
    st$visited       <- rep(0, length(infestation))
    st$running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0),current_score=c(Inf))
    
    st$known_infested  <- read.csv(text="idx")
    st$site_scores     <- params[["site_score"]](infestation)
    st$ordered_indices <- order(st$site_scores, decreasing=TRUE)
    #TODO: assumes that the ordering is unchanging.  instead should find max of note visited
    st$next_idx       <- 1        
    return(st)
  }
  
  initial_cost <- tail(st$running_stats[["total_cost"]], 1)
  next_stat <- tail(st$running_stats, 1)
  while(next_stat$total_visited < st$total_squares & 
          next_stat$total_cost    < initial_cost + max_actions) {
    next_site <- NULL
    #if(next_suspected_nb < dim(suspected_nbs)[1]){
    #  browser()
    #}      
    while (is.null(next_site) & st$next_idx <= st$total_squares) {
      next_site <- st$ordered_indices[st$next_idx]
      if (st$visited[next_site] == 0) {
        st$next_idx <- st$next_idx + 1
        break
      }
      st$next_idx <- st$next_idx + 1
    }
    if(is.null(next_site)) {
      break
    }
    next_stat$current_score <- st$site_scores[next_site]
    next_stat$total_cost    <- next_stat$total_cost + 1
    next_stat$total_visited <- next_stat$total_visited + 1
    st$visited[next_site] <- 1.0
    if(infestation[next_site] > 0) {
      st$visited[next_site] <- 1.0 #dim(st$running_stat)[1]/10.0
      st$known_infested     <- rbind(st$known_infested, next_site)
      next_stat$total_found <- next_stat$total_found + 1
      #print(sprintf("1, lon=%d, lat=%d", next_site$lon, next_site$lat))
    } else {
      #print(0)
    }
    #st$running_stats <- rbind(st$running_stats, next_stat)
    st$running_stats <- rbind(st$running_stats, next_stat)
    next_site <- NULL
  }
  st$running_stats$infestation_lower_bound <- st$running_stats$total_found/st$total_squares
  st$running_stats$infestation_estimate    <- st$running_stats$total_found/st$running_stats$total_visited
  row.names(st$running_stats) <- seq(dim(st$running_stats)[1])
  if(params[["verbose"]]) {
    print(tail(st$running_stats,1))
  }
  return(st)
}

plot_grid <- function(m=NULL, visited=NULL) {
  if(is.null(m)) {
    m <- matrix(rnorm(20),5)
  }
  infest_m <- melt(m) #uses dimensions Var1, Var2, value
  names(infest_m)<-c("lat", "lon", "infest")
  infest_m <- infest_m[infest_m$infest > 0,]

  if(is.null(visited)){ 
    visited <- matrix(0, nrow=dim(m)[1], ncol=dim(m)[2])
    visited[5,6] <- 1 #FIXME
  } else {
    stopifnot(dim(visited) == dim(m))
  }
  visited_m <- melt(visited)
  names(visited_m)<-c("lat", "lon", "visited")
  #visited_m <- visited_m[visited_m$visited > 0,]
  #all_data$visited <- melt(visited)$value
  pl <- ggplot(data=visited_m, aes(lon,lat)) + geom_raster(data=visited_m, aes(fill=visited, alpha=1.0)) + geom_point(data=infest_m, aes(lon,lat))
  #pl <- pl + geom_point(data=all_data, aes(shape=factor(visited)))
  #pl <- pl + geom_point(data=visited_m, aes(shape=factor(visited)))
  #pl <- pl + 
  
  print(pl)
}

ring_search1 <-function(infestation, st=NULL, params=NULL) {
  params2 <- params
  params2$ring <- 1
  return(ring_search(infestation=infestation, st=st, params=params2))
}

random_search <- function(infestation, st=NULL, max_actions=Inf, params=NULL, random_search=FALSE) {
  return(ring_search(infestation=infestation, st=st, max_actions=max_actions, params=params, random_search=TRUE))
}

ring_search <- function(infestation, st=NULL, max_actions=Inf, params=NULL, random_search=FALSE) {
#main search function. st = state.  random_search puts the code into purely random search (for benchmarking)
  if(is.null(st)) {
    st <- ring_search_initialize(infestation=infestation, params=params)
  }
  
  initial_cost <- tail(st$running_stats[["total_cost"]], 1)
  continue_inspection <- params[["continue_inspection"]]
  if(is.null(continue_inspection)) {
    continue_inspection <- function(infestation, latest, st) {
      return(latest$total_cost < params$max_cost & dim(latest)[1] < st$total_squares);
    }
  }
  next_stat <- tail(st$running_stats, 1)
  while(next_stat$total_visited < st$total_squares & 
        next_stat$total_cost    < initial_cost + max_actions & 
        continue_inspection(infestation, next_stat, st)) {
    #next_stat$step <- next_stat$step + 1
    next_site <- NULL
    #if(next_suspected_nb < dim(suspected_nbs)[1]){
    #  browser()
    #}      
    while (st$next_suspected_nb <= dim(st$suspected_nbs)[1] & (! random_search)) {
      next_nb <- st$suspected_nbs[st$next_suspected_nb,]
      st$next_suspected_nb <- st$next_suspected_nb + 1
      if(st$visited[next_nb$lat, next_nb$lon] == 0) {
        next_site <- next_nb
        break
      } else {
      }
    }
    while (is.null(next_site) & st$next_random_idx <= st$total_squares) {
      next_site <- st$randomized_sites[st$next_random_idx,]
      if (st$visited[next_site$lat, next_site$lon] == 0) {
        st$next_random_idx <- st$next_random_idx + 1
        break
      }
      st$next_random_idx <- st$next_random_idx + 1
    }
    if(is.null(next_site)) {
      break
    }
    next_stat$total_cost    <- next_stat$total_cost + 1
    next_stat$total_visited <- next_stat$total_visited + 1
    st$visited[next_site$lat, next_site$lon] <- 1.0
    if(infestation[next_site$lat, next_site$lon] > 0) {
      st$visited[next_site$lat, next_site$lon] <- 2.0 #dim(st$running_stat)[1]/10.0
      st$known_infested     <- rbind(st$known_infested, next_site)
      next_stat$total_found <- next_stat$total_found + 1
      neighbors             <- ring_search_get_unvisited_nbs(next_site, st$visited, ring=params$ring)
      st$suspected_nbs      <- rbind(st$suspected_nbs, neighbors)
      rownames(st$suspected_nbs) <- seq(dim(st$suspected_nbs)[1])
      #print(sprintf("1, lon=%d, lat=%d", next_site$lon, next_site$lat))
    } else {
      #print(0)
    }
    #st$running_stats <- rbind(st$running_stats, next_stat)
    st$running_stats <- rbind(st$running_stats, next_stat)
    next_site <- NULL
  }
  st$running_stats$infestation_lower_bound <- st$running_stats$total_found/st$total_squares
  st$running_stats$infestation_estimate    <- st$running_stats$total_found/st$running_stats$total_visited
  row.names(st$running_stats) <- seq(dim(st$running_stats)[1])
  if(params[["verbose"]]) {
    print(tail(st$running_stats,1))
  }
  return(st)
}

ring_search_get_unvisited_nbs <- function(house, visited, ring=2) {
  nbs <- read.csv(text="lon,lat")
  for(x in seq(house$lon-ring, house$lon+ring)) {
    if (x < 1 | x > params$ncol) {
      next
    }
    for(y in seq(house$lat-ring, house$lat+ring)) {
      if (y < 1 | y > params$nrow) {
        next
      }
      nb = list(lon=x,lat=y)  #wishlist: prioritize by range
      if (all(nb == house)) {
        next
      }
      if (visited[nb$lat,nb$lon] > 0) {
        next
      }
      nbs <- rbind(nbs, nb)
    }
  }
  return(nbs)
}

ring_search_initialize <- function(infestation, params=NULL) {
  st <- list()
  st$total_squares <- dim(infestation)[1] * dim(infestation)[2]
  st$visited <- matrix(0, nrow=dim(infestation)[1], ncol=dim(infestation)[2])
  st$randomized_sites <- melt(st$visited)
  st$randomized_sites <- st$randomized_sites[sample(st$total_squares),]
  names(st$randomized_sites)<-c("lat", "lon", "val")
  st$randomized_sites$val <- NULL
  st$running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0))
  
  st$known_infested <- read.csv(text="lon,lat")
  st$suspected_nbs  <- read.csv(text="lon,lat")
  st$next_suspected_nb <- 1
  st$next_random_idx <- 1    
  
  return(st)
}

z_infestation_search_fixed_algorithm <- function(random_search = FALSE) {
  infestation_test_m <- generate_infestation(params)
  search_stats       <- ring_search(infestation_test_m, st=NULL, max_actions=1, params=params)
  search_stats       <- ring_search(infestation_test_m, st=search_stats, max_actions=0, params=params)
  
  #search_stats       <- ring_search(infestation_test_m, st=NULL, params=params)
  visited            <- search_stats$visited
  #print(search_stats$running_stats)
  plot_grid(infestation_test_m, visited=visited)
  #plot_grid(infestation_test_m, visited=NULL)  
}

#z_infestation_fixed_algorithm(random_search = FALSE)
#z_infestation_fixed_algorithm(random_search = TRUE)

z_bandit_grid_search <- function(test_time=100, params=NULL) {
  source("bandit.R")
  
  infestation_test_m <- generate_infestation(params)
  n_arms <- 2
  search_stats <- list()
  
  #Arm 1 (random search) has much worse search ability than Arm #2 (ring search)
  funcs <- list(random_search, ring_search)
  
  search_stats[[1]] <- funcs[[1]](infestation_test_m, st=NULL, max_actions=0, params=params)
  search_stats[[2]] <- funcs[[2]](infestation_test_m, st=NULL, max_actions=0, params=params)
  
  
  pull_arm <- function(chosen_arm, search_stats) {
    st     <- search_stats[[chosen_arm]]
    new_st <- funcs[[chosen_arm]](infestation_test_m, st=st, max_actions=1, params=params)
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  search_reward <- function(new_st) {
    last_two <- tail(new_st$running_stats$total_found, 2)
    return(last_two[2] - last_two[1])
  }
  
  print("Stochastic search with Reinforcement-Comparison (rc) bandit");
  times   <- seq(1,test_time)
  arms    <- rep(0,length(times))
  rewards <- rep(0,length(times))
  bandit <- initialize_rc(n_arms=n_arms)
  print(bandit)
  for (trial in times) {
    ps <- NULL
    for (arm_idx in seq(n_arms)) {
      ps <- c(ps, probability_arm(bandit, arm_idx))
    }
    cat("trial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm <- next_arm(bandit)
    search_stats <- pull_arm(chosen_arm, search_stats)
    reward       <- search_reward(search_stats[[chosen_arm]])
    cat("  arm: ", chosen_arm, "  reward", reward, "\n")
    bandit <- update_bandit_rc(bandit, chosen_arm, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
  }
  results <- data.frame(T=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards))
  write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)  
}

manzanas_simulate_initial <- function(param) {
  #FIXME: construct a simulated infestation of blocks
  return(NULL)
}

time_stamp <- function(prefix="", suffix="", outmsg=TRUE) {
  #generates a unique time stamp for every data generated
  t <- format(Sys.time(), "%Y-%m-%d__%H-%M-%S");
  s <- as.integer(runif(1, max=1000))
  filename <- paste(prefix, t, s, suffix, sep="")
  if (outmsg) {
    print(filename)
  }
  return(filename)
}


z_bandit_ranked_search <- function(test_time=100, params=NULL) {
#application of MAB for ranked search of houses
  source("bandit.R")
  
  infestation_test_m <- manzanas_simulate_initial(params)
  n_arms <- 2
  search_stats <- list()
  
  #Arm 1 (random search) has much worse search ability than Arm #2 (ring search)
  funcs <- list(random_search, ring_search)

  #initialization
  search_stats[[1]] <- funcs[[1]](infestation_test_m, st=NULL, max_actions=0, params=params)
  search_stats[[2]] <- funcs[[2]](infestation_test_m, st=NULL, max_actions=0, params=params)
  
  
  pull_arm <- function(chosen_arm, search_stats) {
    st     <- search_stats[[chosen_arm]]
    new_st <- funcs[[chosen_arm]](infestation_test_m, st=st, max_actions=1, params=params)
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  search_reward <- function(new_st) {
    last_two <- tail(new_st$running_stats$total_found, 2)
    return(last_two[2] - last_two[1])
  }
  
  print("Stochastic search with Reinforcement-Comparison (rc) bandit");
  times   <- seq(1,test_time)
  arms    <- rep(0,length(times))
  rewards <- rep(0,length(times))
  bandit <- initialize_rc(n_arms=n_arms)
  print(bandit)
  for (trial in times) {
    ps <- NULL
    for (arm_idx in seq(n_arms)) {
      ps <- c(ps, probability_arm(bandit, arm_idx))
    }
    cat("trial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm <- next_arm(bandit)
    search_stats <- pull_arm(chosen_arm, search_stats)
    reward       <- search_reward(search_stats[[chosen_arm]])
    cat("  arm: ", chosen_arm, "  reward", reward, "\n")
    bandit <- update_bandit_rc(bandit, chosen_arm, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
  }
  results <- data.frame(T=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards))
  write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)  
}

z_data_stats <- function(params) {
  casas <- read.csv(params[["bases_tiabaya_casa_path"]])
  polys<-table(casas$polygon)
  block_sizes <- as.numeric(polys)
  hst <- hist(as.numeric(polys), breaks=c(0, 5, 10, 15, 20, max(block_sizes)))
  #block_sizes[order(block_sizes)]
  #manzana 334 has 204 houses!
  cat("Total number of manzanas: ", length(polys),"\n")  #361
  cat(hst$breaks, "  counts:\n")
  cat(hst$counts)
}

z_analyze_manzanas <- function(params) {
  #copied from Corentin
  #byManz$Den1<-byManz$id_manz %in% byHouse$id_manz[byHouse$nbDen1]
  byManz <- read.csv(params[["byManz_fullEID_path"]])
  
  #byHouse <- read.csv(params[["byHouse_fullEID_path"]])
  # A<-glm((Den1>0)~log(nbPos+1)*log(Unspray+1)*ageAP+log(Total+1)+as.factor(D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",A$aic,"\n")
  manz.model <-glm(Den1~log(nbPos+1)*log(Unspray+1)*ageAP+log(Total+1)+as.factor(D),data=byManz,family=binomial());
  print(manz.model)
  print("Odds Ratios:")
  #browser()
  cat("N:",dim(byManz)[1],"aic:",manz.model$aic,"\n")
  
  byManz$DenPredicted <- mean(manz.model$fitted.values)  #WARNING:  imputation
  prob_denuncias_predicted <- manz.model$fitted.values #predict.glm(manz.model, type="response")
  byManz[names(prob_denuncias_predicted),]$DenPredicted <- prob_denuncias_predicted
  byManz$LogDenPredicted <- -log(byManz$DenPredicted)
  
  #ADD IN RANKINGS - Highest density gets highest rank
  byManz$DenPredRank <- rank(byManz$DenPredicted, ties.method="first")
  
  write.csv(byManz, time_stamp("output/byManz", ".csv"))

  #diagnostics
  roccurve<-roc(manz.model$y~manz.model$fitted.values)
  auc(roccurve)

  #table(manz.model$y, ifelse(manz.model$fitted.values>0.14, 1, 0))
  #print(exp(coef(manz.model)))
  #print(exp(cbind(OR = coef(manz.model), confint(manz.model))))
  #print(summaryOR(manz.model),digits=3);

  nullmod <- glm(manz.model$y~1, family="binomial")
  mcfadden_pr2 <- 1-logLik(manz.model)/logLik(nullmod)
  print(mcfadden_pr2)
  # fine and probably easier to explain
}

z_sim_opt_in_manzana <- function(params) {
  
}

z_low_prevalence_experiments <- function() {
  n_arms=2
#execution with changed prevalence between S1 and S2
#   params <- list(
#     n=1000,
#     v=c(0.020, 0.015),
#     verbose=TRUE,
#   ) 
# alg_params <- list(
#   A=list(
#     site_score = function(x) {jitter(x + rpois(length(x), lambda=5))},
#     verbose=TRUE
#   ),
#   B=list(
#     site_score = function(x) {jitter(x + rpois(length(x), lambda=5))}, #same
#     verbose=TRUE
#   )
# )

  #another option - fixed prevalence but changed scoring
  params <- list(
    n=1000,
    v=c(0.02, 0.02),
    verbose=TRUE
  )
  alg_params <- list(
    A=list(
      site_score = function(x) {jitter(x + rpois(length(x), lambda=4))},
      verbose=TRUE
    ),
    B=list(
      site_score = function(x) {jitter(x + rpois(length(x), lambda=5))},
      verbose=TRUE
    )
  )
  arm_names <- c("A", "B")
  
  infestations <- list(
    A=rbinom(n=params[["n"]], size=1, prob=params[["v"]][1]),
    B=rbinom(n=params[["n"]], size=1, prob=params[["v"]][2]))
  
  print("Stochastic search with Reinforcement-Comparison (rc) bandit");
  times   <- seq(1,params[["n"]])
  arms    <- rep(0,length(times)) #which arm in which round
  rewards <- rep(0,length(times))
  pA      <- rep(0,length(times)) #probability of arm A
  mR      <- rep(0,length(times)) #mean reward across the algorithm, discounted

  bandit <- initialize_rc(n_arms=2, discount_factor=0.01, learning_rate=0.1)
  print(bandit)
  
  funcs <- list("A"=ranked_search, "B"=ranked_search)
  
  search_stats <- list(
    A=funcs[["A"]](infestations[["A"]], st=NULL, max_actions=1, params=alg_params[["A"]]),
    B=funcs[["B"]](infestations[["B"]], st=NULL, max_actions=1, params=alg_params[["B"]])
  )
  pull_arm <- function(chosen_arm, search_stats) {
    st     <- search_stats[[chosen_arm]]
    new_st <- funcs[[chosen_arm]](infestations[[arm_names[chosen_arm]]], st=st, max_actions=1, params=params)
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  
  search_reward <- function(new_st) {
    last_two <- tail(new_st$running_stats$total_found, 2)
    return(last_two[2] - last_two[1])
  }
  
  for (trial in times) {
    ps <- c()
    for (arm_idx in seq(n_arms)) {
      ps[arm_names[arm_idx]] <- probability_arm(bandit, arm_idx)
    }
    pA[trial] = ps["A"]
    cat("\n\nTrial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm   <- next_arm(bandit)
    search_stats <- pull_arm(chosen_arm, search_stats)
    reward       <- search_reward(search_stats[[arm_names[chosen_arm]]])
    cat("  arm: ", arm_names[chosen_arm], "  reward", reward, "\n")
    bandit       <- update_bandit_rc(bandit, chosen_arm, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")
    cat("  mean reward:", paste(bandit$mean_reward), "\n")
    mR[trial] <- bandit$mean_reward
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
  }
  wins.by.A <- search_stats[["A"]]$running_stats$total_found
  wins.by.A <- c(wins.by.A, Inf)  - c(0, wins.by.A)
  wins.by.A.timeline <- rep(-1, length(times))
  wins.by.A.timeline[which(arms == which(c("A", "B") == "A"))] <- wins.by.A[2:length(wins.by.A)-2]
  results <- data.frame(times=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards),
                        pA=pA,  mR=mR,  wins.by.A.timeline=wins.by.A.timeline)
  write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  
  #df_long <- melt(data.frame(times=results$times, pA=results$pA, wins.by.A.timeline=results$wins.by.A.timeline), 
  #                          id="times") # convert to long format
  #pl <- ggplot(df_long)
  #WORKING pl <- pl + geom_area(aes(color=variable, fill=variable), position = 'stack')
  
  pl <- ggplot()
  #pl <- pl + geom_line(aes(times, value, color=variable)) + xlab("Time")
  pl <- pl + geom_line(data=results, aes(x=times, y=pA), color="red") + xlab("Time") + ylab("")
  pl <- pl + geom_line(data=results, aes(x=times, y=mR), color="blue") + xlab("Time") + ylab("")
  pl <- pl + geom_point(data=results, aes(x=times, y=wins.by.A.timeline), color="black", size=0.5) + xlab("Time")
  pl <- pl + guides(color=FALSE, fill=FALSE)
  #pl <- pl +  theme(title=element_text(size=14,face="bold")) + 
  #        labs(variable="bold") + 
  #        scale_x_date(labels = date_format("%b%y")) + 
  #        theme(axis.title=element_text(size=15,face="bold"), 
  #        #          legend.position = c(.75, .7))
  #        legend.position = c(.75, .95))
  show(pl)
  ggsave(time_stamp("output/MAB_sim", ".png"), width=4, height=3, units="in")
  ggsave(time_stamp("output/MAB_sim", ".pdf"), width=4, height=3, units="in")
  
  #browser()
  return(results)  
}

z_low_prevalence_experiments()

#RUN RANKING CODE
#z_data_stats(params=params_block)
z_analyze_manzanas(params=params_block)



#z_infestation_search_fixed_algorithm()
#z_bandit_grid_search(params=params_grid_sim)



#z_sim_opt_in_manzana(params=params_manzana_sim)
