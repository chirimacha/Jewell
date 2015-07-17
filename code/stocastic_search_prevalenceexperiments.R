'
################################
This code:
  -Calls Bandit Code 
  -Simulates infestations
  -Implements experiments looking at 
  bandit switching under various prevalence conditions
- implements experiments looking at 
  bandit switching under "blocked" results
  (having multiple searches within an arm 
  before returning results to the bandit)

CREATOR: SN based on code from SG
######################################
'
library(reshape2)
library(ggplot2)
library(pROC)
library(sm)

number <- sample(1:10000, 1)
print(number)
set.seed(number)

TimeNow <- function() {
  x <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")
  return(x)
}


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

#SET OUTPUT PATHS
if(! file.exists("output")) {
  dir.create("output")  
}

#CALL BANDIT CODE 
source("Bandit.R")


params_block <- list(
  sim_prob_shape1=0.8,
  sim_prob_shape2=0.8
)

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

params.bandit <- list(
  learning =0.1,
  discount=0.01
)

params.bandit.2 <-list (
  learning=0.2,
  discount=0.1
)

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
  ncol=80,
  nrow=50,
  num_days=365,
  num_starts=1,
  p_hop=0.6,
  p_skip=0.3,
  p_jump=0.1,
  pMovePerDay=0.01, #per site 
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  ring=1,
  verbose=FALSE
)

params_grid_sim2 <- list(
  max_cost=4000,
  ncol=30,
  nrow=30,
  num_days=365,
  num_starts=2,
  p_hop=0.6,
  p_skip=0.3,
  p_jump=0.1,
  pMovePerDay=0.01, #per site 
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  ring=1,
  verbose=FALSE
)


#set parameters for number of arms and prevalence in each arm
params.arm.1 <- list (
  prev.arm1=.05,
  prev.arm2=.01
)

params.arm.2 <- list (
  prev.arm1=.02,
  prev.arm2=.03
)

params.arm.3 <- list (
  prev.arm1=.03,
  prev.arm2=.02,
  prev.arm3=.02
)

params.arm.4 <- list (
  prev.arm1=.025,
  prev.arm2=.02
)

params.arm.5 <- list (
  prev.arm1  =.02,
  prev.arm2  =.02,
  prev.arm3  =.02,
  prev.arm4  =.02,
  prev.arm5  =.02,
  prev.arm6  =.02,
  prev.arm7  =.02,
  prev.arm8  =.02,
  prev.arm9  =.02,
  prev.arm10 =.02,
  prev.arm11 =.02,
  prev.arm12 =.02,
  prev.arm13 =.03,
  prev.arm14 =.03,
  prev.arm15 =.03
)

params.arm.6 <- list (
  prev.arm1  =.02,
  prev.arm2  =.02,
  prev.arm3  =.02,
  prev.arm4  =.02,
  prev.arm5  =.02,
  prev.arm6  =.02,
  prev.arm7  =.02,
  prev.arm8  =.02,
  prev.arm9  =.02,
  prev.arm10 =.02,
  prev.arm11 =.02,
  prev.arm12 =.02,
  prev.arm13 =.02,
  prev.arm14 =.02,
  prev.arm15 =.02,
  prev.arm16 =.02,
  prev.arm17 =.02,
  prev.arm18 =.02,
  prev.arm19 =.02,
  prev.arm20 =.02,
  prev.arm21 =.02,
  prev.arm22 =.02,
  prev.arm23 =.02,
  prev.arm24 =.02,
  prev.arm25 =.02,
  prev.arm26 =.02,
  prev.arm27 =.02,
  prev.arm28 =.02,
  prev.arm29 =.02,
  prev.arm30 =.02,
  prev.arm31 =.02,
  prev.arm32 =.02,
  prev.arm33 =.03,
  prev.arm34 =.03,
  prev.arm35 =.03,
  prev.arm36 =.03,
  prev.arm37 =.03,
  prev.arm38 =.03,
  prev.arm39 =.03,
  prev.arm40 =.03
)

params.arm.7 <-list (
  prev.arm1 =.04,
  prev.arm2 =.02,
  prev.arm3=.02,
  prev.arm4=.02
)


GenerateInfestation <- function(params,prevalence_rule) {
  #infestation simulator  
  #matrix with 0=not infested, x>0= infested with x bugs
  #1. some houses are initially infested.  
  #2. each turn there is a probability of starting a new infestation
  #assumptions: growth is linear in the infestation size (future: use Hussl and Riker growth models)
  #for (p in names(alt_params)) {
  #  stopifnot(p %in% names(params_def_liberia))
  #  params_liberia[p] <- alt_params[p]
  #}
  print("generating infestation ...")
  MakeMove <- function(name, house, nrow, ncol) {
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
  
  #CREATE INITIAL INFESTATION
  infested_coords <- data.frame()
  for(trial in 1:params$num_starts) {
    new_house <- c(sample(1:params$nrow, 1), sample(1:params$ncol, 1))
    infested_coords <- rbind(infested_coords, new_house)
  }
  names(infested_coords) <- c("lat", "lon")
  
  infested_coords_prime <- infested_coords
 #could also create "time-based" infestation
  #for(t in 1:params$num_days) {
  
  #create prevalence-based infestation: 
    #This infestation stops once a certain prevalence is reached
  prevalence=0 #intialize prevalence 
  while(prevalence<prevalence_rule)  {
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
        new_house <- MakeMove("hop", house, params$nrow, params$ncol)
      } else if (r < params$p_hop + params$p_skip) {
        new_house <- MakeMove("skip", house, params$nrow, params$ncol)
      } else {
        new_house <- MakeMove("jump", house, params$nrow, params$ncol)
      }
      infested_coords_prime <- rbind(infested_coords_prime, new_house)
    }
    infested_coords <- infested_coords_prime
    
    #CALCULATE PREVALENCE (#INFECTED HOUSES/TOTAL HOUSES)
    prevalence <- ((dim (unique(infested_coords))[1])/(params$nrow*params$ncol))
    #print(prevalence)
  }
  #print(infested_coords)
  num_infested <- dim(infested_coords)[1]
  m <- matrix(0, nrow=params$nrow, ncol=params$ncol)
  for(house_num in 1:num_infested) {
    house <- unlist(infested_coords[house_num,])
    #print(house)
    #weighted
    m[house[1], house[2]] <- m[house[1], house[2]] + 1
    
    #0,1
    #m[house[1], house[2]] <- 1
  }
  return(m)
}


# RankedSearch <- function(infestation, st=NULL, max_actions=Inf, params=NULL) {
#   '
#   ranked search function to be used by the bandit.
#   1. the function completes up to max_actions searches
#   2. the function updates the state variable st.  st is passed to the function the next time this arm is pulled.
#   - infestation is the unchanging state of the infestation, which is not fully known to the function
#   - on the first pull, st==NULL.  later st is updated.  ie. like a rabbit, the function "eats its own output"
#   '
#   #visited sites based on their ranks
#   if(is.null(st)) {  #only called when the st==NULL
#     st <- list()
#     st$total_squares <- length(infestation)
#     st$visited       <- rep(0, length(infestation))
#     st$running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0),current_score=c(Inf))
#     
#     st$known_infested  <- read.csv(text="idx")
#     
#     #the scores of the sites.  the algorithm always selects the site with the max score
#     st$site_scores     <- params[["site_score"]](infestation)
#     return(st)
#   }
#   
#   initial_cost <- tail(st$running_stats[["total_cost"]], 1)
#   next_stat <- tail(st$running_stats, 1)
#   while(next_stat$total_visited < st$total_squares & 
#         next_stat$total_cost    < initial_cost + max_actions) {
#     #if(next_suspected_nb < dim(suspected_nbs)[1]){
#     #  browser()
#     #}      
#     next_site_score <- 0 
#     while (next_site_score != -Inf) {
#       next_site <- which.max(st$site_scores)
#       next_site_score <- st$site_scores[next_site]
#       if (st$visited[next_site] == 0) {
#         break
#       }
#       st$site_scores[next_site] <- -Inf  #never selected again
#     }
#     if(next_site_score == -Inf) {
#       break #all scores are -Inf.  This means that finished search
#     }
#     next_stat$current_score <- st$site_scores[next_site]
#     next_stat$total_cost    <- next_stat$total_cost + 1
#     next_stat$total_visited <- next_stat$total_visited + 1
#     st$visited[next_site] <- 1.0
#     if(infestation[next_site] > 0) {
#       st$visited[next_site] <- 1.0 #dim(st$running_stat)[1]/10.0
#       st$known_infested     <- rbind(st$known_infested, next_site)
#       next_stat$total_found <- next_stat$total_found + 1
#       #print(sprintf("1, lon=%d, lat=%d", next_site$lon, next_site$lat))
#     } else {
#       #print(0)
#     }
#     #st$running_stats <- rbind(st$running_stats, next_stat)
#     st$running_stats <- rbind(st$running_stats, next_stat)
#     next_site <- NULL
#   }
#   st$running_stats$infestation_lower_bound <- st$running_stats$total_found/st$total_squares
#   st$running_stats$infestation_estimate    <- st$running_stats$total_found/st$running_stats$total_visited
#   row.names(st$running_stats) <- seq(dim(st$running_stats)[1])
#   if(params[["verbose"]]) {
#     print(tail(st$running_stats,1))
#   }
#   return(st)
# }
# 
# PlotGrid <- function(m=NULL, visited=NULL) {
#   if(is.null(m)) {
#     m <- matrix(rnorm(20),5)
#   }
#   infest_m <- melt(m) #uses dimensions Var1, Var2, value
#   names(infest_m)<-c("lat", "lon", "infest")
#   infest_m <- infest_m[infest_m$infest > 0,]
#   
#   if(is.null(visited)){ 
#     visited <- matrix(0, nrow=dim(m)[1], ncol=dim(m)[2])
#     visited[5,6] <- 1 #FIXME
#   } else {
#     stopifnot(dim(visited) == dim(m))
#   }
#   visited_m <- melt(visited)
#   names(visited_m)<-c("lat", "lon", "visited")
#   #visited_m <- visited_m[visited_m$visited > 0,]
#   #all_data$visited <- melt(visited)$value
#   pl <- ggplot(data=visited_m, aes(lon,lat)) + geom_raster(data=visited_m, aes(fill=visited, alpha=1.0)) + geom_point(data=infest_m, aes(lon,lat))
#   #pl <- pl + geom_point(data=all_data, aes(shape=factor(visited)))
#   #pl <- pl + geom_point(data=visited_m, aes(shape=factor(visited)))
#   #pl <- pl + 
#   
#   print(pl)
# }
# 



RingSearch <- function(infestation, st=NULL, max_actions=Inf, params=NULL, random_search=FALSE) {
  '
  toy search function.
  1. randomly selects sites until makes a hit
  2. explores the ring around the known hit.  if exhausts all known rings, reverts to random search
  
  st = state.  on first call st==NULL which causes this to initialize using RingSearchInitialize
  random_search=TRUE puts the code into purely random search, i.e. not uses rings (for benchmarking)
  
  '
  if(is.null(st)) {
    st <- RingSearchInitialize(infestation=infestation, params=params)
  }
  
  initial_cost <- tail(st$running_stats[["total_cost"]], 1)
  ContinueInspection <- params[["ContinueInspection"]]
  if(is.null(ContinueInspection)) {
    ContinueInspection <- function(infestation, latest, st) {
      return(latest$total_cost < params$max_cost & dim(latest)[1] < st$total_squares);
    }
  }
  next_stat <- tail(st$running_stats, 1)
  while(next_stat$total_visited < st$total_squares & next_stat$total_cost < initial_cost + max_actions & 
        ContinueInspection(infestation, next_stat, st)) {
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
      } else {
        next_site <-NULL
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
      next_stat$total_bugs  <- next_stat$total_bugs + (infestation[next_site$lat, next_site$lon])
      neighbors             <- RingSearchGetUnvisitedNbs(next_site, st$visited, ring=params_grid_sim$ring, params=params)
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
  st$running_stats$unfoundprevalence <- (st$true.prevalence-st$running_stats$total_found)/(st$total_squares-st$running_stats$total_visited)
  
  if(params[["verbose"]]) {
    print(tail(st$running_stats,1))
  }
  if (exists("testing_parameters")) {
    testing_parameters<-rbind(testing_parameters,tail(st$running_stats,1))
    assign("testing_parameters",testing_parameters,envir=.GlobalEnv)
  }
  else {
    assign("testing_parameters",tail(st$running_stats,1),envir=.GlobalEnv)
  }
  
  return(st)
}

#Find unvisited houses in the ring around the infected house
RingSearchGetUnvisitedNbs <- function(house, visited,ring,params) {
  nbs <- read.csv(text="lon,lat")
  for(x in seq(house$lon-ring, house$lon+ring)) {
    if (x < 1 | x > params$ncol) {
      next
    }
    for(y in seq(house$lat-ring, house$lat+ring)) {
      if (y < 1 | y > params$nrow) {
        next
      }
      print(x)
      print(y)
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

#Initialize parameters for ring search 
RingSearchInitialize <- function(infestation, params=NULL) {
  st <- list()
  st$total_squares <- dim(infestation)[1] * dim(infestation)[2]
  st$true.prevalence <- sum(colSums(infestation !=0))
  st$visited <- matrix(0, nrow=dim(infestation)[1], ncol=dim(infestation)[2])
  st$randomized_sites <- melt(st$visited)
  st$randomized_sites <- st$randomized_sites[sample(st$total_squares),]
  names(st$randomized_sites)<-c("lat", "lon", "val")
  st$randomized_sites$val <- NULL
  st$running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0),total_bugs=c(0))
  
  st$known_infested <- read.csv(text="lon,lat")
  st$suspected_nbs  <- read.csv(text="lon,lat")
  st$next_suspected_nb <- 1
  st$next_random_idx <- 1    
  
  return(st)
}


ZPrevalenceRingTest <- function(params=NULL,prevalence_rule=NULL, num.trials=NULL) {
  '
  function for benchmarking how well the ring search does at given 
  prevalance levels
  '
  #rm(testing_parameters_all, envir=.GlobalEnv)
  
  infestation_test_1 <- GenerateInfestation(params=params,prevalence_rule=prevalence_rule)
  for (x in 1:num.trials) {
    for (t in 1:100) {
      test<-RingSearch(infestation_test_1, st=NULL, max_actions=100, params=params)
    }
    if (x==1) {
      nobugs <- sum(testing_parameters$infestation_estimate==0)
      testing_parameters_all<-testing_parameters
      testing_parameters_all$x<-x
      rm(testing_parameters)
    } else {
      nobugs <- rbind(nobugs,sum(testing_parameters$infestation_estimate==0))
      testing_parameters$x<-x
      testing_parameters_all<-rbind(testing_parameters_all,testing_parameters)
      rm(testing_parameters)
    }
  print(x)
  }
  d <- density(nobugs)
  plot(d)
  sm.density.compare(testing_parameters_all$infestation_estimate, testing_parameters_all$x, xlab=)
  assign("nobugs",nobugs,envir=.GlobalEnv)
  return(testing_parameters_all)
}


'function for applying the bandit to a given infestation, 
  using ring search and banditing on location 
@params test.time:        number of times to run the bandit
@params params:           grid parameters
@params params.arm:       arms and prevalence parameters
@params new.infestation:  generate new infestations or start with old ones [for benchmarking]
'
ZBanditGridProblem <- function(test.time=NULL, params=NULL, params.arm=NULL, new.infestations=NULL,
                               params.bandit=NULL) {

  
  n_arms <- length(params.arm)
  
  search_stats <- list()
  
  #GENERATE INFESTATIONS
  #generates new infestations for each arm of the bandit
  if (new.infestations=="Yes") {
  infestations<- lapply(params.arm,GenerateInfestation,params=params) 
  assign("static.infestations",infestations,envir=.GlobalEnv)
  }
  
  if (! exists("infestations")) {
    infestations <- static.infestations
    #print("newfromold")
  }
  
  #SET UP BANDIT ARMS; CURRENTLY USING RING SEARCH STRATEGY
  #initialize RingSearch function for all infestations 
  
  search_stats <- lapply(infestations, RingSearch, st=NULL, max_actions=0, params=params)
  
  #function to pull bandit
  pull_arm <- function(chosen_arm, search_stats) {
    st     <- search_stats[[chosen_arm]]
    new_st <- RingSearch(data.frame(infestations[[chosen_arm]]),st=st,max_actions=1,params=params) 
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  #function to find reward for bandit
  search_reward <- function(new_st) {
    last_two <- tail(new_st$running_stats$total_found, 2)
    return(last_two[2] - last_two[1])
  }
  #function to find remaining prevalence (for testing)
  UnfoundPrevalence <- function(new_st) {
    unfound <-tail(new_st$running_stats$unfoundprevalence,1)
    return(unfound)
  }
  print("Stochastic search with Reinforcement-Comparison (rc) bandit");
  times       <- seq(1,test.time)
  arms        <- rep(0,length(times))
  rewards     <- rep(0,length(times))
  mR          <- rep(0,length(times))
  unfoundprev <- rep(0,length(times))
  ps_holding <- matrix(nrow=length(times),ncol=length(params.arm))
  bandit <- initialize_rc(n_arms=n_arms, learning_rate=params.bandit$learning, discount_factor=params.bandit$discount)
  print(bandit)
  
  for (trial in times) {
    ps <- NULL
    for (arm_idx in seq(n_arms)) {
      ps <- c(ps, probability_arm(bandit, arm_idx))
    }
    cat("trial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm <- next_arm(bandit)
    #print(chosen_arm)
    search_stats <- pull_arm(chosen_arm, search_stats)
    reward       <- search_reward(search_stats[[chosen_arm]])
    unfound      <- UnfoundPrevalence(search_stats[[chosen_arm]])
    cat("  arm: ", chosen_arm, "  reward", reward, " unfound", unfound, "\n")
    bandit <- update_bandit_rc(bandit, chosen_arm, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
    ps_holding[trial,] <-ps
    mR[trial] <- bandit$mean_reward
    unfoundprev[trial] <- unfound
  }
  colnames(ps_holding)=params.arm
  results <- data.frame(ps_holding,T=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards), MeanReward=mR, UnfoundPrev.ChosenArm=unfoundprev)
  #write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)  
}

#function that tells you how many times the bandit was pulled in a given turn
GetResults <- function(params.arm=NULL,data=NULL) {
  i<-sum(data$ChosenArm==params.arm)
  return(i)
}

'
#See ZBanditGridProblem. This is an extension with the functionality 
to conduct a search multiple times in the same arm. Also adds functionality
to count total rewards/bugs/houses available and found. 

@params block.size: number of searches conducted each time an arm is pulled
'
ZBanditGridProblemBlockedResults <- function(test.time=NULL, params=NULL, 
                                              params.arm=NULL, new.infestations=NULL, block.size=NULL,
                                              params.bandit=NULL) {
  

  n_arms <- length(params.arm)
  
  search_stats <- list()
  
  #GENERATE INFESTATIONS
  #generates new infestations for each arm of the bandit
  if (new.infestations=="Yes") {
    infestations<- lapply(params.arm,GenerateInfestation,params=params) 
    assign("static.infestations",infestations,envir=.GlobalEnv)
  }
  
  if (! exists("infestations")) {
    infestations <- static.infestations
    #print("newfromold")
  }
  
  log10plus <- function(x) {
    log10(1+x)
  }
    
  #infestation stats
  total.houses.available <-params$nrow*params$ncol*Reduce("+",params.arm) #total infested for benchmarking
  total.bugs.available <- Reduce("+",lapply(infestations,sum)) #total bugs on chart
  total.rewards.available <- Reduce("+",lapply(lapply(infestations,log10plus),sum)) #total rewards if reward=log10(1+bugs)
  
  #SET UP BANDIT ARMS; CURRENTLY USING RING SEARCH STRATEGY
  #initialize RingSearch function for all infestations 
  search_stats <- lapply(infestations, RingSearch, st=NULL, max_actions=0, params=params)
  
  #function to pull bandit
  #RingSearch searches "block.size" on each pull (vs. 1 in ZBanditGridProblem)

  pull_arm <- function(chosen_arm, search_stats,block.size) {
    st     <- search_stats[[chosen_arm]]
    new_st <- RingSearch(data.frame(infestations[[chosen_arm]]),st=st,max_actions=block.size,params=params) 
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  
  #function to find reward for bandit
  #UPDATED: calculates rewards for a vector of results [number of searches]
  #UPDATED: rewards are log10(1+number of bugs found)
  
  search_reward <- function(new_st, block.size) {
    last <- tail(new_st$running_stats$total_bugs, (block.size+1)) #using total_bugs instead of total_found 
    last.reward <- rep(0,block.size)
    for (k in 2:(block.size+1)) {
      j=k-1
      last.reward[j] <- log10(1+last[k]-last[j])
    }
    return(last.reward)
  }
  
  #new function to determine total bugs found and total infested houses found in this turn
  BugsFound <-function(new_st, block.size) {
    
    last.bug      <- tail(new_st$running_stats$total_bugs,(block.size+1))
    last.house    <- tail(new_st$running_stats$total_found,(block.size+1))
    last.bug      <-last.bug[block.size+1]-last.bug[1]
    last.house    <-last.house[block.size+1]-last.house[1]
    bugs.found    <- (cbind(last.bug,last.house))
    return(bugs.found)
  }
  
  #function to find remaining prevalence (for testing)
  UnfoundPrevalence <- function(new_st) {
    unfound <-tail(new_st$running_stats$unfoundprevalence,1)
    return(unfound)
  }
  print("Stochastic search with Reinforcement-Comparison (rc) bandit");
  times       <- seq(1,test.time)
  arms        <- rep(0,length(times))
  rewards     <- rep(0,length(times))
  mR          <- rep(0,length(times))
  unfoundprev <- rep(0,length(times))
  blocking <- rep(0,length(times))
  ps_holding <- matrix(nrow=length(times),ncol=length(params.arm))
  bugs.houses <-matrix(nrow=length(times),ncol=2)
  static.houses.available <-rep(total.houses.available,length(times))
  static.bugs.available <-rep(total.bugs.available,length(times))
  static.rewards.available <-rep(total.rewards.available,length(times))
  
  bandit <- initialize_rc(n_arms=n_arms, learning_rate=params.bandit$learning, discount_factor=params.bandit$discount)
  print(bandit)
  
  for (trial in times) {
    ps <- NULL
    for (arm_idx in seq(n_arms)) {
      ps <- c(ps, probability_arm(bandit, arm_idx))
    }
    cat("trial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm <- next_arm(bandit)
    #print(chosen_arm)
   
    search_stats <- pull_arm(chosen_arm, search_stats,block.size)
    reward       <- search_reward(search_stats[[chosen_arm]],block.size) 
    unfound      <- UnfoundPrevalence(search_stats[[chosen_arm]])
    bug         <- BugsFound(search_stats[[chosen_arm]],block.size)
    
    cat("  arm: ", chosen_arm, "  reward", reward, " unfound", unfound, "total bugs&houses", bug, "\n")
    
    #Update the bandit for each time an arm was searched 
    #FIXME [SN]: NEED TO DECIDE THE ORDER IN WHICH TO PRESENT THE RESULT VECTOR TO THE BANDIT (RIGHT NOW IN ORDER OF SEARCH)
    for (block in block.size) {
      print(reward[block])
      bandit <- update_bandit_rc(bandit, chosen_arm, reward[block])
    }
    cat("  preferences:", paste(bandit$preferences), "\n")
    
    #Bandit reward is updated according to RC algorithm
    #rewards reported are the average reward in the set of pulls
    rewards[trial] <- sum(reward)/block.size
    arms[trial]    <- chosen_arm
    ps_holding[trial,] <-ps
    mR[trial] <- bandit$mean_reward
    unfoundprev[trial] <- unfound
    blocking[trial] <-block.size
    bugs.houses[trial,] <-bug
  }
  colnames(ps_holding)=params.arm
  colnames(bugs.houses)=c("total.bugs.found","total.infested.houses.found")

  results <- data.frame(ps_holding,bugs.houses,
                        T=times, ChosenArm=arms, BlockAvgReward=rewards, CumulativeAvgReward=cumsum(rewards*blocking), 
                        MeanReward=mR, UnfoundPrev.ChosenArm=unfoundprev, BlockSize=blocking,
                        static.rewards.available,static.bugs.available,static.houses.available)
  results$CumulativeHousesFound<-cumsum(results$total.infested.houses.found)
  results$CumulativeBugsFound<-cumsum(results$total.bugs.found)
  #write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)  
}

'run simulations on the bandit'
ZBanditSimulations <- function(NumSim=NULL, params.arm=NULL, params=NULL, 
                               new.infestations=NULL, test.time=NULL, graph_prev=NULL, params.bandit=NULL) {
  
  for (i in 1:NumSim) {
    if (i==1) {
      test <- ZBanditGridProblem(test.time=test.time, params=params, params.arm=params.arm,new.infestations="Yes", params.bandit=params.bandit)
      test.results <- lapply(seq_along(params.arm), GetResults, data=test)
      pl <- list()
      pl_prev <- list()
    } else {
      test <- ZBanditGridProblem(test.time=test.time, params=params, params.arm=params.arm,new.infestations=new.infestations, params.bandit=params.bandit)
      test.results.add <- lapply(seq_along(params.arm), GetResults, data=test)
      test.results <- rbind(test.results, test.results.add)
    }
    
    list <- colnames(test[,1:length(params.arm)])
    
    for (x in 1:length(params.arm)) {
      list[x]<-paste("Pref","Arm", x, list[x], sep="_")
    }
    
    test2 <- reshape(test,varying=1:length(params.arm), v.names="Pref",timevar="Arms",times=list, direction="long")
    test$ChosenLabel <- as.factor(paste("Prev","Arm", test$ChosenArm, as.numeric(params.arm[test$ChosenArm]), sep="_"))
    test$ChosenArm <- as.factor(test$ChosenArm)
    pl[[i]] <- ggplot() 
    pl[[i]] <- pl[[i]]+ geom_line(data=test2, aes(x=T, y=Pref,group=Arms,color=Arms)) + xlab("Time") + ylab("")
    pl[[i]] <- pl[[i]] + geom_point (data=test,aes(x=T, y=MeanReward,color="Mean Reward"))
    if (graph_prev=="Yes") {
      pl[[i]] <- pl[[i]] + geom_line (data=test, aes(x=T, y=UnfoundPrev.ChosenArm, group=ChosenArm, color=ChosenLabel))
    } else {
      pl_prev[[i]] <- ggplot()
      pl_prev[[i]] <- pl_prev[[i]] + geom_line (data=test, aes(x=T, y=UnfoundPrev.ChosenArm, group=ChosenArm, color=ChosenLabel))+ xlab("Time") +ylab("Unfound Prevalence")
    }  
    #pl[[i]] <- print(pl[[i]])
    #   pl <- pl + geom_point(data=test, aes(x=times, y=wins.by.A.timeline), color="black", size=0.5) + xlab("Time")
    #   pl <- pl + guides(color=FALSE, fill=FALSE)
  }
  colnames(test.results) <- params.arm
  rownames(test.results) <-NULL
  timestamp <- TimeNow()
  write.csv(test.results, paste("output/bandit_arm_results_", timestamp, ".csv", sep=""))
  pdf(paste("output/bandit_arm_plots_",timestamp,".pdf",sep=""))
  bquiet = lapply(pl, print)
  if (graph_prev !="Yes") {
    bquiet2 = lapply(pl_prev, print)
  }
  dev.off()
  return(test)
}

'run simulations on the blocked bandit'
ZBanditBlockedSimulations <- function(NumSim=NULL, params.arm=NULL, params=NULL, 
                               new.infestations=NULL, test.time=NULL, graph_prev=NULL, block.size=NULL, params.bandit=NULL) {
  
  for (i in 1:NumSim) {
    if (i==1) {
      test <- ZBanditGridProblemBlockedResults(test.time=test.time, params=params, params.arm=params.arm,new.infestations="Yes",block.size=block.size, params.bandit=params.bandit)
      test.results <- lapply(seq_along(params.arm), GetResults, data=test)
      test.info <- tail(test,1)
      test.info.df <-data.frame(T=test.info$T,BlockSize=test.info$BlockSize,
                          CumulativeAvgReward=test.info$CumulativeAvgReward,
                          CumulativeBugsFound=test.info$CumulativeBugsFound,
                          CumulativeHousesFound=test.info$CumulativeHousesFound,
                          Static.Rewards.Available=test.info$static.rewards.available,
                          Static.Bugs.Available=test.info$static.bugs.available,
                          Static.Houses.Available=test.info$static.houses.available
                          )
      pl <- list()
      pl_prev <- list()
    } else {
      test <-ZBanditGridProblemBlockedResults(test.time=test.time, params=params, params.arm=params.arm,new.infestations=new.infestations,block.size=block.size, params.bandit=params.bandit) 
      test.results.add <- lapply(seq_along(params.arm), GetResults, data=test)
      test.info.add <- tail(test,1)
      test.info.add.df <-data.frame( T=test.info.add$T,BlockSize=test.info.add$BlockSize,
                                    CumulativeAvgReward=test.info.add$CumulativeAvgReward,
                                    CumulativeBugsFound=test.info.add$CumulativeBugsFound,
                                    CumulativeHousesFound=test.info.add$CumulativeHousesFound,
                                    Static.Rewards.Available=test.info.add$static.rewards.available,
                                    Static.Bugs.Available=test.info.add$static.bugs.available,
                                    Static.Houses.Available=test.info.add$static.houses.available
                                    )   
      test.results <- rbind(test.results, test.results.add)
      test.info.df <-rbind(test.info.df,test.info.add.df)
      print(dim(test.info.df))
    }
    
    list <- colnames(test[,1:length(params.arm)])
    
    for (x in 1:length(params.arm)) {
      list[x]<-paste("Pref","Arm", x, list[x], sep="_")
    }
    
    test2 <- reshape(test,varying=1:length(params.arm), v.names="Pref",timevar="Arms",times=list, direction="long")
    test$ChosenLabel <- as.factor(paste("Prev","Arm", test$ChosenArm, as.numeric(params.arm[test$ChosenArm]), sep="_"))
    test$ChosenArm <- as.factor(test$ChosenArm)
    pl[[i]] <- ggplot() 
    pl[[i]] <- pl[[i]]+ geom_line(data=test2, aes(x=T, y=Pref,group=Arms,color=Arms)) + xlab("Time") + ylab("")
    pl[[i]] <- pl[[i]] + geom_point (data=test,aes(x=T, y=MeanReward,color="Mean Reward"))
    if (graph_prev=="Yes") {
      pl[[i]] <- pl[[i]] + geom_line (data=test, aes(x=T, y=UnfoundPrev.ChosenArm, group=ChosenArm, color=ChosenLabel))
    } else {
      pl_prev[[i]] <- ggplot()
      pl_prev[[i]] <- pl_prev[[i]] + geom_line(data=test, aes(x=T, y=UnfoundPrev.ChosenArm, group=ChosenArm, color=ChosenLabel))+ xlab("Time") +ylab("Unfound Prevalence")
    }  
    #pl[[i]] <- print(pl[[i]])
    #   pl <- pl + geom_point(data=test, aes(x=times, y=wins.by.A.timeline), color="black", size=0.5) + xlab("Time")
    #   pl <- pl + guides(color=FALSE, fill=FALSE)
  }
  colnames(test.results) <-params.arm
  rownames(test.results) <-NULL
  test.results <- cbind(test.results,test.info.df)
  test.results2 <-data.frame(lapply(test.results, as.numeric), stringsAsFactors=FALSE)
  timestamp <- TimeNow()
  write.csv(test.results2, paste("output/bandit_arm_results_", timestamp, ".csv", sep=""))
  pdf(paste("output/bandit_arm_plots_",timestamp,".pdf",sep=""))
  bquiet = lapply(pl, print)
  if (graph_prev !="Yes") {
    bquiet2 = lapply(pl_prev, print)
  }
  dev.off()
  return(test.results2)
}

# #RUN SIMULATIONS 
# # NOTE:   YOU MUST RUN EACH SIMULATION MORE THAN ONCE OR YOU WILL GET AN ERROR DUE TO
# #         FAILING TO MEET IN THE MINIMUM DIMENSIONS IN THE OUTPUT MATRIX       
# ZBanditSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.1,new.infestations="No", graph_prev="Yes")
# ZBanditSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.2,new.infestations="No", graph_prev="Yes")
# ZBanditSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.3,new.infestations="No", graph_prev="Yes")
# ZBanditSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.4,new.infestations="No", graph_prev="Yes")
# ZBanditSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No")
# ZBanditSimulations(NumSim=10,test.time=1000, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No")
# 
# ZBanditBlockedSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.1,new.infestations="No", graph_prev="Yes", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=50, params=params_grid_sim, params.arm=params.arm.2,new.infestations="No", graph_prev="Yes", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=50, params=params_grid_sim, params.arm=params.arm.3,new.infestations="No", graph_prev="Yes", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=400, params=params_grid_sim, params.arm=params.arm.4,new.infestations="No", graph_prev="Yes", block.size=10)
# 
# ZBanditBlockedSimulations(NumSim=10,test.time=100, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=20, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No", block.size=50)
# ZBanditBlockedSimulations(NumSim=10,test.time=30, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No", block.size=50)
# ZBanditBlockedSimulations(NumSim=10,test.time=50, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No", block.size=50)
# 
# 
# ZBanditSimulations(NumSim=10,test.time=1000, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No")
# ZBanditBlockedSimulations(NumSim=10,test.time=100, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=200, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=20, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No", block.size=50)
# ZBanditBlockedSimulations(NumSim=10,test.time=40, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No", block.size=50)
# 
# #Head to head size comparisons
# ZBanditSimulations(NumSim=10,test.time=1000, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No")
# ZBanditSimulations(NumSim=10,test.time=1000, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No")
# ZBanditBlockedSimulations(NumSim=10,test.time=100, params=params_grid_sim, params.arm=params.arm.5,new.infestations="No", graph_prev="No", block.size=10)
# ZBanditBlockedSimulations(NumSim=10,test.time=100, params=params_grid_sim2, params.arm=params.arm.6,new.infestations="No", graph_prev="No", block.size=10)
# 

ZBanditBlockedSimulations(NumSim=10,test.time=50, params=params_grid_sim2, params.arm=params.arm.7, params.bandit=params.bandit.2, new.infestations="No", graph_prev="Yes", block.size=10)

