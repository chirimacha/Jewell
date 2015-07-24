'
################################
Bandit Simulations for Optimization

This code:
  -Calls Bandit Code 
  -Simulates infestations
  -Simulates Searches
  -Simulates Bandit based on the parameters of the bandit
  -Returns "Rewards" sum(log(1-bugs))


To call this code you need to:

Set prevalence/arm parameters

######################################
'

library(reshape2)
library(ggplot2)
library(pROC)
library(sm)

# number <- sample(1:10000, 1)
# print(number)
# set.seed(number)

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

#CALL BANDIT CODE 
source("Bandit.R")

#Infestation Generator
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

#Toy Search Function: Ring Search
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

'function for applying the bandit to a given infestation, 
using ring search and banditing on location 
@params test.time:        number of times to run the bandit
@params params:           grid parameters
@params params.arm:       arms and prevalence parameters
@params new.infestation:  generate new infestations or start with old ones [for benchmarking]


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
  
  bandit <- initialize_rc(n_arms=n_arms, learning_rate=params.bandit[1], discount_factor=params.bandit[2])
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

    randomizer <- runif(block.size) # IMPORTANT: we are randomizing the order in which chiris are presented
    random.chiris <-data.frame(cbind(reward,randomizer))
    random.chiris <- random.chiris[order(randomizer),]
    # print(random.chiris)
    
    for (block in block.size) {
      # print(random.chiris$reward[block])
      bandit <- update_bandit_rc(bandit, chosen_arm, random.chiris$reward[block])
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
  
  #Only care about the number of bugs
#   results <- data.frame(ps_holding,bugs.houses,
#                         T=times, ChosenArm=arms, BlockAvgReward=rewards, CumulativeAvgReward=cumsum(rewards*blocking), 
#                         MeanReward=mR, UnfoundPrev.ChosenArm=unfoundprev, BlockSize=blocking,
#                         static.rewards.available,static.bugs.available,static.houses.available)
#   results$CumulativeHousesFound<-cumsum(results$total.infested.houses.found)
#   results$CumulativeBugsFound<-cumsum(results$total.bugs.found)
  results<-sum(rewards*blocking)
  
  #write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)  
}

