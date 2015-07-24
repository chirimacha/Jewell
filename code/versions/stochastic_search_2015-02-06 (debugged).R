#this code is for some new algorithms to simulate stochastic search
library(reshape2)
library(ggplot2)

setwd("~/Dropbox/chagas_models_aim3_spatial_uncertainty/code")

source("bandit.R")

params <- list(
  max_cost=100,
  ncol=20,
  nrow=20,
  num_starts=2,
  p_hop=0.6,
  p_skip=0.3,
  p_jump=0.1,
  pMovePerDay=0.01, #per site 
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  num_days=365
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
  print("infestation")
  print(infested_coords)
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
  pl <- ggplot(data=visited_m, aes(lat,lon)) + geom_raster(data=visited_m, aes(fill=visited, alpha=1.0)) + geom_point(data=infest_m, aes(lat,lon))
  #pl <- pl + geom_point(data=all_data, aes(shape=factor(visited)))
  #pl <- pl + geom_point(data=visited_m, aes(shape=factor(visited)))
  #pl <- pl + 
  
  print(pl)
}

random_ring <- function(infestation, params) {
  visited <- matrix(0, nrow=dim(infestation)[1], ncol=dim(infestation)[2])
  
  get_unvisited_nbs <- function(house, visited) {
    nbs <- read.csv(text="lon,lat")
    for(x in c(house$lon-1, house$lon, house$lon+1)) {
      if (x < 1 | x > params$ncol) {
        next;
      }
      for(y in c(house$lat-1, house$lat, house$lat+1)) {
          if (y < 1 | y > params$nrow) {
            next;
          }
          nb = list(lon=x,lat=y)
          if (all(nb == house)) {
            next;
          }
          if (visited[nb$lat,nb$lon] > 0) {
            next;
          }
          nbs <- rbind(nbs, nb)
        }
      }
    return(nbs)
  }
  #total_found <- 0
  #total_visited <- 0
  #total_cost <- 0
  #step <- 0
  total_squares <- dim(infestation)[1] * dim(infestation)[2]
  randomized_sites <- melt(visited)
  randomized_sites <- randomized_sites[sample(total_squares),]
  names(randomized_sites)<-c("lat", "lon", "val")
  randomized_sites$val <- NULL
  continue_inspection <- params[["continue_inspection"]]
  if(is.null(continue_inspection)) {
    continue_inspection <- function(infestation, latest) {
      return(latest$total_cost < params$max_cost & dim(latest)[1] < total_squares);
    }
  }
  known_infested <- read.csv(text="lon,lat")
  suspected_nbs  <- read.csv(text="lon,lat")
  next_suspected_nb <- 1
  running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0))
  next_stat <- tail(running_stats, 1)
  next_random_idx <- 1
  #browser()
  while(next_stat$total_visited < total_squares & continue_inspection(infestation, next_stat)) {
    #next_stat$step <- next_stat$step + 1
    next_site <- NULL
    #if(next_suspected_nb < dim(suspected_nbs)[1]){
    #  browser()
    #}      
    while (next_suspected_nb < dim(suspected_nbs)[1]) {
      next_nb <- suspected_nbs[next_suspected_nb,]
      next_suspected_nb <- next_suspected_nb + 1
      if(visited[next_nb$lat, next_nb$lon] == 0) {
        next_site <- next_nb
        break
      }
    }
    while (is.null(next_site) & next_random_idx < total_squares) {
      next_site <- randomized_sites[next_random_idx,]
      if (visited[next_site$lat, next_site$lon] == 0) {
        next_random_idx <- next_random_idx + 1
        break
      }
      next_random_idx <- next_random_idx + 1
    }
    if(is.null(next_site)) {
      break
    }
    next_stat$total_cost    <- next_stat$total_cost + 1
    next_stat$total_visited <- next_stat$total_visited + 1
    visited[next_site$lat, next_site$lon] <- 1
    if(infestation[next_site$lat, next_site$lon] > 0) {
      visited[next_site$lat, next_site$lon] <- 2
      known_infested <- rbind(known_infested, next_site)
      next_stat$total_found <- next_stat$total_found + 1
      neighbors <- get_unvisited_nbs(next_site, visited)
      suspected_nbs <- rbind(suspected_nbs, neighbors)
      sprintf("1, lon=%d, lat=%d", next_site$lon, next_site$lat)
    } else {
      print(0)
    }
    running_stats <- rbind(running_stats, next_stat)
    next_site <- NULL
  }
  running_stats$infestation_lower_bound <- running_stats$total_found/total_squares
  running_stats$infestation_estimate    <- running_stats$total_found/running_stats$total_visited
  print(tail(running_stats,1))
  return(list(running_stats=running_stats,visited=visited))
}

infestation_test_m <- generate_infestation(params)
search_stats       <- random_ring(infestation_test_m, params)
visited            <- search_stats$visited
search_stats$running_stats
plot_grid(infestation_test_m, visited=visited)
#plot_grid(infestation_test_m, visited=NULL)