#####################################################################

#implementation of bandit heuristics

#1. Reinforcement comparison methods Sutton and Barto (1998) 
#2. UCB family of algorithms by Auer, Cesa-Bianchi & Fisher (2002)

# version 2014-10-25 by Sasha Gutfraind 
#   based on some testing code from White et al.
# ref: "Algorithms for the multi-armed bandit problem" / Precup, Kuleshov

##################
#  usage guide
##################
#something like
#bandit <- initialize_[BANDIT_NAME](<args>) #must be called before any other function
#cycle #1
#arm_idx <- next_arm(bandit)
#r <- gaussian_reward_toy(my_system,arm_idx)
#bandit <- update_bandit(bandit,arm_idx,r)
#cycle #2
#arm_idx <- next_arm(bandit)
#r <- gaussian_reward_toy(system,arm_idx)
#bandit <- update_bandit(bandit,arm_idx,r)
#...etc

##################
#  note
##################
# arm is the same as strategy.  arms are numbered 1..n_arms

#####################################################################

library(reshape2)
library(plyr)
library(ggplot2)

#TODO: set this to your code, or create an environment variable
#setwd("~/Dropbox/chagas_models_aim3_spatial_uncertainty/code")
setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep=""))

timeNow <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")

#n_arms<-3
test_time <- 100

#the test reward function: the reward is stochastic, often 0, and when positive averages
#U(arm1) = 1
#U(arm2) = 2
#...
#true_best_arm = n_arms

gaussian_reward_toy <- function(arm_idx) { 
  stopifnot(arm_idx >= 1)
  stopifnot(arm_idx <= n_arms)
  if(runif(1) < 0.5) {
    return(0)
  } else {
    return(rnorm(1, mean=arm_idx));
  }
  return(rnorm(1, mean=arm_idx));
}

initialize_rc <- function (n_arms, learning_rate=0.1, discount_factor=0.1) {
#Reinforcement comparison bandit
#learning_rate=in [0,1] increases the tendency to prefer recently successful strategies
#discount_factor in [0,1];  greater for more restless bandits
  bandit <- list("name"="rc", "n_arms"=n_arms, learning_rate=learning_rate, discount_factor=discount_factor)
  bandit$preferences <- rep(0, n_arms)
  bandit$mean_reward <- 0.
  #print(bandit)
  return(bandit)
}

initialize_ucb1 <- function(n_arms, reward_func) {
  #Upper Confidence Bound (1) bandit
  bandit <- list("name"="ucb1", "n_arms"=n_arms)
  bandit$mean_reward <- rep(0, n_arms)
  for(arm_idx in seq(n_arms)) {
    bandit$mean_reward[arm_idx] <- reward_func(arm_idx)
  }
  bandit$num_pulls   <- rep(1,n_arms)
  #print(bandit)
  return(bandit)
}

next_arm <- function(bandit) {
#returns the randomly-drawn next arm
  if(bandit$name=="rc") {
    return(next_arm_rc(bandit));
  } else {if (bandit$name=="ucb1") {
    return(next_arm_ucb1(bandit));
  } else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } }
}  

next_arm_rc <- function(bandit) {
  #returns the randomly-drawn next arm
  p_vals <- exp(bandit$preferences)
  return(sample.int(bandit$n_arms, size=1, replace=TRUE, prob=p_vals))
}

next_arm_ucb1 <- function(bandit) {
  #returns the next arm.  this is a deterministic algorithm
  total_pulls <- sum(bandit$num_pulls)
  values <- list()
  for(arm_idx in seq(bandit$n_arms)) {
    values <- c(values, bandit$mean_reward[arm_idx] + sqrt( (2*log(total_pulls))/bandit$num_pulls[arm_idx] ))  
  }
  best_arm <- which.max(values)  #in cases of ties, it returns the first arm
  #print(paste(values))
  #cat("best arm: ", best_arm, "\n")
  return(best_arm)
}

probability_arm <- function(bandit, arm_idx) {
#returns the current probability of selecting arm with index arm_idx
  if(bandit$name=="rc") {
    return(probability_arm_rc(bandit, arm_idx));
  } else {if (bandit$name=="ucb1") {
    return(probability_arm_ucb1(bandit, arm_idx));
  } else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } } 
}

probability_arm_rc <- function(bandit, arm_idx) {
#compute the probability for an arm
  return(exp(bandit$preferences[arm_idx]) / sum(exp(bandit$preferences)))
}

probability_arm_ucb1 <- function(bandit, arm_idx) {
  #compute the probability for an arm
  cat("WARNING: this is a deterministic algorithm\n")
  best_arm <- next_arm_ucb1(bandit)
  if(arm_idx == best_arm){
    return(1.0)
  } else {
    return (0.0)
  }
}

update_bandit <- function(bandit, chosen_arm, reward) {
#called after the reward was returned
  if(bandit$name=="rc") {
    return(update_bandit_rc(bandit, chosen_arm, reward))
  } else {if (bandit$name=="ucb1") {
    return(update_bandit_ucb1(bandit, chosen_arm, reward))
  } else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } } 
}  


update_bandit_rc <- function(bandit, chosen_arm, reward) {
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\t", bandit$preferences, "\n")
  #cat("ZZZ  preference: ", bandit$preferences[chosen_arm], "  arm: ", chosen_arm, "   reward: ", reward, "  mean_reward: ", bandit$mean_reward, "\n")
  bandit$preferences[chosen_arm] <- bandit$preferences[chosen_arm]                + bandit$learning_rate   * (reward - bandit$mean_reward)
  bandit$mean_reward             <- (1-bandit$discount_factor)*bandit$mean_reward + bandit$discount_factor * reward 
  #bandit$mean_reward <- bandit$mean_reward + bandit$discount_factor*(reward - bandit$mean_reward)
  #cat("ZZZ  preference: ", bandit$preferences[chosen_arm], "\n")
  return(bandit)
}

update_bandit_ucb1 <- function(bandit, chosen_arm, reward) {
  #cat("ZZZ:", paste(bandit$num_pulls), "\n")
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\n")
  #TODO: should the first term be multiplied by (1-bandit$learning_rate)
  np <- bandit$num_pulls[chosen_arm]
  bandit$mean_reward[chosen_arm] <- (np*bandit$mean_reward[chosen_arm] + reward)/(np+1)
  bandit$num_pulls[chosen_arm]   <- np + 1
  #cat("ZZZ:", paste(bandit$num_pulls), "\n")
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\n")
  return(bandit)
}

#############################################
#  PLOT and COMPARE
#############################################


zz_test_bandit_rc <- function(n_arms=3, test_time=100) {
  print("testing rc bandit");
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
    reward     <- gaussian_reward_toy(chosen_arm)
    cat("  arm: ", chosen_arm, "  reward", reward, "\n")
    bandit <- update_bandit(bandit, chosen_arm, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
  }
  results <- data.frame(T=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards))
  write.csv(results, paste("output/rc_results_", timeNow, ".csv", sep=""))
  return(results)
}

zz_test_bandit_ucb1 <- function(n_arms=3, test_time=100) {
  print("testing ucb1 bandit");
  times   <- seq(1,test_time)
  arms    <- rep(0,length(times))
  rewards <- rep(0,length(times))
  bandit <- initialize_ucb1(n_arms=n_arms, reward_func=gaussian_reward_toy)
  print(bandit)
  for (trial in times) {
    cat("trial: ", trial, "\n")
    chosen_arm <- next_arm(bandit)
    reward     <- gaussian_reward_toy(chosen_arm)
    cat("  arm: ", chosen_arm, "  reward", reward, "\n")
    bandit <- update_bandit(bandit, chosen_arm, reward)
    cat("  num_pulls:    ", paste(bandit$num_pulls), "\n")    
    cat("  mean_rewards: ", paste(bandit$mean_reward), "\n")    
    rewards[trial] <- reward
    arms[trial]    <- chosen_arm
  }
  results <- data.frame(T=times, ChosenArm=arms, Reward=rewards, CumulativeReward=cumsum(rewards))
  write.csv(results, paste("output/ucb1_results_", timeNow, ".csv", sep=""))
  return(results)
}

zz_comparison_plot <- function(true_best_arm, rc.results=NULL, ucb1.results=NULL) {
  if(class(rc.results) == class(NULL)) {
    rc.results <- read.csv("returns/rc_results_May-22-2013_11-37-20.csv", header = TRUE)
    #names(rc.results) <- c("Sim", "T", "ChosenArm", "Reward", "CumulativeReward")
  }
  rc.results$CumCorrect <- cumsum(rc.results$ChosenArm == true_best_arm)
  rc.results <- transform(rc.results, Algorithm = "ReinfComp")
  if(class(ucb1.results) == class(NULL)) {    
    ucb1.results <- read.csv("returns/ucb1_results_May-22-2013_11-37-20.csv", header = TRUE)
    #names(ucb1.results) <- c("Sim", "T", "ChosenArm", "Reward", "CumulativeReward")
  }
  ucb1.results$CumCorrect <- cumsum(ucb1.results$ChosenArm == true_best_arm)
  ucb1.results <- transform(ucb1.results, Algorithm = "UCB1")
  results <- rbind(rc.results, ucb1.results)
  # Plot average reward as a function of time.
  stats <- ddply(results,
                 c("Algorithm", "T"),
                 function (df) {df$CumulativeReward/df$T})
  #ylim(0, 1) +
  #return(results)
  pl<-ggplot(stats, aes(x = T, y = V1, group = Algorithm, color = Algorithm)) +
    geom_line() +
    xlab("Time") +
    ylab("Average Reward") +
    ggtitle("Performance of Different Algorithms")
  print(pl)
  ggsave("output/simple_comparisons_average_reward.pdf")

  # Plot frequency of selecting correct arm as a function of time.
  # In this instance, 5 is the correct arm.
  stats <- ddply(results,
                 c("Algorithm", "T"),
                 function (df) {df$CumCorrect/df$T})
  #ylim(0, 1) +
  #rownames(stats) <- stats[,"Algorithm"]
  #stats[,"Algorithm"] <- NULL
  #stats <- stats / seq(test_time)
  pl2<-ggplot(stats, aes(x = T, y = V1, group = Algorithm, color = Algorithm)) +
    geom_line() +
    xlab("Time") +
    ylab("Probability of Selecting Best Arm") +
    ggtitle("Accuracy of Different Algorithms")
  print(pl2)
  ggsave("output/simple_comparisons_average_accuracy.pdf")
  
  # Plot variance of chosen arms as a function of time.
  #stats <- ddply(results,
  #               c("Algorithm", "T"),
  #               function (df) {var(df$ChosenArm)})
  #pl3<-ggplot(stats, aes(x = T, y = V1, group = Algorithm, color = Algorithm)) +
  #  geom_line() +
  #  xlab("Time") +
  #  ylab("Variance of Chosen Arm") +
  #  ggtitle("Variability of Different Algorithms")
  #print(pl3)
  #ggsave("output/simple_comparisons_variance_choices.pdf")
  
  # Plot cumulative reward as a function of time.
  stats <- ddply(results,
                 c("Algorithm", "T"),
                 function (df) {mean(df$CumulativeReward)})
  pl4<-ggplot(stats, aes(x = T, y = V1, group = Algorithm, color = Algorithm)) +
    geom_line() +
    xlab("Time") +
    ylab("Cumulative Reward of Chosen Arm") +
    ggtitle("Cumulative Reward of Different Algorithms")
  print(pl4)
  ggsave("output/simple_comparisons_cumulative_reward.pdf")
}

z_demo <- function() {
  rc.results   <- zz_test_bandit_rc(n_arms=3, test_time=100)
  ucb1.results <- zz_test_bandit_ucb1(n_arms=3, test_time=100)
  browserxyz<-zz_comparison_plot(true_best_arm=3, rc.results=rc.results, ucb1.results=ucb1.results)
}

