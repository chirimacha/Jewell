#
#
DOCUMENTATION <- '

1. def_params below contains configuration of the software

2. use command line parameters to change the parameters at will.  
e.g. 
Rscript run_bandit.R --exec_mode="initialize"  --arm_names=TYABAYA,MELGAR         
Rscript run_bandit.R --exec_mode="update"      --arm=TYABAYA --rewards=0,0,0,1,0 --last_state_path="data/bandit_2015-07-02.csv" 

3. Notes:
* arm_names and rewards should be listed WITHOUT SPACES
* you can revise the same file by controlling the output_path:
Rscript run_bandit.R --exec_mode="initialize"  --arm_names="TYABAYA","MELGAR" --output_path=output/my_bandit.Rsave
Rscript run_bandit.R --exec_mode="update"      --rewards=0,0,0,1,10,0  --arm=TYABAYA  --last_state_path=output/my_bandit.Rsave --output_path=output/my_bandit.Rsave


Details with -h flag
'
TODO <- '
'

library(getopt)

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

if(! file.exists("output")) {
  dir.create("output")  
}
source("bandit.R")

#options(warn=2)

bandit_initialize <- function(params) {
  bandit <- initialize_rc(n_arms=length(params$arm_names))  
  #TODO: adjust the prior probabilities
  #TODO: adjust the parameters
  arm_names      <- params$arm_names  
  
  scores <- data.frame(site_num=integer(), arm=integer(), arm_name=character(), reward=double(), total_reward=double())
  
  recommended_arm <- next_arm_rc(bandit)
  cat("\n")
  cat("RECOMMEND ARM:", arm_names[recommended_arm], "\n")
  cat("\n")
  
  bandit_state <- list(bandit=bandit,
                       arm_names=arm_names,
                       recommended_arms=c(recommended_arm),
                       scores=scores)
  
  if(is.null(params[["output_path"]])) {
    fpath <- time_stamp(paste("output/bandit_state_", sep=""), ".RSave")
  } else {
    fpath <- params[["output_path"]]
  }
  save(bandit_state, file=fpath)
  cat("Saved bandit to:", fpath, "\n")
}

bandit_update <- function(params) {
  #load the bandit
  load(params[["last_state_path"]], newenv<-new.env())
  bandit_state <- get("bandit_state", newenv)   

  bandit <- bandit_state[["bandit"]]
  recommended_arms <- bandit_state[["recommended_arms"]]
  arm_names <- bandit_state[["arm_names"]]
  
  arm_idx <- which(arm_names == params$arm)
  scores  <- bandit_state[["scores"]]
  num_sites   <- length(params$rewards)
  
  site_num <- dim(scores)[1] + 1
  if(site_num == 1) {
    total_reward <- 0
  } else {
    total_reward <- scores[(site_num-1), "total_reward"]
  }
  for (trial in 1:num_sites) {
    cat("trial: ", trial, "\n")
    reward     <- as.numeric(params$rewards[trial])
    cat("  arm: ", params$arm, "  reward", reward, "\n")
    bandit <- update_bandit(bandit, arm_idx, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")

    total_reward <- total_reward + reward
    scores <- rbind(scores, list(site_num=site_num, arm=arm_idx, arm_name=params$arm, reward=reward, total_reward=total_reward))
    site_num <- site_num + 1
  }
  
  print(scores)
  
  recommended_arm <- next_arm_rc(bandit)
  cat("\n")
  cat("RECOMMEND ARM:", arm_names[recommended_arm], "\n")
  cat("\n")
  
  
  bandit_state <- list(bandit=bandit,
                       arm_names=arm_names,
                       recommended_arms=c(recommended_arms, recommended_arm),
                       scores=scores)
  if(is.null(params[["output_path"]])) {
    fpath <- time_stamp(paste("output/bandit_state_", sep=""), ".RSave")
  } else {
    fpath <- params[["output_path"]]
  }
  save(bandit_state, file=fpath)
  cat("Saved bandit to:", fpath, "\n")
}


#Rscript run_bandit.R --exec_mode="update"  --rewards=0,0,0,1  --arm=TYABAYA  --last_state_path=output/bandit_state_2015-07-07__21-51-37256.RSave 


parse_cmdl <- function(params=def_params, alt_params=list()) {
  cat("-------------------------------------------------------------------\n")
  cat("-------------- Multi-armed Bandit for Stochastic Search -----------\n")
  cat("-------------------------------------------------------------------\n")
  cat("\n")
  cat("Use -h to see options\n")
  for (p in names(alt_params)) {
    stopifnot(p %in% names(def_params))
    params[p] <- alt_params[p]
  }  
  spec <- matrix(c(
    'exec_mode',       'm', 1, "character", "One or more from [initialize,update]",
    'arm',             'p', 1, "character", " (update only) Name of arm pulled",
    'arm_names',       'a', 1, "character", " (initialization only) Two or more names of arms like c(\"name1\",\"name2\")",
    'rewards',         'r', 1, "character", " (update only) Numerical rewards separated by commas (0=nothing found, 1=infestation)",
    'output_path',     'o', 1, "character", "path for output of result",
    'last_state_path', 'l', 1, "character",  "path to last state of the bandit",
    'verbose',         'v', 2, "integer",   "NOT USED",
    'help'   ,         'h', 0, "logical",   "Writes this message."
  ), byrow=TRUE, ncol=5)
  opt <- getopt(spec)
  
  if ( !is.null(opt$exec_mode) ) {
    params[['exec_mode']] <- opt$exec_mode
  } else {
    params[['exec_mode']] <- def_params$exec_mode
  }
  cat("Exec mode:", params[['exec_mode']], "\n")

  if (params$exec_mode == "initialize") {
    if ( !is.null(opt$arm_names)) {
      params[['arm_names']] <- unlist(strsplit(opt$arm_names, ",")[[1]])
    } else {
      print("Warning: Using default arm names")
      params[['arm_names']] <- def_params$arm_names
    }
    cat("Arm names:", params[['arm_names']], "\n")
  }
  
  if (params$exec_mode == "update") {
    if ( !is.null(opt$rewards)) {
      params[['rewards']] <- unlist(strsplit(opt$rewards, ",")[[1]])
      cat("rewards:", params$rewards, "\n")      
    } else {
      print("Warning: missing rewards argument")
      q(1)
    }
    if ( !is.null(opt$arm)) {
      params[['arm']] <- opt$arm
      cat("Arm:", params$arm, "\n")      
    } else {
      print("Warning: missing arm argument")
      q(1)
    }
  }
  
  
  if ( !is.null(opt$last_state_path) ) {
    params[['last_state_path']] <- opt$last_state_path
    cat("State file:", opt$last_state_path, "\n")
    if(! file.exists(opt$last_state_path)) {
      cat("last_state_path:", opt$last_state_path, " does not exist!\n")
      q(1)
    }
  }
  if ( !is.null(opt$output_path) ) {
    params[['output_path']] <- opt$output_path
    cat("Output path:", opt$output_path, "\n")
    if(file.exists(opt$output_path)) {
      cat("Warning: output file:", opt$output_path, "already exists, and will be overwritten!\n")
    }
  }
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  }
  return(params)
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


############################################
#             MAIN PROGRAM
############################################
def_params <- list(  #overriden by some command line arguments
  exec_mode = "initialize",  #update
  arm_names = c("arm1", "arm2"),
  some_other_param=NULL
)

browser()
params <- parse_cmdl(params=def_params, alt_params=list())

if(params$exec_mode == "initialize") {
  bandit_initialize(params=params)
} else if(params$exec_mode == "update") {
  bandit_update(params=params)
} else {
  print("Unrecognized execution mode:")
  print(params$exec_mode)
}
