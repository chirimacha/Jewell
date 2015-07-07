#
#
DOCUMENTATION <- '

1. def_params below contains configuration of the software

2. use command line parameters to change the parameters at will.  
e.g. 
Rscript run_bandit.R --exec_mode="initialize"  --arms=2   --OTHER_PARAM=TODO_PARAM_DATA 
Rscript run_bandit.R --exec_mode="update"      --sites=20 --hits=10 --last_state_path="data/bandit_2015-07-02.csv" 

Details with -h flag
'
TODO <- '
'

library(getopt)
setwd(paste(Sys.getenv("UPTAKE"), "/northshore/wisca_apm/code", sep=""))

#options(warn=2)

#source("paths.R")
source("bandit.R")

bandit_initialize <- function(params) {
  bandit <- initialize_rc(arms=params$arms)  #TODO: actually call the bandit initialize
  
  #output_path
}

bandit_update <- function(params) {
  #load the bandit
  
  #
  
  #output_path
}

parse_cmdl <- function(params=def_params, alt_params=list()) {
  cat("-------------------------------------------------------------------\n")
  cat("-------------- Antibiotic Prescription Manager prototype ----------\n")
  #cat("-------- (c) Upt^ke Technologies 2015 + NorthShore Health----------\n")
  cat("-------------------------------------------------------------------\n")
  cat("\n")
  cat("Use -h to see options\n")
  for (p in names(alt_params)) {
    stopifnot(p %in% names(def_params))
    params[p] <- alt_params[p]
  }  
  spec <- matrix(c(
    'exec_mode',       'm', 1, "character", "One or more from [initialize,update]",
    'output_path',     'o', 1, "character", "path for output of WISCA or pWISCA result",
    'last_state_path', 'p', 1, "character",  "path to last state of the bandit",
    'verbose',         'v', 2, "integer",   "NOT USED",
    'help'   ,         'h', 0, "logical",   "Writes this message."
  ), byrow=TRUE, ncol=5)
  opt <- getopt(spec)
  
  if ( !is.null(opt$exec_mode) ) {
    params[['exec_mode']] <- opt$exec_mode
    cat("Exec mode:", opt$exec_mode, "\n")
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
  arms      = 2,
  some_param=NULL
)

params <- z_parse_cmdl(params=def_params, alt_params=list())

if(params$exec_mode == "initialize") {
  bandit_initialize(params=params)
} else if(params$exec_mode == "initialize") {
  bandit_update(params=params)
} 
