#
#
options(warn=2)

DOCUMENTATION <- '

1. def_params below contains configuration of the software

2. use command line parameters to change the parameters at will.  
e.g. 
Rscript run_bandit.R --exec_mode="initialize"  --arm_names=TYABAYA,MELGAR         
Rscript run_bandit.R --exec_mode="update"      --arm=TYABAYA --chiris=0,0,0,1,0 --last_state_path="data/bandit_2015-07-02.csv" 

3. Notes:
* arm_names and rewards should be listed WITHOUT SPACES
* generally, last_state_path and output_path should be different files like
--last_state_path=Bandit2015_07_17.RSave
--output_path=Bandit2015_07_18.RSave
* if output_path is not specified, the file last_state_path would NOT be updated, instead a new file will be created.
* you can revise the same file by controlling the output_path:
Rscript run_bandit.R --exec_mode="initialize"  --arm_names="TYABAYA","MELGAR" --output_path=output/my_bandit.Rsave
Rscript run_bandit.R --exec_mode="update"      --chiris=0,0,0,1,2.1,0  --arm=TIABAYA  --last_state_path=output/my_bandit.Rsave --output_path=output/my_bandit.Rsave


Details with -h flag

Rpackages required:
getopt
reshape2
plyr
ggplot2

Español
1. def_params: configuración del codigo

2. se puede usar el "command line" para cambiar los parametros en cualquier momento.  
e.g. 
Rscript run_bandit.R --exec_mode="initialize"  --arm_names=TIABAYA,MELGAR         
Rscript run_bandit.R --exec_mode="update"      --arm=TIABAYA --chiris=0,0,0,1,0 --recomendaciones=1 --last_state_path="data/bandit_2015-07-02.csv" 

3. Notas:
* Entra los arm_names and chiris SIN ESPACIOS 
* "last_state_path" y "output_path" deben ser archivos diferentes. Por ejemplo: 
--last_state_path=Bandit2015_07_17.RSave
--output_path=Bandit2015_07_18.RSave
* Si no se entra un "output_path," el archivo last_state_path no cambiará y el programa hará archivo nuevo (con otro nombre). 
* Si quieres revisar un archivo puedes hacerlo con el output_path. Por ejemplo: 
Rscript run_bandit.R --exec_mode="initialize"  --arm_names=TIABAYA,"MELGAR" --output_path=output/my_bandit.Rsave
Rscript run_bandit.R --exec_mode="update"      --chiris=0,0,0,1,2.1,0  --arm=TIABAYA  --last_state_path=output/my_bandit.Rsave --output_path=output/my_bandit.Rsave


Hay más detalles con "-h" 

Se necesita estos Rpackages:
getopt
reshape2
plyr
ggplot2


## TODO <-

[DONE SG] Update to make rewards =log10(1+chiris)
[TODO SN] Create running record of things entered. Potentially think about pulling in files rather than typing
          in rewards. 
[TODO SN/SG] Optimize parameters 

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
  
  bandit.run <-1 #BANDIT COUNTER
  
  bandit <- initialize_rc(n_arms=length(params$arm_names))  
  #TODO: adjust the prior probabilities
  #TODO: adjust the parameters
  arm_names      <- params$arm_names  

  scores <- data.frame(site_num=integer(), arm=integer(), arm_name=character(), chiris=double(), reward=double(), total_reward=double())
  
  if (is.na(params$recomendaciones)) {
    params$recomendaciones=1
  }
  for (pulls in 1:params$recomendaciones) {
    recommended_arm <- next_arm_rc(bandit)
    cat("\n")
    cat("RECOMMEND ARM:", arm_names[recommended_arm], "\n")
    cat("\n")
    
    if (pulls==1) {
      all.recommended <-recommended_arm 
    } else {
      all.recommended <-c(all.recommended,recommended_arm)
    }
  }
  
  flip<-as.matrix(all.recommended)
  recs<-data.frame(arm.recommended=flip)
  recs$bandit.run <-bandit.run
  
  
  bandit_state <- list(bandit=bandit,
                       arm_names=arm_names,
                       recommended_arms=c(all.recommended),  
                       scores=scores,
                       bandit.run=bandit.run,
                       all.recs=recs,
                       recomendaciones=params$recomendaciones)
  
  if(is.null(params[["output_path"]])) {
    fpath <- time_stamp(paste("output/bandit_state_","run_",bandit.run,"_", sep=""), ".RSave")
  } else {
    fpath <- params[["output_path"]]
  }
  fpathcsv <-time_stamp(paste("output/bandit_state_","run_",bandit.run, sep=""), ".csv")
  save(bandit_state, file=fpath)
  
  cat("Saved bandit to:", fpath, "\n")
}

bandit_update <- function(params) {
  #load the bandit
  load(params[["last_state_path"]], newenv<-new.env())
  bandit_state <- get("bandit_state", newenv)
  bandit.run <-bandit_state[["bandit.run"]]+1 #Update counter
  
  bandit <- bandit_state[["bandit"]]
  recommended_arms <- bandit_state[["recommended_arms"]]
  arm_names <- bandit_state[["arm_names"]]
  all.recs <- bandit_state[["all.recs"]]

  arm_idx <- which(arm_names == params$arm)
  if(length(arm_idx)==0) {
    cat("Arm: ", params$arm, "is not recognized!\n.  Use one of:", arm_names)
    cat("Existing!\n")
    q()
  }
  scores  <- bandit_state[["scores"]]
  num_sites   <- length(params$chiris)
  
  randomizer <- runif(num_sites) # IMPORTANT: we are randomizing the order in which chiris are presented
  found.chiris <-params$chiris
  random.chiris <-data.frame(cbind(found.chiris,randomizer))
  #print(random.chiris)
  random.chiris <- random.chiris[order(randomizer),]
  #print(random.chiris)
  found.chiris <- as.vector(random.chiris$found.chiris)
  #print(found.chiris)
  
  site_num <- dim(scores)[1] + 1
  if(site_num == 1) {
    total_reward <- 0
  } else {
    total_reward <- scores[(site_num-1), "total_reward"]
  }
  for (trial in 1:num_sites) {
    cat("trial: ", trial, "\n")
    # chiris     <- as.numeric(params$chiris[trial])  
    chiris  <-as.numeric(found.chiris[trial]) #RANDOMIZED ORDER CHIRIS
    reward   <- log10(1+chiris) #IMPORTANT: we set rewards as log10(1+chiris)
    cat("  arm: ", params$arm, "  chiris", chiris, "  reward", reward, "\n")
    bandit <- update_bandit(bandit, arm_idx, reward)
    cat("  preferences:", paste(bandit$preferences), "\n")

    total_reward <- total_reward + reward
    scores <- data.frame(rbind(scores, list(bandit.run=bandit.run,site_num=site_num, arm=arm_idx, arm_name=params$arm, 
                                 chiris=chiris, reward=reward, total_reward=total_reward)))
    scores$arm_name <- as.character(scores$arm_name)  #works around initialization bug

    site_num <- site_num + 1
  }
  
  row.names(scores)<-NULL 
  print(scores)
  pref<-data.frame(bandit$preferences)
  pref<-t(pref)
  row.names(pref)<-NULL
  colnames(pref)<-arm_names
  scores.pref<-cbind(scores,pref)
  
  
  if (is.na(params$recomendaciones)) {
    params$recomendaciones=1
  }
  for (pulls in 1:params$recomendaciones) {
    recommended_arm <- next_arm_rc(bandit)
  
    cat("\n")
    cat("RECOMMEND ARM:", arm_names[recommended_arm], "\n")
    cat("\n")
  
    if (pulls==1) {
      all.recommended <-recommended_arm 
    } else {
      all.recommended <-c(all.recommended,recommended_arm)
    }
  }
  
  
  flip<-as.matrix(all.recommended)
  recs<-data.frame(arm.recommended=flip)
  recs$bandit.run <-bandit.run
  
  all.recs <-rbind(all.recs,recs)
  
  bandit_state <- list(bandit=bandit,
                       arm_names=arm_names,
                       recommended_arms=c(recommended_arms, all.recommended),
                       scores=scores, 
                       bandit.run=bandit.run,
                       all.recs=all.recs,
                       recomendaciones=params$recomendaciones)
  
    
  
    
#   max <- max(bandit_state$recommended_arms)
#   total_rec <-sum(recomendaciones) 
#   recs<-matrix(NA,nrow=bandit.run,ncol=max)
#   
#   for (i in 1:total_rec) {
#     if i=1 {
#       rec <-bandit.state$recomendaciones[i]
#       all.rec <-bandit.state$recomendaciones[i]
#     } else {
#       r
#     for (j in 1:rec) {
#       recs[i,j] <-       
  
  if(is.null(params[["output_path"]])) {
    fpath <- time_stamp(paste("output/bandit_state_","run_",bandit.run,"_", sep=""), ".RSave") 
  } else {
    fpath <- params[["output_path"]]
  }
  save(bandit_state, file=fpath)
  write.csv(scores.pref,paste("output/bandit_arm_results_","run_",bandit.run,"_", TimeNow(), ".csv", sep=""),row.names =FALSE)
  write.csv(all.recs,paste("output/bandit_arm_recs_","run_",bandit.run,"_", TimeNow(), ".csv", sep=""),row.names =FALSE)
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
    'arm',             'p', 1, "character", " (update only, manual only) Name of arm pulled",
    'arm_names',       'a', 1, "character", " (initialization only, manual only) Two or more names of arms like c(\"name1\",\"name2\")",
    'chiris',          'c', 1, "character", " (update only, manual only) Number of bugs (chiris) separated by commas (0=nothing found, >0 size of infestation (#bugs))",
    'output_path',     'o', 1, "character", "path for output of result",
    'last_state_path', 'l', 1, "character",  "path to last state of the bandit",
    'recomendaciones', 'r', 1, "integer",  "# de recommendaciones/calls to the bandit (optional)" ,
    'iniciar',         'i', 1, "character", '(initialization only, auto only) path to intialization data',
    'inspecciones_ruta','ins', 1, "character", '(update only, auto only) path to inspecciones data',
    'dia' ,            'd', 2, "integer", '(update only, auto only) day from which to pull data',
    'mes' ,            'mes', 2, "integer", '(update only, auto only) month from which to pull data',
    'ano' ,            'n', 4, "integer", '(update only, auto only) month from which to pull data',
    'verbose',         'v', 2, "integer",   "NOT USED",
    'help'   ,         'h', 0, "logical",   "Writes this message."
  ), byrow=TRUE, ncol=5)
  
  cmdArgs<-commandArgs(TRUE)
  #for debugging from Rstudio
  #cmdArgs<-c('--exec_mode=update', '--arm=MELGAR', '--chiris=0,0,1,0,0,0,0,0,0,0', 
  #           '--output_path=/tmp/1.rs','--last_state_path=/tmp/1.rs')
  #cmdArgs<-c('--exec_mode=initialize', '--arm_names=TIABAYA,ASA,MELGAR,HUNTER', 
  #           '--output_path=/tmp/1.rs')
  

  opt <- getopt(spec, opt=cmdArgs)
  
  if ( !is.null(opt$exec_mode) ) {
    params[['exec_mode']] <- opt$exec_mode
  } else {
    params[['exec_mode']] <- def_params$exec_mode
  }
  cat("Exec mode:", params[['exec_mode']], "\n")

  if (params$exec_mode == "initialize") {
    if ( !is.null(opt$iniciar)) {
     assignments <- read.csv(opt$iniciar)
     unique_arms <-as.list(as.character(unique(assignments$sector)))
     params[['arm_names']] <-unlist(unique_arms)
    } else {
      if ( !is.null(opt$arm_names)) {
        params[['arm_names']] <- unlist(strsplit(opt$arm_names, ",")[[1]])
      } else {
        print("Warning: Using default arm names")
        params[['arm_names']] <- def_params$arm_names
      }
    }
    cat("Arm names:", params[['arm_names']], "\n")
  }
  
  if (params$exec_mode == "update") {
    if (! is.null(opt$inspecciones_ruta)) {
      #auto mode 
      inspect <- read.csv(opt$inspecciones_ruta) #read data
      if (! is.null(opt$dia) & ! is.null(opt$ano) & !isnull
    
    if ( !is.null(opt$chiris)) {
      params[['chiris']] <- unlist(strsplit(opt$chiris, ",")[[1]])
      cat("chiris:", params$chiris, "\n")      
    } else {
      print("Warning: missing chiris argument")
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
  if ( !is.null(opt$recomendaciones)) {
    params[['recomendaciones']] <- opt$recomendaciones
    cat("recomendaciones:", opt$recomendaciones, "\n")
  } else {
    params[['recomendaciones']] <- def_params$recomendaciones
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

TimeNow <- function() {
  #generates the current time
  x <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")
  return(x)
}


############################################
#             MAIN PROGRAM
############################################
def_params <- list(  #overriden by some command line arguments
  exec_mode = "initialize",  #update
  arm_names = c("arm1", "arm2"),
  recomendaciones = 1, 
  some_other_param=NULL
)

params <- parse_cmdl(params=def_params, alt_params=list())

if(params$exec_mode == "initialize") {
  bandit_initialize(params=params)
} else if(params$exec_mode == "update") {
  bandit_update(params=params)
} else {
  print("Unrecognized execution mode:")
  print(params$exec_mode)
}
