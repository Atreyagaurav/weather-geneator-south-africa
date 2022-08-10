rm(list=ls())
library(depmixS4) # HMMs fit
library(rapportools)
library(markovchain)
library(rebmix)
library(moments)
library(MASS)
library(abind)
library(forecast)
library(biwavelet)
library(parallel)


#-- NOTE: Run auxiliary function "fit.NHMMs.depmix" beforehand -- #
#Aux#
#// Auxiliary function for NHMMs run
fit.long.NHMMs.depmix <- function(my.nstates,
                             my.num.sim,
                             my.synoptic.pcs,
                             my.dates.vec,
                             num.iteration.hmms,
                             modeltype,
                             init.pars) {
  
  #'  ---------- # s-NHMM for Annual WRs Simulations -------------------------------
  #'  Fits a non-homogeneous HMM with Seasonality  using 
  #'  depmixs4 package and simulates 
  #'  markov chains and response models with Seasonality 
  #'  
  #'  Arguments: 
  #' 
  #'  my.nstates : Number of hidden states of Weather Regimes
  #'  my.num.sim : Number of simulation of WRs (1 set of t (number of years), 2 set of t, etc)
  #'  my.synoptic.pcs : "Observed" series of PCs for geopotential heights 
  #'  my.dates.vec : vector of dates 
  #'  num.iteration.msar : number of simulation of Markov Chain 
  #'  modeltype : annual, seasonal, interannual


  # Prepare simulation output: [number of days (t) * number of sets of simulated WRs (entire stretch, t)] x number of ensemble members
  ts.length <- length(my.dates.vec)*my.num.sim
  matrix.hmms.seq.states <- array (NA,c(ts.length,num.iteration.hmms))
  
  lst.tpm.probs <- list()
  lst.real.eigen.vectors <- list()
  
  # count number of days assigned to the states #
  sum.state.days <- matrix(NA,nrow = my.nstates,ncol = num.iteration.hmms)
  
  # Seasonality 
  if(modeltype == "annual"){my.period = 365}
  if(modeltype == "seasonal"){my.period = 30}
  if(modeltype == "interannual"){my.period = 365}
  
  #define day of year for all dates
  tsteps <- as.numeric(format(my.dates.vec, "%j")) #seq_len(nrow(my.synoptic.pcs))
  
  #------- transition model -----------------------------
  # Seasonal Covariates for hidden states 
  if(modeltype == "annual"){
    omegaT <- (2*pi*tsteps)/my.period
    my.covar.trans <- data.frame(CosT = cos(omegaT), SinT = sin(omegaT)) 
    
    fo <- as.formula(" ~ -1 + CosT + SinT")
    #fo <- as.formula(" ~ 1 + CosT + SinT") 
    
  }else if(modeltype == "seasonal") 
  {
    omegaT <- (2*pi*tsteps)/my.period
    my.covar.trans <- data.frame(CosT = cos(omegaT), SinT = sin(omegaT)) 
    fo <- as.formula("~ 1 + CosT + SinT")
    
  } else if(modeltype == "interannual") 
  {
    omegaT <- (2*pi*tsteps)/my.period
    my.covar.trans <- data.frame(CosT = cos(omegaT), SinT = sin(omegaT),
                                 t.trend = tsteps) 
    fo <- as.formula("~ 1 + CosT + SinT + t.trend")
  }
  
  #A list of transition models, each created by a call to transInit. The length
  #of this list should be the number of states of the model.
  transition <- list() 
  if(is.null(init.pars)){   
    for(i in 1:my.nstates){
      transition[[i]] <- transInit(formula = fo, 
                                   data = my.covar.trans, 
                                   nstates = my.nstates)
    }
  } else {
    # initial transition parameters are provided 
    for(i in 1:my.nstates){
      set.seed(i*950)
      par.cosT <- runif(my.nstates,0.01,0.9)
      par.sinT <- 1 - par.cosT
      pstart. <- c(par.cosT,par.sinT)#init.pars$transition[[i]], # intercepts 
      
      pstart.[which(is.infinite(pstart.))] <- 0.9
      pstart.[which(is.nan(pstart.))] <- 0.01 
      
      transition[[i]] <- transInit(formula = fo, 
                                   data = my.covar.trans,  
                                   nstates = my.nstates,
                                   pstart = pstart.) 
      
      
    }
  }
  
  #----------------- prior --------------------------  
  if(is.null(init.pars)){
    my.prior <- transInit(~1,ns = my.nstates,data=data.frame(1),
                          ps = runif(my.nstates))
  }else {
    if(any(init.pars$prior == 0)){
      my.prior <- transInit(~1,ns = my.nstates,data=data.frame(1),
                            ps = runif(my.nstates))
    }else{
      my.prior <- transInit(~1,ns = my.nstates,data=data.frame(1),
                            ps = init.pars$prior)
    }
    
  }
  
  # #-----State-wise GLM response model ----------------
  my.response.models <- list() 
  if(is.null(init.pars)){
    for(i in 1:my.nstates){
      my.response.models[[i]] <- list(
        GLMresponse(formula = PC1~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC2~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC3~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC4~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC5~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC6~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC7~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC8~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC9~1,data = data.frame(my.synoptic.pcs ),family = gaussian()),
        GLMresponse(formula = PC10~1,data = data.frame(my.synoptic.pcs ),family = gaussian())
      )
    }
  } else {
    for(i in 1:my.nstates){
      my.response.models[[i]] <- list(
        GLMresponse(formula = PC1~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[1]]),
        GLMresponse(formula = PC2~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[2]]),
        GLMresponse(formula = PC3~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[3]]),
        GLMresponse(formula = PC4~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[4]]),
        GLMresponse(formula = PC5~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[5]]),
        GLMresponse(formula = PC6~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[6]]),
        GLMresponse(formula = PC7~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[7]]),
        GLMresponse(formula = PC8~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[8]]),
        GLMresponse(formula = PC9~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[9]]),
        GLMresponse(formula = PC10~1,data = data.frame(my.synoptic.pcs ),family = gaussian(),pstart = init.pars$response[[i]][[10]])
      )
    }
  }  
  
  #---------------- Initialize ---------------
  mod <- makeDepmix(response=my.response.models,
                    transition=transition,
                    prior=my.prior,
                    ntimes= nrow(my.synoptic.pcs),
                    homogeneous=FALSE)
  
  # Parameter Estimation : ----------------------------
  #if(is.null(init.pars)){
  tmp.mod.list <- list()
  for(j in 1:10){ 
    #
    set.seed(j*950)
    fmod.depmix <- fit(mod,emc = em.control(random = FALSE),
                       verbose = F) #conrows = conr.nh)  # )#
    # 
    tmp.mod.list[[j]] <-  fmod.depmix
    print(paste("Initial run",j))
    rm(fmod.depmix)
    #
  }
  
  # check the convergence message here before selecting the model 
  all.msgs <- sapply(tmp.mod.list, function(x){
    if(class(x) == "depmix.fitted"){
      return(x@message)
    }else {
      return(c())
    }
  })
  
  logLike.list <- as.numeric(unlist(lapply(tmp.mod.list,logLik)))  
  
  index.converged <- stringr::str_match(all.msgs, "Log likelihood converged") 
  
  if(all(is.na(index.converged))){
    print("EM did not converge with inital values. Recalculating with random starting values")
    
    tmp.mod.list <- list()
    for(j in 1:10){ 
      #
      set.seed(j*950)
      fmod.depmix <- fit(mod,emc = em.control(random = TRUE),
                         verbose = F) #conrows = conr.nh)  # )#
      
      tmp.mod.list[[j]] <-  fmod.depmix
      print(paste("Initial run",j))
      rm(fmod.depmix)
      #
    }
    
    all.msgs <- sapply(tmp.mod.list, function(x){
      if(class(x) == "depmix.fitted"){
        return(x@message)
      }else {
        return(c())
      }
    })
    
    logLike.list <- as.numeric(unlist(lapply(tmp.mod.list,logLik)))  
    
    index.converged <- stringr::str_match(all.msgs, "Log likelihood converged") 
    
  } else{
    print("Comparing Loglikelihood...")
  }
  
  
  logLike.list[which(is.na(index.converged))] <- 99999
  
  mod.num <- which.min(logLike.list) 
  # 
  fmod.depmix <- tmp.mod.list[[mod.num]]
  
  init.seed <- mod.num*1991 
  
  
  # -------------------------------
  prob.state <- forwardbackward(fmod.depmix)$gamma
  
  # 'Observed' or 'TRUE" state sequence 
  # state sequnce (using the Viterbi algorithm) #
  seq.state <- posterior(fmod.depmix,type='viterbi')$state
  # Probability of each state 
  delta.probs <- posterior(fmod.depmix,type='viterbi') %>% dplyr::select(-state)
  
  # # ---------------------------------------------------------- #
  # ------ Simulation ------------------------------------------ #
  # # ---------------------------------------------------------- #
  for (it.cnt in 1:num.iteration.hmms)
  {
    sim.fmod <- depmixS4::simulate(fmod.depmix,nsim = my.num.sim, seed = it.cnt)
    
    sim.seq.state <- sim.fmod@states # different from viterbi sequence 
    
    matrix.hmms.seq.states[,it.cnt] <- sim.seq.state
    
    prct.done <- round(it.cnt/num.iteration.hmms*100, digits = 2)
    if(prct.done%%5 == 0){print(paste("---finishing ",it.cnt,"out of",num.iteration.hmms," simulations -----:",prct.done,"%"))}
  }
  
  return(list(matrix.hmms.seq.states=matrix.hmms.seq.states,
              fitted.model = fmod.depmix,
              viterbi.seq = seq.state, 
              init.seed = init.seed))
  
}

#1#
#// Load processed 500mb-GPH hgt region data and dates
hgt.synoptic.region <- readRDS(file='hgt.500.Pacific.NorthAmer.synoptic.region_19480101_20211231.rds')

# for a specific WR number
# e.g, 10 PCs
#start_date="1948-01-01"; end_date="2021-12-31" # for TID/HFAM
start_date="1950-01-01"; end_date="2013-12-31" # for CA Watershed Study
first.date.weather <- as.Date(start_date); last.date.weather <- as.Date(end_date)
dates.weather <- seq(as.Date(start_date),as.Date(end_date), by="days")

start_date_synoptic="1948-01-01"; end_date_synoptic="2021-12-31"
dates.synoptic <- seq(as.Date(start_date_synoptic),as.Date(end_date_synoptic), by="days")

identical.dates.idx <- dates.synoptic%in%dates.weather
hgt.synoptic.region <- hgt.synoptic.region[identical.dates.idx,]


#2#
##/ Use PCA beforehand
hgt.synoptic.region.pca <- prcomp(hgt.synoptic.region,center=T,scale=T)
num_eofs <- 10
synoptic.pcs <- hgt.synoptic.region.pca$x[,1:num_eofs]


#3#
##/ HMMs runs followed by s-NHMMs
# for a specific WR number
num.states = 13 # WRs number
number.years.long <- 1000 # e.g., 1000 years; 2000 years, etc
my.num.sim = ceiling(number.years.long/length(unique(format(dates.synoptic[identical.dates.idx],'%Y')))) # number of chunks of historical periods; e.g., 1 is one set of simulation equal to the historical
number.years.long2 <- my.num.sim*length(unique(format(dates.synoptic[identical.dates.idx],'%Y')))
num.iteration.hmms = 5 # number of iteration to simulate WRs
lst.WRs.sNHMMs.states <- list()
# HMMs outputs
# e.g, 10 PCs
modHMMs <- depmix(list(PC1~1,PC2~1,PC3~1, PC4~1, PC5~1, PC6~1, PC7~1, PC8~1,
                       PC9~1, PC10~1),
                  nstates = num.states,
                  family=list(gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),gaussian(),
                              gaussian(),gaussian()),
                  ntimes =  nrow(synoptic.pcs),
                  data = data.frame(synoptic.pcs))
fit.modHMMs.depmix <- fit(modHMMs)
# synoptic.state.assignments <- posterior(fit.modHMMs.depmix)$state # state sequence (using the Viterbi algorithm) #
# weather.state.assignments <- synoptic.state.assignments[dates.synoptic%in%dates.weather]   #state assignments associated with the local weather variables (in case they span a different time period)

#s-NHMMs run
hhmod <- fit.modHMMs.depmix
hhpars <- c(unlist(getpars(hhmod)))
hhconMat <- hhmod@conMat
init.pars <- list()
init.pars[['transition']] <- lapply(hhmod@transition,
                                    function(x) x@parameters$coefficients)
init.pars[['prior']]  <- hhmod@prior@parameters$coefficients #  prob.kmeans.list[[p]]  #
init.pars[['response']] <- lapply(hhmod@response,
                                  function(x) lapply(x, function(y) unlist(y@parameters)))
init.pars[['conMat']] <- hhconMat

# s-NHMM on PCs with seasonality
tmp.list <- fit.long.NHMMs.depmix(my.nstates=num.states,
                                  my.num.sim=my.num.sim,
                                  my.synoptic.pcs=synoptic.pcs,
                                  my.dates.vec=dates.synoptic[identical.dates.idx],
                                  num.iteration.hmms,
                                  modeltype='annual',
                                  init.pars = init.pars)
lst.WRs.sNHMMs.states <- tmp.list
gc()
print(paste("-- finishing: ",num.states))


#saveRDS(lst.WRs.sNHMMs.states,
#        file = paste0("/home/fs01/nn289/WRs.large.ensemble_NHMMs/myOutput/lst.long.",number.years.long2,".yrs.WRs.sNHMMs.",num.states,".states.",num.iteration.hmms,".iter_v3.rds"))
saveRDS(lst.WRs.sNHMMs.states,
        file = paste0("out-lst.long.",number.years.long2,".yrs.WRs.sNHMMs.",num.states,".states.",num.iteration.hmms,".iter_long.CA.WshdStd.rds"))

