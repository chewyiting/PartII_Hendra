####################################################
# 	MODELLING HENDRA VIRUS DYNAMICS (SIR)
# 	------------------------------------------------
#   EFFECT OF RECOVERY RATE ON PATHOGEN PERSISTENCE
# 	------------------------------------------------
####################################################

# In this series I vary:
# - recovery rate (gamma, range of values stored in g_r)
# - transmission rate (beta, range of values stored in b_r)

# Constants:
# - R0 = 4 (set_R0)
# - m = 0.14 

source("~/birthpulse.R")

#----------Stochastic SIR model----------
library(adaptivetau) #Implement Cao et al's adaptive tau-leaping algorithm
# Transition matrix
transitions = list(c(S=+1), #birth
                   c(S=-1), #death.S
                   c(I=-1), #death.I
                   c(R=-1), #death.R
                   c(S=-1,I=+1,X=+1), #infection, transmission events recorded in demivariable X
                   c(I=-1,R=+1)) #recovery

# Rates of each transition (rateFunc input)
probabilities <- function(x,p,t){
  if(x["I"]!=0){ 
  return(c(
    (sum(x["S"],x["I"],x["R"]))*(birthpulse(t,p)),            #birth
    p$m*x["S"],                                               #death.S
    p$m*x["I"],                                               #death.I
    p$m*x["R"],                                               #death.R
    p$b*x["S"]*x["I"]*(1/(1e-10+sum(x["S"],x["I"],x["R"]))),  #infection
    p$g*x["I"]                                                #recovery
  ))} else { #If I = 0, halt simulation by making all rates = 0
    return(c(0,0,0,0,0,0))
  } 
}

#----------Stochastic simulations----------
library(compiler)  

# Parameters 
ntrials <- 1000       # Number of simulations for given set of conditions
tmax <- 10*365        # 10 years
set.seed(3)           # For reproducible simulated values

set_phi <- 0
set_m <- 0.14/365
set_s <- 130

g_r <- (seq(1,12,1))/365

# Calculate initial.v from popsize.R

# SIR_v is a function returning a matrix containing all simulations for explored range of recovery rates
# Matrix dimensions defined by length of g_r and ntrials

SIR_v <- function(recovery,initial.v){
  n.cond1 <- length(recovery)
  l <- vector("list",n.cond1*ntrials)
  dim(l) <- c(n.cond1,ntrials)
  for (g_ind in 1:n.cond1){
    init.values=c(S=initial.v,I=1,R=0,X=0)
    print(paste(c("gamma=",recovery[g_ind]))
    params=list(s=set_s,m=set_m,phi=set_phi,
                b=b_r[g_ind],g=g_r[g_ind])
    for(i in 1:ntrials){
      l[g_ind,i] <- list(ssa.adaptivetau(
        init.values,transitions,probabilities,params,
        tf=tmax,tl.params=list(maxtau=15/24))) # Max tau limited to 15 hours, to capture small changes in birth rate.
    }
  }
  l
}
enableJIT(1) 
system.time(SIR_ccs.1 <- SIR_v(g_r,initial.v)) 
saveRDS(SIR_ccs.1,"SIR_ccs.1.rds")
