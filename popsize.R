####################################################
# 	POPULATION MODEL FOR HeV DYNAMIICS
# 	------------------------------------------------
#   CALCULATING AVERAGE ANNUAL POPULATION SIZE (v) 
#   FOR A GIVEN SET OF DEMOGRAPHIC PARAMETERS (s,m,phi)
# 	------------------------------------------------
####################################################

# In this series I vary:
# - phase parameter (phi, range of values stored in phi_r)

# Constants:
# - s = 130
# - m = 0.14yr-1
# - v = 50000

source('~/birthpulse.R')

#----------Deterministic Demographic Model----------
dNdt <-function(t,y,p){
  phi <- p$phi
  m <- p$m
  s <- p$s
  N <- y
  births <- birthpulse(t,p)
  dN <- N*(births-m)
  return(list(dN))
}

#-----Ensuring that the average annual B(t) is equal to the constant death rate (m) -----
set_m <- 0.14/365                          # Avg lifespan is 7.14 yrs, mortality rate = 0.14yr^-1
phi_r <- seq(-pi/2,pi/2,length.out = 13)
set_s <- 130
timesteps <- seq(0,1*365,1e-1)             # Deterministic solution for 1 year
set_v <- 50000
B.det <- lapply(1:length(phi_r),function(phi){
  birthpulse(timesteps,
                    list(s=set_s,m=set_m,phi=phi))})
m.det <- sapply(1:length(phi_r),function(phi){
  mean(B.det[[phi]])
  })

#-----Calculating N0 for a given v (avg population)-----
N.det <- lapply(1:length(phi_r), function(phi_ind){
                    params=list(s=set_s,m=m.det[phi_ind],phi=phi_r[phi_ind])
                    ode(set_v, timesteps, dNdt, params)
                    })
N.avg <-sapply(1:length(phi_r),function(phi_ind){
  mean(N.det[[phi_ind]][,2])
})
v_norm <- N.avg/set_v
v2_r <- set_v/v_norm
v3_r <- round(v2_r,0)

#-----Checking if N0 (stored in v3_r) gives set_v-----
N.det.2<-lapply(1:length(phi_r),function(phi_ind){
  params=list(s=set_s,m=m.det[phi_ind],phi=phi_r[phi_ind])
  ode(v3_r[phi_ind],timesteps,dNdt,params)
})
N.avg.2 <- sapply(1:length(phi_r),function(phi_ind){
  mean(N.det.2[[phi_ind]][,2])
})
B.det.2 <- lapply(1:length(phi_r),function(phi_ind){ #Using m.det calculated from B.det. and same as N.det.2
  birthpulse(timesteps,list(s=set_s,m=m.det[phi_ind],phi=phi_r[phi_ind]))
})
