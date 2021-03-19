####################################################
# 	MODELLING HENDRA VIRUS DYNAMICS (SIR)
# 	WITHIN A HOST POPULATION 
#   WITH AN ANNUAL, SEASONAL BIRTH PULSE 
# 	AND A CONSTANT DEATH RATE
####################################################

# Yi-Ting Chew, 2021

# An annual, seasonal birth rate function (developed by Peel et al., 2014)
# Birth rate function uses parameter estimates from Peel et al., 2014, who fitted the function to a dataset of P. poliocephalus in Gordon, NSW obtained by Peggy Eby
# Constant death rate is kept equal to annual average birth rate for stable average inter-annual population size (v), such that dN/dt = 0.

# Demographic structure :
# s = 130 (synchrony parameter, corresponds to 95% of births occurring in 28days)
# m = 0.14 (per capita death rate, corresponds to an average lifespan of 7.14 years)
# phi (phase parameter, ranges from -pi/2 to pi/2 to change timing of pathogen introduction relative to birth pulse)

#----------Birth rate function----------
library(deSolve)
birthpulse <- function(t,p){
  phi <- p$phi
  m <- p$m
  s <- p$s
  if(s==0){
    birth <- m
  } else {
    ki <- m/besselI(s/2,0,TRUE) 
    cos_term <- exp(-s*(cos(-phi+pi*(t/365)))^2)
    birth <- ki*cos_term
    }
  return(birth)
}
