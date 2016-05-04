######################################################
## Description: This script simulates the direct and indirect (via body size) effects of temperature on consumer-resource dynamics
## Code author(s): Matt Barbour
## Email: barbour@zoology.ubc.ca
######################################################

## load required libraries ----
library(deSolve)

## functions for simulating consumer-resource dynamics ----
# general assumptions:
# (1) logistic growth for the resource(s)
# (2) type 1 functional response of consumer 
# (3) density-independent mortality for consumer 
# function is in the necessary format "function(t,y,p)" for solving used the ode() function in the R package "deSolve"
rmcr_tempsize <- function(t,y,p) {
  # t,y,p is the necessary format for solving using the ode() function in R
  R <- y[1]
  C <- y[2]
  with(as.list(p), {
    dR.dt <- rMT * R * (1 - R / KMT) - aMT * C * R 
    dC.dt <- eM * aMT * C * R - mMT * C # assuming e is independent of temperature, Peters 1983, but depends on body size which is indirectly affected by temperature...
    return(list(c(dR.dt,dC.dt)))
  })
}

## function for plotting consumer-resource dynamics over time ----
dynamic_matplot <- function(param.vector, # vector of parameters for the dynamical model
                            init.state, # initial state variables for the model
                            sim.length, # length of simulation
                            model, # dynamical model
                            ylim, # limits of state variables for plotting
                            ...){
  
  # Run the experiment
  df <- ode(init.state, 1:sim.length, model, param.vector)
  
  # plot the results. 
  matplot(df[ ,"time"], df[ ,names(init.state)],
          type = "l", ylim=ylim, ...)
}

## temperature-dependent parameters ----
# from Gilbert et al. 2014, Ecol. Letts.
rT <- exp(-EB/k*TR) 
KT <- K0*exp(EB/k*TR - ES/k*TS) # delta K is thought to be less than 1 (Kratina et al. 2012, aquatic microcosms; Roemmich and McGowan 1995, marine)
mT <- m0*exp(-Em/k*TC) 
aT <- a0*sqrt(( v0C^2 * exp(-2*EvC/k*TC) + v0R^2 * exp(-2*EvR/k*TR ))) # dependent on type of foragin interaction (Dell et al. 2013). For example, active-capture depends on body velocity of both consumer and resource, sit-and-wait depends primarily on body velocity of resource, grazing depends primarily on body velocity of consumer

## mass-dependent parameters ----
# from DeLong et al. 2015, Am. Nat.
rM <- r0*MR^s.r 
KM <- K0*MR^s.K
aM <- a0*MC^s.a
eM <- e0*MC^s.e

## effects of temperature on average adult body size ----
# don't know if it is accurate that I want to be calculating the change in temperature because I don't have the change in temperature anywhere else in the model...
# need help writing this function
MR <- M0R + bR*M0R*(TR - T0R) # function needs to be fine tuned. (TR - T0R) describes change in temperature; M0R is baseline adult body size; bR is the slope of the relationship between adult body size and temperature
MC <- M0C + bC*M0C*(TC - T0C) # function needs to be fine tuned, see above notes.

## combining temperature- and mass-dependent parameters ----
rMT <- r0*(MR^s.r)*exp(-EB/k*TR)
KMT <- K0*(MR^s.K)*exp(EB/k*TR - ES/k*TS)
aMT <- a0*(MC^s.a)*exp(-Em/k*TC)
eM # no appartent direct effect of temperature...but now, it appears there could be an indirect effect of temperature mediated by changes in body size

## temperature-independent parameters ----
## base parameters
r0 <- 2 # intrinsic growth rate; From Gilbert et al. 2014, Fig. 3 legend
e0 <- 0.15 # conversion efficiency; From Gilbert et al. 2014, Fig. 3 legend
K0 <- 100 # carrying capacity; From Gilbert et al. 2014, Fig. 3 legend
a0 <- 0.1 # attack rate; From Gilbert et al. 2014, Fig. 3 legend 
m0 <- 0.6 # mortality rate; From Gilbert et al. 2014, Fig. 3 legend

## activation energies
EB <- 0.32 # activation energy of the metabolic rate; EB = 0.53 - 0.85 Savage et al. 2004, diverse ectotherms; exact value from Gilbert et al. 2014, Fig. 3 legend
Em <- 0.65 # activation energy of mortality rate; Em = 0.65 Dell et al. 2011 across diverse taxa; Em = 0.45 Savage et al. 2004 for fish
ES <- 0.9 # activation energy of nutrient supply. From Gilbert et al. 2014, Fig. 3 legend
EvC <- 0.46 # activation energy of consumer velocity. 0.46 corresponds to the mean activatiaion energy for body velocity across diverse taxa and trophic groups Dell et al. 2014, J. Animal Ecology.
EvR <- 0.46 # activation energy of resource velocity. 0.46 corresponds to the mean activatiaion energy for body velocity across diverse taxa and trophic groups Dell et al. 2014, J. Animal Ecology.

## other parameters. Note that if body velocities for consumers and resources are about equal, then this represents an active-capture foraging strategy.
v0C <- 1 # consumer body velocity. arbitratily set to 1
v0R <- 1 # resource body velocity. arbitrarily set to 1
k <- 8.62*10^-5 # Boltzmann's constant in eV/K, where K is temperature in Kelvin

## temperatures. assuming they are all the same for right now
TR <- 15 # body temperature of resource in C, but don't I need it in Kelvin to be consistent with Boltmann's constant
TC <- TR # body temperature of consumer in C
TS <- TC <- TR # temperature of supply rate? in C

## exponents for mass-scaling. From DeLong et al. 2015, Am. Nat. Table 1 which was for a protist consumer and algae resource
s.r <- -0.2 # mass-scaling for intrinsic growth rate
s.K <- -0.81 # mass-scaling for carrying capacity
s.a <- 1 # mass-scaling for attack rate
s.e <- -0.50 # mass-scaling for conversion efficiency
s.m <- -0.29 # mass-scaling for mortality rate

## parameters for calculating changes in size due to increasing temperature.
bR # 
bC #
M0R #
M0C #