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
rM <- r0*Mr^s.r 
KM <- K0*Mr^s.K
aM <- a0*Mc^s.a
eM <- e0*Mc^s.e

## effects of temperature on average adult body size ----
Mr # need this function
Mc # need this function 

## combining temperature- and mass-dependent parameters ----
rMT <- r0*(Mr^s.r)*exp(-EB/k*TR)
KMT <- K0*(Mr^s.K)*exp(EB/k*TR - ES/k*TS)
aMT <- a0*(Mc^s.a)*exp(-Em/k*TC)
eM # no appartent direct effect of temperature...but now, it appears there could be an indirect effect of temperature mediated by changes in body size

## temperature-independent parameters ----
## activation energies
EB <- 0.65 # activation energy of the metabolic rate; EB = 0.53 - 0.85 Savage et al. 2004, diverse ectotherms
Em <- 0.65 # activation energy of mortality rate; Em = 0.65 Dell et al. 2011 across diverse taxa; Em = 0.45 Savage et al. 2004 for fish
ES <- 0.65 # activation energy of nutrient supply, just set to 0.65 because I don't know of possible values
EvC # activation energy of consumer velocity
EvR # activation energy of resource velocity

## other parameters
v0C # consumer body velocity
v0R # resource body velocity
k <- 8.62*10^-5 # Boltzmann's constant in eV/K, where K is temperature in Kelvin
TR # body temperature of resource in Kelvin
TC # body temperature of consumer in Kelvin
TS # temperature of supply rate? in Kelvin

## exponents for mass-scaling
s.r # mass-scaling for intrinsic growth rate
s.K # mass-scaling for carrying capacity
s.a # mass-scaling for attack rate
s.e # mass-scaling for conversion efficiency