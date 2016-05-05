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

## model with only direct effect of temperature and no indirect effect via body size
# should replicate Gilbert et al. 2014
rmcr_temp <- function(t,y,p) {
  # t,y,p is the necessary format for solving using the ode() function in R
  R <- y[1]
  C <- y[2]
  with(as.list(p), {
    dR.dt <- rT * R * (1 - R / KT) - aT * C * R 
    dC.dt <- e * aT * C * R - mT * C # assuming e is independent of temperature, Peters 1983
    return(list(c(dR.dt,dC.dt)))
  })
}

## model incorporating temperature-size relationships
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

## state variable values ----
# initial values at beginning of "experiments"
R <- 0.6
C <- 0.1
i.state <- c(R=0.6,C=0.1) 

## temperature-independent parameters ----
## base parameters
r0 <- 2 # intrinsic growth rate; From Gilbert et al. 2014, Fig. 3 legend
e <- e0 <- 0.15 # conversion efficiency; From Gilbert et al. 2014, Fig. 3 legend
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
TR <- 15 + 273.15 # body temperature of resource in Kelvin, starting temperature as in Gilbert et al. 2014 Fig. 3
TC <- TR # body temperature of consumer in Kelvin
TS <- TC <- TR # temperature of supply rate? in Kelvin

## exponents for mass-scaling. Assuming the are all the same for right now
# parameters from DeLong et al. 2015, Am. Nat. Table 1 (protist consumer and algae resource) are in between notes
s.r <- 0.75 # -0.2 # mass-scaling for intrinsic growth rate
s.K <- 0.75 #-0.81 # mass-scaling for carrying capacity
s.a <- 0.75 # 1 # mass-scaling for attack rate
s.e <- 0.75 #-0.50 # mass-scaling for conversion efficiency
s.m <- 0.75 # -0.29 # mass-scaling for mortality rate

## parameters for calculating changes in size due to increasing temperature.
bR <- -5*10^-4 # arbitrarily chose a negative number 
bC <- bR # assuming slope for consumer is the same as for the resource
M0R <- 1 # base adult body size of resource, arbitrarily set to one.
M0C <- 3 # assuming base adult body size of consumer is 3x greater than resource

## effects of temperature on average adult body size ----
# need help writing this function
MR <- M0R + bR*TR # linear function, probably needs to be fine tuned. M0R is baseline adult body size; bR is the slope of the relationship between adult body size and temperature
MC <- M0C + bC*TC # linear function, probably needs to be fine tuned, see above notes.

## mass-dependent parameters ----
# from DeLong et al. 2015, Am. Nat.
rM <- r0*MR^s.r 
KM <- K0*MR^s.K
aM <- a0*MC^s.a
eM <- e0*MC^s.e

## temperature-dependent parameters ----
# from Gilbert et al. 2014, Ecol. Letts.
rT <- r0*exp(-EB/(k*TR)) # Gilbert didn't have r0...why not?
KT <- K0*exp(EB/(k*TR) - ES/(k*TS)) # delta K is thought to be less than 1 (Kratina et al. 2012, aquatic microcosms; Roemmich and McGowan 1995, marine)
mT <- m0*exp(-Em/(k*TC)) 
aT <- a0*sqrt(( v0C^2 * exp((-2*EvC)/(k*TC)) + v0R^2 * exp((-2*EvR)/(k*TR )))) # dependent on type of foragin interaction (Dell et al. 2013). For example, active-capture depends on body velocity of both consumer and resource, sit-and-wait depends primarily on body velocity of resource, grazing depends primarily on body velocity of consumer

## combining temperature- and mass-dependent parameters ----
rMT <- r0*(MR^s.r)*exp(-EB/(k*TR))
KMT <- K0*(MR^s.K)*exp(EB/(k*TR) - ES/(k*TS))
aMT <- a0*(MC^s.a)*exp(-Em/(k*TC))
eM # no appartent direct effect of temperature...but now, it appears there could be an indirect effect of temperature mediated by changes in body size

## Let's try running the simulation... ----
Time <- 1000 # time steps for simulation

# Just looking at direct effects of temperature, not working...
dynamic_matplot(param.vector = c(rT = rT, KT = KT, aT = aT, e = e, mT = mT), init.state = i.state, sim.length = Time, model = rmcr_temp, ylim = c(0,KT))

# Looking at direct and indirect effects of temperature, not working...
dynamic_matplot(param.vector = c(rMT = rMT, KMT = KMT, aMT = aMT, eM = eM, mMT = mMT), init.state = i.state, sim.length = Time, model = rmcr_tempsize, ylim = c(0,KMT))

# the simulation below shows that the code should work. It uses the parameters from Gilbert et al. 2014, Fig. 3
p.rmcr_temp <- c(rT = 2, KT = 100, aT = 0.1, e = 0.15, mT = 0.6) # define parameters
Time <- 1000 
rm_test_data <- ode(i.state,1:Time, rmcr_temp, p.rmcr_temp) # run simulation and get data
rm_test_data[Time, "C"]/rm_test_data[Time, "R"] # C:R biomass ratio = 0.3, but in Gilbert et al. 2014, Fig. 3, the C:R biomass ration should equal 0.9... 

dynamic_matplot(param.vector = p.rmcr_temp, init.state = i.state, sim.length = Time, model = rmcr_temp, ylim = c(0,100))


