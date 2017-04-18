# Copyright 2017 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.


###############################################################################
# Fit a Bayesian model for the replicated observations at the same species
# using data cloning.
#
# This function can be used for different distributions 
# by selectively extracting the code with the tag number at the right

model {
  # Input data for logConc - log-concentration
  #                SpeciesNum - species number to link replicate observations
  #                NData - number of data values
  #                NSpecies - number of unique species
  #
  # Parameters for distributions vary depending on the distribution
  #    lnorm   -   mu, sigma - intercept and slope (on log scale) - lognormal
  
  # Parameters for observal error
  #       sigma_obs - sd of observational error
  #
  
  
  # observed data comes from underlying mean with a lognormal distribution
  # this is true regardless of the distribution of the mu_species whose
  # distributions can change
  # We want separate latent variables for each clone because this sample size
  # is NOT increasing as we increase the number of chemicals tested.
  for(k in 1:K){
     for(i in 1:NData){
        logConc[i,k]  ~ dnorm(log_mu_species[SpeciesNum[i],k], tau_obs)
     } 
  }
  # Contributions to the likelihood from the observed data#
#  for(i in 1:NData){
#     likc.data[i] <- dnorm(logConc[i,1], mu_species[SpeciesNum[i]], tau_obs)
#  }
#  lik.data <- prod(likc.data[])  # total contribuiton to the likelihood from data
  # priors and derived variables for the observation process
  tau_obs <- 1/(sigma_obs*sigma_obs)
  sigma_obs ~dunif(.01, 5)
  
  #--------------------------------------------------------------------
  # LOGNORMAL distribution for the "mean" of each species 
  # fit a distribution to the mean of each observation
  for(k in 1:K){
     for(i in 1:NSpecies){
        mu_species[i,k] ~ dlnorm(mu, tau)
        log_mu_species[i,k] <- log(mu_species[i,k])
#        likc.latent[i] <- dlnorm(mu_species[i], mu, tau) # contribution from the likelihoo
     }
  }
#  lik.latent <- prod(likc.latent[])  # total contribution to the likelihood from latent species value
  tau <- 1/(sigma*sigma)
  sigma ~ dunif(.01,5)
  
  mu ~ dnorm(0, .001)
  
  # estimate the average ranking of the species for each clone
  for(k in 1:K){
     rank_mu_species[1:NSpecies,k] <- rank(mu_species[1:NSpecies,k])
  }
  
  # estimate the HC[5]
  hc     <- qlnorm(.05, mu, tau)  # the HC[5] 
  hc.log <- log(hc)
  #------------------------------------------------------------------
  
  # total complete-data likelihood
#  lik <- lik.data * lik.latent
}  # end of model.jags

