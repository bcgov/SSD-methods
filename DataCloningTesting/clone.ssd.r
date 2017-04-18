# Fit SSD to CCME data using data cloning and compare results to straight likelihood method
# I will ignore for now the repeated measurements on the same species
# This just illustrates the congruence of results between the two methods

# This is a Bayesian analysis using data cloning to get the MLEs (see Lele's paper)

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


set.seed(23423423)

library(dclone)
library(fitdistrplus)
library(ggplot2)
library(grid)
library(plyr)
library(rjags)

#load test data
endpoints<-read.csv("Cd_BC_Multi.csv",header=TRUE,as.is=TRUE, strip.white=TRUE)
endpoints$Source <- "BC Cd"
head(endpoints)

sd(log(endpoints$Conc))

###############################################################################
# Fit a Bayesian model for the SSD using a
# log-normal model

model.clone <- dclone::custommodel("
model {
  # Input data for Conc   - concentration
  #                Ndata - number of data values
  #
  # Parameters for lognormal distribution
  #   mu, sigma - intercept and slope (on log scale)

  
  # likelihood function for observed data for each clone
  for(k in 1:K){
     for(i in 1:Ndata){
       Conc[i,k]  ~ dlnorm(mu, tau)
     }
  }
  
  # compute the log-likelihood (only need to do for the original set of data in column 1)
  for(i in 1:Ndata){
      loglikc[i] <- log(dlnorm(Conc[i,1], mu, tau))
  }
  
  loglik <- sum(loglikc[])  #  loglikelihood?
  
  # priors and derived variables
  tau <- 1/(sigma*sigma)
  sigma  ~ dunif(.01, 5)
  
  mu ~ dnorm(0, .001)
  
  # estimate the HC[5]
  hc.log <- qnorm(.05, mu, tau)  # the HC[5] on the log scale
  hc     <- exp(hc.log)
  
}  # end of model.jags
")



# Next create the data.txt file.

# data for entrainment
Ndata <- nrow(endpoints)

# The datalist will be passed to JAGS with the names of the data values.
data.list.clone <- list(Ndata   = Ndata,
                       Conc    = dcdim(data.matrix(endpoints$Conc)),
                       K=1)
data.list.clone

# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

init.list <- list(
  list(mu=mean(logConc), tau=1/var(logConc), tau_obs=1/mean(species.sd$sd.log, na.rm=TRUE)^2)
)  # end of list of lists of initial values



# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("mu","sigma","loglik",
                  "hc","hc.log") # parameters to monitor

# Finally, the actual call to JAGS
fit.clone <- dc.fit(data  =data.list.clone,
                    params=monitor.list,
                    model =model.clone,
                    n.clones=c(1,2,4,8,16),
                    unchanged=c("Ndata"),
                    multiply="K")


summary(fit.clone)
coef(fit.clone)
dcsd(fit.clone)
vcov(fit.clone)
confint(fit.clone)


# use fitdistrplus to estimate the parameters
library(fitdistrplus)

fit <- fitdist( endpoints$Conc, "lnorm")
fit
logLik(fit)
quantile(fit, probs=.05)
fit.boot <- bootdist(fit, niter=200)
quantile(fit.boot, probs=.05)
