# simple example of data cloning taken from their web pages (with slight modification)

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
n <- 25
p <- .3
y <- rbinom(n=1, size=n, p=p)


library(dclone)
library(rjags)


model <- dclone::custommodel("
model {
   y ~ dbin(p, n)  # likelihood
 
   p ~ dunif(.001, .999)  # prior

   # compute the likelihood
   loglik <- log( pbin(y, p, n))
}
")


data.list <- list(y=y, n=n)
data.list


monitor.list <- c("p","loglik")

fit <- jags.fit(data   =data.list, 
                params =monitor.list,
                model  = model)

summary(fit)
plot(fit)


#---------------------------------------------------
# Now for the data cloning

model.clone <- dclone::custommodel("
model {
  for(k in 1:K){  # cloning loop
    y[1,k]   ~ dbin(p, n)  # likelihood
    loglikc[1,k] <- log( dbin(y[1,k], p, n))
  }                          
  
  loglik <- mean(loglikc[1,1:K])  # the mean loglikelihood?
  p ~ dunif(.001, .999)  # prior
 }
")

data.list.clone <- list( y=dcdim(data.matrix(y)), n=n, K=1)
data.list.clone

monitor.list <- c("p", "loglik")

fit.clone <- dc.fit(data  =data.list.clone,
                    params=monitor.list,
                    model =model.clone,
                    n.clones=c(1,2,4,8,16),
                    unchanged="n",
                    multiply="K")

summary(fit.clone)
phat <- summary(fit.clone)$statistics["p","Mean"]
phat
# actual log likelihood
dbinom(y, n, phat, log=TRUE )



coef(fit.clone)
phat <- y/n
cat("MLE ", phat, "\n" )
dcsd(fit.clone)
cat('Theoretical SE ', sqrt(p*(1-p)/n), "\n")
cat("Empirical   SE ", sqrt(phat*(1-phat)/n), "\n")

vcov(fit.clone)
confint(fit.clone)
