# Fit SSD to CCME data for several chemicals but allowing multiple measurements for each 
# endpoint

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

library(dclone)
library(fitdistrplus)
library(ggplot2)
library(grid)
library(gridGraphics)
library(plyr)
library(rjags)
source("../CCME-fits/mygofstat.R")

# Create the plotting directory
dir.create(file.path("Plots"))


#load test data
endpoints<-read.csv("Cd_BC_Multi.csv",header=TRUE,as.is=TRUE, strip.white=TRUE)
head(endpoints)

xtabs(~Source, data=endpoints, exclude=NULL, na.action=na.pass)
xtabs(~Species+Source, data=endpoints, exclude=NULL, na.action=na.pass)

# How variable are the repeated readings on the same species
species.sd <- plyr::ddply(endpoints, c("Source","Species"), plyr::summarize,
                    n   =length(Conc),
                    mean=mean(Conc),
                    mean.log=mean(log(Conc)),
                    sd  =sd  (Conc),
                    sd.log=sd(log(Conc))
                    )
# how many species are repeated
xtabs(~Source+n, data=species.sd, exclude=NULL, na.action=na.pass)

# plot the sd vs the mean 
ggplot(data=species.sd, aes(x=mean, y=sd))+
  ggtitle("SD of repeated readings vs mean")+
  geom_point()+
  facet_wrap(~Source)

ggplot(data=species.sd, aes(x=mean.log, y=sd.log))+
  ggtitle("SD of repeated readings vs mean - log scale")+
  geom_point()+
  facet_wrap(~Source)
# Good evidence that we need to work on the log-scale

endpoints$logConc <- log(endpoints$Conc)


#-----------------------------------------------------------------
#  Fit the distributions to the (geometric) mean of repeated measurements
#
endpoint.mean <- plyr::ddply(endpoints, c("Chemical","Species"), summarize, gmean = exp(mean(log(Conc))))
endpoint.mean$Conc <- endpoint.mean$gmean
endpoint.mean


# create list of possible distributions
fit.list <- list( lnorm  = list(data=NULL, distr="lnorm",  method="mle")) #,
#                  llog   = list(data=NULL, distr="llog",   method="mle"), # log logistic
#                  gomp   = list(data=NULL, distr="gompertz",method="mle"),
#                  lgumbel= list(data=NULL, distr="lgumbel",method="mle"), # log gumbel
#                  gamma  = list(data=NULL, distr="gamma",  method="mle"),                         
#                  pareto = list(data=NULL, distr="pareto", method="mle"),
#                  weibull= list(data=NULL, distr="weibull",method="mle")
#                  burr   = list(data=NULL, distr="burr",   method="mle")

# Fit each distribution from the list to each chemical
fit.all <- plyr::dlply(endpoint.mean, "Chemical", function(x, probs=c(.05,.10), nboot=1000){
      cat("Analyzing ", x[1,"Chemical"],"\n")
      # fit all of the distributions in the list to this data. If
      # a distribution fitting does not converge, then it returns NULL
      dist.fits <- plyr::tryapply(fit.list, function(distr, x){
           distr$data <- x$Conc
           if(distr$distr=="burr"){
              distr$start <-list(shape1=4, shape2=1, rate=1)
              distr$method<- "mme"
              distr$order <- 1:3
              distr$memp  <- function (x, order){sum(x^order)}
           }
           if(distr$distr=='gamma'){
             distr$start <- list(scale=var(x$Conc)/mean(x$Conc), 
                                    shape=mean(x$Conc)^2/var(x$Conc)^2)
           }
           if(distr$distr == 'gompertz'){
              # use the vgam to get the parameters of the fit
              fit <- vglm(Conc~1, gompertz, data=x)
              distr$start <- list(shape=exp(unname(coef(fit)[2])), scale=exp(unname(coef(fit)[1])) )
            }
           if(distr$distr=='lgumbel'){
             distr$start <- list(location=mean(log(x$Conc)), scale=pi*sd(log(x$Conc))/sqrt(6))
           }
           if(distr$distr=='pareto'){ #use the vgam package to estimate the starting values
             fit<- vglm(Conc~1, paretoff, data=x)
             distr$start  <- list(shape=exp(unname(coef(fit)))) 
             distr$fix.arg<- list(scale=fit@extra$scale)
           }
           if(distr$distr=="llog" ){ # log=logistic
              distr$start <- list(shape=mean(log(x$Conc)), scale=pi*sd(log(x$Conc))/sqrt(3))
           }
           cat("   Fitting ", distr$distr,"\n")
           res <- do.call("fitdist", distr)
           print(res)
           res
      },x=x)
      
      # Get the aic table by hand. 
      aic.table <- plyr::ldply(dist.fits, function(x){
        # extract the distribution name and AIC value
        aic <- x$aic
        distname <- x$distname
        aicc <- AICcmodavg::AICc(x) # conflict with VGAM package
        nparm <- length(x$estimate)
        data.frame(distname=distname, aic=aic, k=nparm, aicc=aicc)
      })
      aic.table <- aic.table[order(aic.table$aicc),]
      aic.table$delta.aicc <- aic.table$aicc - min(aic.table$aicc)
      aic.table$weight     <- round(exp(-aic.table$delta.aicc/2)/sum(exp(-aic.table$delta.aicc/2)),3)
      
      # make predictions for HC5 based on all of the distributions
      # Include a bootstrap approximation to the se of the estimates
      pred.table <- plyr::ldply(dist.fits, function(x, probs=.05, nboot=1000){
        distname <- x$distname
        pred     <- quantile(x, probs)
        fit.boot <- bootdist(x, niter=nboot)
        #browser()
        se       <- apply(quantile(fit.boot, probs=probs)$bootquant, 2, sd, na.rm=TRUE)
        lcl      <- quantile(fit.boot, probs=probs)$quantCI[1,]
        ucl      <- quantile(fit.boot, probs=probs)$quantCI[2,]
        #cat("pred.table", distname, pred, se, lcl, ucl, "\n")
        data.frame(distname=distname, quantile=probs, pred=unlist(pred$quantiles), se=se, lcl=unlist(lcl), ucl=unlist(ucl) )
      }, probs=probs, nboot=nboot)
      
      pred.table <- merge(pred.table, aic.table[,c("distname","weight")])
      
      # compute the model averaged quantiles, the model averaged lcl and ucl, and unconditional se
      q.modavg <- plyr::ddply(pred.table, "quantile", function(x){
        avg.pred <- sum(x$pred * x$weight)
        avg.lcl  <- sum(x$lcl  * x$weight)
        avg.ucl  <- sum(x$ucl  * x$weight)
        # get the unconditional se
        avg.u.se <- sqrt(sum(x$weight * (x$se^2 + (x$pred-avg.pred)^2)))
        data.frame(avg.pred=avg.pred, avg.u.se=avg.u.se, avg.lcl=avg.lcl, avg.ucl=avg.ucl)
      })
      q.modavg
      
      
      # make the cdf plot comparing all of the fits
      # capture the plot and save it
      # get the list of distributions that were fit
      dist.names <- plyr::laply(dist.fits, function(x) x$distname)
      cdfcomp(dist.fits, xlogscale=TRUE, legendtext=dist.names,
              main=paste("Empirical and theoretical CDF's - ",x[1,"Chemical"],sep=""))
      gridGraphics::grid.echo()
      cdf.comp.plot <- grid::grid.grab()
      
      #compute the goodmess of fit statistics - NOTE this does not work for burr distribution, not sure why...
      # must check to see if more than one distribution was fit - groan.
      # Also had to fix an error in compute p-value for chi2 test when df=0
      if(length(dist.fits) ==1) {gof.stat <- mygofstat(dist.fits[[1]], fitnames=dist.names)}
      if(length(dist.fits) >1 ) {gof.stat <- mygofstat(dist.fits,fitnames =dist.names)}
      
      # return the entire fitting stuff for this distibution
      list(Chemical=x[1,"Chemical"], dist.fits=dist.fits, aic.table=aic.table, pred.table=pred.table, 
           q.modavg=q.modavg,
           cdf.comp.plot =cdf.comp.plot,
           gof.stat =gof.stat)
})

display.results <- function(x){
   # display the results from the fitting functions
   cat("\n\n\n\n\n\n\n\n\n *******************************************************************\n")
   cat("\n\nResults when applied to ", x$Chemical, "\n")
   cat("\n\nAICc table \n")
   print(x$aic.table)
   
   cat("\n\nPredictions of endpoints \n")
   print(x$pred.table)
   
   cat("\n\nModel averaged endpoint \n")
   print(x$q.modavg)
   
   cat("\n\nCDF comparative plot\n")
   grid.newpage()
   grid.draw(x$cdf.comp.plot)
   
   cat("\n\nGoodness of fit statistics\n")
   print(x$gof.stat)
   
}


# display all of the results from all of the Chemicals
plyr::l_ply(fit.all, display.results)

# Save the final plot to the plotting directory
# save the comparative plot
plyr::l_ply(fit.all, function(x){
   file.name=file.path("Plots",paste(x$Chemical,"-comparative-plot.png",sep=""))
   png(file.name, h=6, w=6, units="in", res=300)
      grid.newpage()
      grid.draw(x$cdf.comp.plot)
   dev.off()
})


fit.all[[1]]


#-------------------------------------------------------------------------------------------
# Do the fit allowing for measurement error
#


# Get the data ready for jags
# We need to create a species number because jags doesn't deal with character indices
endpoint.mean <- plyr::ddply(endpoints, c("Species"), summarize, gmean = exp(mean(log(Conc))))
endpoint.mean
species.list <- endpoint.mean$Species
endpoints$SpeciesNum <- match(endpoints$Species, species.list)
xtabs(~Species+SpeciesNum, data=endpoints, exclude=NULL, na.action=na.pass)


# get the BUGS model and select for the distribution that is wanted
BUGS.model <- readLines(file.path("fit.distr.clone-lnorm.r"))
model.clone <- dclone::custommodel(BUGS.model)
write.jags.model(model.clone, "model.txt")

# Next create the data.list.

# data for number of data values and the number of unique species
NData <- nrow(endpoints)
NSpecies <- length(unique(endpoints$Species))

# data for density
logConc <- endpoints$logConc
SpeciesNum <- endpoints$SpeciesNum

# The datalist will be passed to JAGS with the names of the data values.
data.list.clone <- list(NData   = NData,
                        NSpecies= NSpecies,
                        SpeciesNum=SpeciesNum,
                        logConc    = dcdim(data.matrix(endpoints$logConc)),
                        K=1)
data.list.clone

# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.


init.list <- list(
  list(sigma_obs=mean(species.sd$sd.log, na.rm=TRUE), mu_species=endpoint.mean$gmean,  # observation process
       mu=mean(logConc), sigma=sd(logConc))   # initial values for log-normal distribution of species endpoints
)  # end of list of lists of initial values
init.list


# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("mu_species","sigma_obs","rank_mu_species", #"lik.data","lik.latent","lik",
                  "mu","sigma",  # for log-normal distribution
                  "hc","hc.log") # parameters to monitor

# Finally, the actual call to JAGS
set.seed(12342354)  # intialize seed for MCMC 

# Finally, the actual call to JAGS
results <- dc.fit(data    = data.list.clone,
                    params= monitor.list,
                    initss= init.list,
                    model = model.clone,
                    n.clones=c(1,2,4,8),
                    unchanged=c("NData","NSpecies","SpeciesNum"),
                    multiply="K")


summary(results)
coef(results)
dcsd(results)
sqrt(diag(vcov(results)))
confint(results)






# make a plot of the fitted curve, the observed data, the estimated mean of obs data, etc
# Extract the mean and variance of the fit
mu    <- coef(results)["mu"]
sigma <- coef(results)["sigma"]

mu.se <- dcsd(results)["mu"]
sigma.se<-dcsd(results)["sigma"]


# Compute the fitted curve from the cloned data
fitted.clone     <- data.frame(logConc=seq(min(log(endpoints$Conc)),max(log(endpoints$Conc)),length.out=100))
fitted.clone$cdf <- plnorm(exp(fitted.clone$logConc), meanlog=mu, sdlog=sigma)

# Compute the fitted curve when the geometric mean was used
names(coef(fit.all[[1]]$dist.fits$lnorm))
fitted.avg     <- data.frame(logConc=fitted.clone$logConc)
fitted.avg$cdf <- plnorm(exp(fitted.avg$logConc), 
                           meanlog=coef(fit.all[[1]]$dist.fits$lnorm)["meanlog"], 
                           sdlog  =coef(fit.all[[1]]$dist.fits$lnorm)["sdlog"])


# extract the log(HC5) and actual HC5 from the cloning MLE
hc.clone        <- coef(results)["hc"]
hc.ci.clone     <- quantile(results, prob=c(.025, .975))[,"hc"]
hc.log.clone    <- coef(results)["hc.log"]
hc.log.ci.clone <- quantile(results, prob=c(.025, .975))[,"hc.log"]


# extract the log(HC5) from the fitted distribution
names(fit.all[[1]]$pred.table)
fit.all[[1]]$pred.table
hc     <- fit.all[[1]]$pred.table[ fit.all[[1]]$pred.table$quantile==0.05, "pred"]
hc.ci  <- fit.all[[1]]$pred.table[ fit.all[[1]]$pred.table$quantile==0.05, c("lcl","ucl")]
hc.log <- log(hc)
hc.log.ci <- log(hc.ci)

# Extract the fitted mean for each species
select<- grepl("^mu_",names(coef(results)))
mu_species <- data.frame(mean=coef(results)[select],
                         t(quantile(results, prob=c(.025, .975))[,select]))
mu_species$logmean <- log(mu_species$mean)
mu_species$logmean.lcl <- log(mu_species$X2.5.)
mu_species$logmean.ucl <- log(mu_species$X97.5.)
mu_species$SpeciesNum <- as.numeric(substr(row.names(mu_species),
                                1+regexpr("[", row.names(mu_species), fixed=TRUE), 
                               -1+regexpr(",", row.names(mu_species), fixed=TRUE)))
select<- grepl("^rank_",names(coef(results)))
mu_species$avgrank <- coef(results)[select]
mu_species$avgecdf <- (mu_species$avgrank-0.5)/length(unique(mu_species$SpeciesNum))
# average over the clones for this latent variable
mu_species <- plyr::ddply(mu_species, "SpeciesNum", summarize,
                          mean=mean(mean),
                          X2.5. = mean(X2.5.),
                          X97.5.= mean(X97.5.),
                          logmean=mean(logmean),
                          logmean.lcl=mean(logmean.lcl),
                          logmean.ucl=mean(logmean.ucl),
                          avgrank = mean(avgrank),
                          avgecdf = mean(avgecdf))
mu_species$ecdf <- ecdf(mu_species$logmean)(mu_species$logmean)


plotdata <- merge(endpoints, mu_species, by="SpeciesNum")
plotdata <- merge(plotdata, endpoint.mean, by="Species")
plotdata


plot1 <- ggplot(data=plotdata, aes(y=ecdf))+
   ggtitle("Fit of log-normal to BC Cd data (red=AVG, blue=MLE)")+
   ylab("Cumulative probability")+xlab("log(Concentration)")+
   geom_point(data=mu_species, aes(x=logmean, y=avgecdf), color="blue")+
 #  geom_errorbarh(aes(xmin=logmean.lcl, x=logmean, xmax=logmean.ucl), 
 #           color="blue",height=.05, size=.1)+

   geom_point(aes(x=logConc),shape=1, color="red",size=2)+
   geom_line(data=fitted.clone, aes(x=logConc, y=cdf), color="blue")+
   geom_line(data=fitted.avg,   aes(x=logConc, y=cdf), color="red")+
   geom_point( aes( x=log(gmean)), shape="X", size=3, color="red")+
   geom_hline(yintercept=.05)+
   geom_segment(data=plotdata, aes(x=log(gmean), y=ecdf, xend=logmean, yend=avgecdf))+
   annotate("text", x=Inf, y=.5, hjust=1, vjust=.5,
           label=paste("log HC5 AVG=",format(round(hc.log,2),nsmall=2),
                          "\n  HC5 =",format(round(hc,2),nsmall=2),
                     "  (",format(round(exp(hc.log.ci[1]),2),nsmall=2),",",
                           format(round(exp(hc.log.ci[2]),2),nsmall=2),")",sep=""), color="red")+
   annotate("text", x=Inf, y=.3, hjust=1, vjust=.5,
           label=paste("log HC5 MLE=",format(round(hc.log.clone,2),nsmall=2),
                          "\n  HC5 =",format(round(hc.clone,2),nsmall=2),
                     "  (",format(round(exp(hc.log.ci.clone[1]),2),nsmall=2),",",
                           format(round(exp(hc.log.ci.clone[2]),2),nsmall=2),")",sep=""), color="blue")+
   geom_errorbarh(data=NULL, aes(xmin=hc.log.ci[1], x=hc.log, xmax=hc.log.ci[2], y=.06), 
            color="red",height=.05, size=1)+
   geom_errorbarh(data=NULL, aes(xmin=hc.log.ci.clone[1], x=hc.log.clone, xmax=hc.log.ci.clone[2], y=.04), 
            color="blue",height=.05, size=1)
plot1
ggsave(plot1, file=file.path("Plots","BC_Cd_results-lnorm.png"), h=6, w=6, unit="in",dpi=300)






