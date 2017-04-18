# Fit SSD to CCME data for several metals

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


# We compare fits for the distributions with AIC
# Refer to https://cran.r-project.org/web/packages/fitdistrplus/vignettes/FAQ.html
# for information on the fitdistrplus package.
# Refer to AICcmodavg package for details on creating the AIC tables and
# how to do model averaging.

# Load the necessary libraries
library(AICcmodavg)
library(actuar) # for the burr distribution
library(FAdist) # for log-logistic
library(fitdistrplus)
library(ggplot2)
library(grid)
library(gridGraphics)
library(plyr)   # for doing a separate analysi by chemical
library(VGAM) # for the gompertz, pareto, gumbel distributions



source("mygofstat.r")  # corrects and error in the code in the fistdistr package (groan)
#source("C:/R-repositories/SSD-fitting/CCME-fits/mygofstat.r")

#load test data
endpoints<-read.csv("CCME data.csv",header=TRUE,as.is=TRUE, strip.white=TRUE)
#endpoints<- read_csv("C:/R-repositories/SSD-fitting/CCME-fits/CCME data.csv")

head(endpoints)

xtabs(~Chemical, data=endpoints, exclude=NULL, na.action=na.pass)

# create a directory for plots
dir.create(file.path("Plots"))


#use descdist function from fitdistrplus package to identify candidate distributions - create a Cullen Frey graph
sumstats <- plyr::ddply(endpoints, "Chemical", function(x, boot){
    # create the individual plots
    png(file.path("Plots",paste(x$Chemical,".png",sep="")), h=6, w=6, units="in", res=300)
       res <- descdist(x$Conc,discrete=FALSE, boot=boot)
    dev.off()
    # return the summary statistics as a dataframe, rather than a list
    attr(res, "class") <- NULL # remove class attribute from results
    data.frame(res)
}, boot=100)

# summary statistics
sumstats

# create list of possible distributions (the pareto and burr distributions have been commented out of the list below because they are not included in SSD Master)
fit.list <- list( lnorm  = list(data=NULL, distr="lnorm",  method="mle"),
                  llog   = list(data=NULL, distr="llog",   method="mle"), # log logistic
                  gomp   = list(data=NULL, distr="gompertz",method="mle"),
                  lgumbel= list(data=NULL, distr="lgumbel",method="mle"), # log gumbel
                  gamma  = list(data=NULL, distr="gamma",  method="mle"),                         
#                  pareto = list(data=NULL, distr="pareto", method="mle"),
                  weibull= list(data=NULL, distr="weibull",method="mle")
#                  burr   = list(data=NULL, distr="burr",   method="mle")
                 )


# define log-gumbel distribution
# These functions are needed because this is a non-standard distribution
dlgumbel <- function(x, location=0, scale=0, log=FALSE){ 
            fx <- dgumbel(log(x), location=location, scale=scale, log=FALSE)/x
            if(log) fx <- log(fx)
            fx}
qlgumbel <- function(p, location=0, scale=0, lower.tail=TRUE, log.p=FALSE){ 
            if(log.p) p<- exp(p)
            if(!lower.tail) p <- 1-p
            exp( qgumbel(p, location=location, scale=scale))}
plgumbel <- function(q, location=0, scale=0, lower.tail=TRUE, log.p=FALSE){ 
             Fq <- pgumbel(log(q), location=location, scale=scale)
             if(!lower.tail) Fq <- 1-Fq
             if(log.p) Fq <- log(Fq)
             Fq}
rlgumbel <- function(n, location=0, scale=0){ exp(rgumbel(n, location=location, scale=scale))}

#endpoints <- endpoints[ endpoints$Chemical=="Boron_CCME",]  # for testing to fit only one chemical

# Fit each distribution from the list to each chemical
fit.all <- plyr::dlply(endpoints, "Chemical", function(x, probs=c(.05,.10), nboot=1000){
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
              main=paste("Empirical and theoretical CDF's - ",x[1,"Chemical"],sep=""),
              fitcol=c("red","blue","darkgreen","black","darkblue","cyan","violet"),
              fitlty=1:3)
      gridGraphics::grid.echo()
      cdf.comp.plot <- grid::grid.grab()
      
      #compute the goodmess of fit statistics - NOTE this does not work for burr distribution, not sure why...
      # must check to see if more than one distribution was fit - groan.
      # Also had to fix an error in compute p-value for chi2 test when df=0
      if(length(dist.fits) ==1) {gof.stat <- mygofstat(dist.fits[[1]], fitnames=dist.names)}
      if(length(dist.fits) >1 ) {gof.stat <- mygofstat(dist.fits,fitnames =dist.names)}
#      The following fails because of a bug in gofstat. See mygofstat file for details
#      if(length(dist.fits) ==1) {gof.stat <- fitdistrplus::gofstat(dist.fits[[1]], fitnames=dist.names)}
#      if(length(dist.fits) >1 ) {gof.stat <- fitdistrplus::gofstat(dist.fits,fitnames =dist.names)}
      
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

# save the comparative plot
plyr::l_ply(fit.all, function(x){
   file.name=file.path("Plots",paste(x$Chemical,"-comparative-plot.png",sep=""))
   png(file.name, h=6, w=6, units="in", res=300)
      grid.newpage()
      grid.draw(x$cdf.comp.plot)
   dev.off()
})

# Make a nice plot of the estimates and ci and compare them to the CCME guidelines

ccme.csv<- textConnection(
"Chemical,n, HC5, lcl, ucl, dummy, dummy, dummy
Boron  , 28 , 1.5  , 1.2 , 1.7 , 1.2  , 0.59 , 3.20
Cadmium  , 36 , 0.09  , 0.04 , 0.24 , 0.14  , 0.06 , 0.34
Chloride  , 28 , 120  , 90 , 150 , 73  , 27 , 198
Endosulfan  , 12 , 0.003  , 0.0007 , 0.01 , 0.010  , 0.0012 , 0.51
Glyphosate  , 18 , 800  , 490 , 1320 , 900  , 459 ,  2301
Silver  , 9 , 0.25  , 0.17 , 0.39 , 0.19  , 0.069 , 0.89
Uranium , 13 , 15  , 8.5 , 25 , 15  , 3.1 , 124")

ccme <- read.csv(ccme.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
ccme$dummy <- NULL
ccme$dummy.1 <- NULL
ccme$dummy.2 <- NULL
ccme$EstType <- "CCME"
ccme$n  <- NULL
ccme

# get the estimates of HC5
fit.hc5 <- plyr::ldply(fit.all, function(x){
    # extract the model averaged
    q <- x$q.modavg[ x$q.modavg$quantile == 0.05,]
    data.frame( HC5=q$avg.pred, lcl=q$avg.lcl, ucl=q$avg.ucl, EstType="ModelAvg")
})

fit.hc5$Chemical <- gsub("_CCME", "", fit.hc5$Chemical)
fit.hc5quuantile <- NULL
fit.hc5

both <- rbind(ccme, fit.hc5)
both

lab_log10 <- function (x){
   # create nice labels for the next plot
   y <- as.character(x)
   y
}


compplot <- ggplot(data=both, aes(x=Chemical, y=HC5, color=EstType))+
  geom_point( position=position_dodge(w=0.4))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.1,position=position_dodge(w=0.4))+
  scale_y_log10(labels=lab_log10)+
  #scale_color_discrete(name="Method")+
  labs(x = "Chemical")+
  ylab(bquote('Log' ~HC[5]))+
  theme(axis.text.x = element_text(angle=45,hjust = 1))+
  scale_color_manual(labels = c("Least Squares", "MLE"), values = c("red", "black"))+
  guides(colour=guide_legend("Method"))
compplot

ggsave("Plots/compplot_mod.png", compplot, width = 5.97, height = 5.97)


#calculating percent difference of HC5 values
both.hc5 <- reshape2::dcast(both, Chemical ~ EstType, value.var="HC5")
 
both.hc5$perdiff <- abs(both.hc5$CCME - both.hc5$ModelAvg)/(both.hc5$CCME + both.hc5$ModelAvg)*2*100
both.hc5
summarize(both.hc5, 
          mean.perdiff = mean(perdiff))


#calculate width ratio (ModelAvg/CCME)
both$ci.width <- both$ucl - both$lcl
both.hc5 <- reshape2::dcast(both, Chemical ~ EstType, value.var="ci.width")

both.hc5$width.ratio <- both.hc5$ModelAvg / both.hc5$CCME
both.hc5
summarize(both.hc5, 
          mean.widthratio = mean(width.ratio))


