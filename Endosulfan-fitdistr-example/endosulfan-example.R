# Example of fitting SSD from the fitdistrplus package
# This example has no censoring and is the simplest of the set

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



library(fitdistrplus)
library(ggplot2)

# (1) load of data and how the first few records
#
data(endosulfan)
head(endosulfan)
xtabs(~Australian+group, data=endosulfan)

# (2) plot and description of data for non Australian fish in base 10 logarithm
#
log10ATV <-log10(subset(endosulfan,(Australian == "no") & (group == "Fish"))$ATV)

plotdist(log10ATV)
descdist(log10ATV,boot=1000)

# (3) fit of a normal and a logistic distribution to data in log10
# (classical distributions used for SSD)
# and visual comparison of the fits
#
fln <- fitdist(log10ATV,"norm")
summary(fln)

fll <- fitdist(log10ATV,"logis")
summary(fll)

cdfcomp(list(fln,fll),legendtext=c("normal","logistic"),
xlab="log10ATV")

denscomp(list(fln,fll),legendtext=c("normal","logistic"),
xlab="log10ATV")

qqcomp(list(fln,fll),legendtext=c("normal","logistic"))
ppcomp(list(fln,fll),legendtext=c("normal","logistic"))

gofstat(list(fln,fll), fitnames = c("lognormal", "loglogistic"))

# (4) estimation of the 5 percent quantile value of 
# logistic fitted distribution (5 percent hazardous concentration  : HC5)
# with its two-sided 95 percent confidence interval calculated by 
# parametric bootstrap 
# with a small number of iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
#
# in log10(ATV)
bll <- bootdist(fll,niter=101)
HC5ll <- quantile(bll,probs = 0.05)
# in ATV
10^(HC5ll$quantiles)
10^(HC5ll$quantCI)

# (5) estimation of the 5 percent quantile value of 
# the fitted logistic distribution (5 percent hazardous concentration  : HC5)
# with its one-sided 95 percent confidence interval (type "greater")
# calculated by 
# nonparametric bootstrap 
# with a small number of iterations to satisfy CRAN running times constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.
# 
# in log10(ATV)
bllnonpar <- bootdist(fll,niter=101,bootmethod = "nonparam")
HC5llgreater <- quantile(bllnonpar,probs = 0.05, CI.type="greater")
# in ATV
10^(HC5llgreater$quantiles)
10^(HC5llgreater$quantCI)

# (6) fit of a logistic distribution 
# by minimizing the modified Anderson-Darling AD2L distance
# cf. ?mgedist for definition of this distance
#

fllAD2L <- fitdist(log10ATV,"logis",method="mge",gof="AD2L")
summary(fllAD2L)
plot(fllAD2L)

