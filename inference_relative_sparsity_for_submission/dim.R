library(latex2exp)
library(numDeriv)
library(Deriv)
library(xtable)
source('get.var.R')
source('grad.R')
source('utils.R')
source('sim.R')
source('center.scale.R')
source('is.R') 
#bottom, left, top, and right.
#par(mar=c(4.1, 4.1, 4.1, 2.1))
dontestb=0
#system("cat dim.R")


#ns=seq(1,1e5,length.out=100)
#delta=1
#par(mfrow=c(2,1))
#plot(ns,lambdan(ns,delta)/sqrt(ns)) #approx 0.005
#plot(ns,lambdan(ns,delta)*ns^((delta-1)/2))# like 27...

# good exp
# n=1e3
# T=2
# M=500

