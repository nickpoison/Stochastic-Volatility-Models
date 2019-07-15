##########################################
##                                      ##
##    2 parameter model                 ##
##  sample parameters one at a time     ##
##########################################

library(astsa)
library(plyr)
library(MASS)
source('R/pgasSV.r')


npart   = 20       # Number of particles used in pgas                        
nmcmc   = 1000     # Number of iterations in the mcmc samplers after burnin       
burnin  = 200      # Number of iterations to burn                               
mcmseed = 90210

set.seed(9753)
num   =  1000
lev   =  0
beta  = .1
Vsd   =  1 
x  =  arima.sim(list(order=c(1,0,0), ar=.92), sd=1.5, n=num) + lev
y  =  beta*exp(x/2)*rnorm(num,0,Vsd)

#  Improper prior and initial values of phi and q
prior    =  1  #  pi(phi, q) \propto  1/q  (Jeffreys prior)
phi_init = .9 
q_init   =  1  # q is sigma^2 


# Run it
 ptm <- proc.time()
u = pgasSV(nmcmc, burnin,  y, prior, phi_init, q_init, npart, mcmseed)
 (time2run = proc.time() - ptm)

########################################
# Pretty pictures
parms = cbind(u$phi, sqrt(u$q))
names = c(expression(phi), expression(sigma))
nobs  = length(y)

culera  =  c(rgb(.85,.30,.12,.7), rgb(.12,.65,.85,.7))
culer   =  c(rgb(.85,.30,.12),    rgb(.12,.65,.85))

dev.new(width=9, height=6) 
layout(matrix(c(1,3, 2,4), 2), heights=c(1,1))
for (i in 1:2){
 tsplot(parms[,i], ylab='trace', xlab='index', ,main='', col=culera[i])
 lines(lowess(parms[,i], f=.05), lwd=2, col=culera[i%%2+1])
 # abline(h=mean(parms[,i]), col=culera[i%%2+1], lwd=2)
 mtext(names[i], 3, line=.25, cex=1)
}
u1 = acf1(parms[,1], 300, plot=FALSE)
u2 = acf1(parms[,2], 300, plot=FALSE)
tsplot(u1, col=culera[1], lwd=2, ylab='ACF', xlab='LAG', ylim=c(-.1,.8))
lines(u2,  col=culera[2], lwd=2)
abline(h=c(0,-sqrt(2/1000),sqrt(2/1000)), col=gray(.5), lty=c(1,2,2))
z  = kde2d(parms[,1], parms[,2], h=c(.01,.15))
plot(parms[,1], parms[,2], pch=20, col=rgb(24,116,205,max=255,alpha=150), 
     xlab=expression(phi), ylab=expression(sigma), cex=1.25, cex.lab=1.25, panel.first=Grid())
contour(z, drawlabels=FALSE, nlevels=15,  add=TRUE, col=gray(.5))
abline(h=mean(parms[,2]), col=gray(.5))
abline(v=mean(parms[,1]), col=gray(.5))

 