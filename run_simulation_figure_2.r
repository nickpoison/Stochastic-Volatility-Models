##########################################
##                                      ##
##        2 parameter model             ##
##    the old crappy method that        ##  
##  samples parameters one at a time    ##
##                                      ##
##########################################

library(astsa)
library(plyr)
library(MASS)
source('R/pgasSV.r')


npart   = 20       # Number of particles used in pgas                        
nmcmc   = 1000     # Number of iterations in the mcmc samplers after burnin       
burnin  = 200      # Number of iterations to burn                               
mcmseed = 90210

#
set.seed(1989)
num  = 1000
mu   =  0
beta =  1
Vsd  =  1
phi2 =  0.95  # Model II
sig2 =  0.35  # Model II
x2   =  arima.sim(list(order=c(1,0,0), ar=phi2), sd=sig2, n=num) + mu    
y2   =  beta*exp(x2/2)*rnorm(num,0,Vsd) 
y    =  as.vector(y2)

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


culer   =  c(rgb(.85,.30,.12),    rgb(.12,.65,.85))

dev.new(width=9, height=6) 
layout(matrix(c(1,3, 2,4), 2), heights=c(1,1))
for (i in 1:2){
 tsplot(parms[,i], ylab='trace', xlab='index', ,main='', col=culer[i])
 mtext(names[i], 3, line=.25, cex=1)
}
u1 = acf1(parms[,1], 300, plot=FALSE)
u2 = acf1(parms[,2], 300, plot=FALSE)
tsplot(u1, col=culer[1], lwd=2, ylab='ACF', xlab='LAG', ylim=c(-.2,1), lty=6)
 lines(u2,  col=culer[2], lwd=2)
 legend('topright', lty = c(6,1), legend=names, col=culer, cex=1.1, bg='white', lwd=2, text.col=culer)
 abline(h=c(0,-2/sqrt(1000),2/sqrt(1000)), col=gray(.5), lty=c(1,2,2))
z  = kde2d(parms[,1], parms[,2], h=c(.075,.15))
plot(parms[,1], parms[,2], pch=20, col=rgb(24,116,205,max=255,alpha=150), 
     xlab=expression(phi), ylab=expression(sigma), cex=1.25, cex.lab=1.25, panel.first=Grid())
 contour(z, drawlabels=FALSE, nlevels=15,  add=TRUE, col=gray(.5))
 abline(h=mean(parms[,2]), col=gray(.5))
 abline(v=mean(parms[,1]), col=gray(.5))

 