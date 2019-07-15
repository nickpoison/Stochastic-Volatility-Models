#########################################################################
##              Application on S&P500 using 3-parm SV Model           ##
##                  3 parameter model                                 ##
##              joint sampling not adaptive                           ##
#########################################################################


source("R/pgasSV_binorm.R")
library(astsa)
library(mcmc)
library(plyr)
library(MASS) 
load('data/sp500.gr.rda')

y = 100*window(sp500.gr, start=2005) 


#-------------------------------------------------------------------
npart   = 10       # Number of particles used in pgas
nmcmc   = 2000     # Number of iterations in the mcmc samplers after burnin       
burnin  = 100      # Number of iterations to burn                               
mcmseed = 90210

#  Bi-normal prior and initial values
init=c(0.85, 0.3, 0.15)
hyper = c(0.85, 0.3, 0.1, 0.1, -0.25)
sigma_MH = 0.05*matrix(c(1,0.25,0.25,1),nrow=2,ncol=2)
lambda_init = 2
alpha.tar = 0.3 

# Run it
ptm <- proc.time()
u3na = pgasSV_binorm(nmcmc, burnin, y, init, hyper, sigma_MH, npart, parms_2 = FALSE, mcmseed)
(time2run = proc.time() - ptm)

# Acceptance Rate
cat("The acceptance rate is", u3na$acp, "%.")


### Pretty pictures
parms = cbind(u3na$phi, sqrt(u3na$q), u3na$mu)
names = c(expression(phi), expression(sigma), expression(mu))
culer = c(rgb(.66,.12,.85), rgb(.12,.66,.85), rgb(.8*.52,.8*.87,.8*.08) )
culerb = c(rgb(.66,.12,.85, .4), rgb(.12,.66,.85,.4), rgb(.8*.52,.8*.87,.8*.08,.4) )

# parameters
dev.new(width=9, height=6)  
IF = sprintf("%.2f", round(cbind(initseq(parms[,1])$var.pos/initseq(parms[,1])$gamma0, initseq(parms[,2])$var.pos/initseq(parms[,2])$gamma0, initseq(parms[,3])$var.pos/initseq(parms[,3])$gamma0), 3 ))
par(mfrow=c(3,3))
for (i in 1:3){
tsplot(parms[,i], ylab='trace', xlab='index', ,main='', col=culer[i])
 abline(h=mean(parms[,i]), col=rgb(0,0,0,.5), lwd=1)
 mtext(names[i], 3, line=.25, cex=1)
 }
  L=-.1; U=1
  acf1(parms[,1], 100, ylim=c(L,U), main='', col=culer[1])
    legend('topright', legend=paste('IF =', IF[1]),  bty='n', cex=1.25)
  acf1(parms[,2], 100, ylim=c(L,U), main='', col=culer[2])
    legend('topright', legend=paste('IF =', IF[2]),  bty='n', cex=1.25)
  acf1(parms[,3], 100, ylim=c(L,U), main='', col=culer[3])
    legend('topright', legend=paste('IF =', IF[3]),  bty='n', cex=1.25)
 for (i in 1:2){
   hist(parms[,i], prob=TRUE, breaks=20, main='', xlab='')
   Grid(nx=0, ny=NULL)
   hist(parms[,i], breaks=20, col=culerb[i], border=culerb[i], prob=TRUE, add=TRUE)
  lines(density(parms[,i], adjust=3))
  abline(v=mean(parms[,i]), col=rgb(0,0,0,.8), lwd=2, lty=1)
} 
  hist(parms[,3], prob=TRUE, breaks=20, main='', xlab='', ylim=c(0,3.6), xlim=c(-1.5,1.5))
   Grid(nx=0, ny=NULL)
  hist(parms[,3], breaks=20, col=culerb[3], border=culerb[3], prob=TRUE, add=TRUE)
  lines(density(parms[,3], adjust=4))
  abline(v=mean(parms[,3]), col=rgb(0,0,0,.8), lwd=2, lty=1)
 


# states

tspar =  tsp(sp500.gr)
mX = ts(apply(u3na$X, 2, mean), start=tspar[1], frequency=tspar[3])
lX = ts(apply(u3na$X, 2, quantile, 0.025), start=tspar[1], frequency=tspar[3])
uX = ts(apply(u3na$X, 2, quantile, 0.975), start=tspar[1], frequency=tspar[3])
culera = c(rgb(.85,.30,.12,.3), rgb(.12,.65,.85,.3))
 
dev.new(width=9, height=6) 
par(mfrow=c(2,1))
tsplot(y, ylab='S&P500 Returns (%)', col='dodgerblue3', lwd=2)
tsplot(mX, col=culer[2], ylim=c(min(lX)-.1, max(uX)+.1), lwd=2, ylab='Log Volatiltiy (%)')
  xx=c(time(mX), rev(time(mX)))
  yy=c(lX, rev(uX))
polygon(xx, yy, border=NA, col=culera[2]) 	   
