#########################################################################
##              Application on S&P500 using 2-parm SV Model            ##
##                  2 parameter model (set mu=0)                       ##
##                  this is old prior method                           ##
#########################################################################


library(astsa)
source("R/pgasSV.R")
library(plyr)
library(MASS) 
library(mcmc)
load('data/sp500.gr.rda')

y = 100*window(sp500.gr, start=2005) 



#-------------------------------------------------------------------
npart   = 20        # Number of particles used in pgas
nmcmc   = 5000      # Number of iterations in the mcmc samplers after burnin       
burnin  = 100       # Number of iterations to burn                               
mcmseed = 90210


prior    =  20   #  pi(phi, q) \propto  1/q  (Jeffreys prior)
phi_init = .8 
q_init   = .05   # q is sigma^2 

ptm <- proc.time()
u = pgasSV(nmcmc, burnin,  y, prior, phi_init, q_init, npart, mcmseed)
(time2run = proc.time() - ptm)




############### Results  ########################
culerb = c(rgb(.85,.30,.12,.3), rgb(.12,.65,.85,.3))
culera = c(rgb(.85,.30,.12,.6), rgb(.12,.65,.85,.7))
culer  = c(rgb(.85,.30,.12), rgb(.12,.65,.85))

### Pretty pictures
parms = cbind(u$phi, sqrt(u$q))
names = c(expression(phi), expression(sigma))
nobs  = length(y)

tspar =  tsp(sp500.gr)
mX = ts(apply(u$X, 2, mean), start=tspar[1], frequency=tspar[3]) 
lX = ts(apply(u$X, 2, quantile, 0.025), start=tspar[1], frequency=tspar[3]) 
uX = ts(apply(u$X, 2, quantile, 0.975), start=tspar[1], frequency=tspar[3]) 

# state and data
dev.new(width=9, height=6) 
par(mfrow=c(2,1))
tsplot(y, ylab='S&P500 Returns (%)', col='dodgerblue3', lwd=2)
tsplot(mX, col=culera[1], ylim=c(min(lX)-.1, max(uX)+.1), lwd=2, ylab='Log Volatiltiy', )
  xx=c(time(mX), rev(time(mX)))
  yy=c(lX, rev(uX))
polygon(xx, yy, border=NA, col=culerb[1]) 	   


# parameters
dev.new(width=9, height=6)  
IF = sprintf("%.2f", round(cbind(initseq(parms[,1])$var.pos/initseq(parms[,1])$gamma0, initseq(parms[,2])$var.pos/initseq(parms[,2])$gamma0),2))
layout(matrix(c(1,3, 2,4), 2), heights=c(1,1))
for (i in 1:2){
tsplot(parms[,i], ylab='trace', xlab='index', ,main='', col=culer[i])
 abline(h=mean(parms[,i]), col=culera[i%%2+1], lwd=2)
 mtext(names[i], 3, line=.25, cex=1)
 }
  u1 = acf1(parms[,1], 100, plot=FALSE)
  u2 = acf1(parms[,2], 100, plot=FALSE)
  L  = min(u1,u2)
  U  = max(c(u1,u2,1))
 tsplot(u2, col=culera[2], lwd=2, ylab='ACF', xlab='LAG', ylim=c(L,U))
 lines(u1, col=culera[1], lwd=2)
 abline(h=0)
 legend('topright', legend=c(IF[1], IF[2]), lwd=2, col=culer, title.col=1, text.col=culer, title="Inefficiency", bg='white', text.font=2)
plot(parms[,1], parms[,2], pch=20, col=rgb(24,116,205,max=255,alpha=150), xlab=expression(phi), ylab=expression(sigma), cex=1.25, cex.lab=1.25, panel.first=Grid() )
abline(h=mean(parms[,2]), col=gray(.5))
abline(v=mean(parms[,1]), col=gray(.5))


