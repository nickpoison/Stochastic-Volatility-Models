##########################################################################
##            MSV Model using PGAS joint modeling not adaptive          ##
##               for banks BOA, Citi and JPM                            ##
##                                                                      ##
##            we pay our debts on time ... sometimes                    ##
##########################################################################


library(astsa)
library(plyr)
library(MASS) 
library(mcmc)
source("R/pgasMSV_binorm.R")

# Load data
load("data/BCJ.rda")
y   = 100*BCJ 


#-------------------------------------------------------------------
npart   = 20       # Number of particles used in pgas                        
nmcmc   = 2000     # Number of iterations in the mcmc sampler after burnin       
burnin  = 500      # Number of iterations to burn                               
mcmseed = 90210

#  Initial values of phi, q, mu 
init      = c(0.9, 0.25)
beta_init = c(1.5, 1.5, 1.5)
hyper     = c(0.9, 0.25, 0.075, 0.1, -0.25)
sigma_MH  = .03 * matrix(c(1,-.25,-.25,1),nrow=2,ncol=2)

# Run it 
ptm <- proc.time()
u_nd = pgasMSV_binorm(nmcmc, burnin, y, init, beta_init, hyper, sigma_MH, npart, mcmseed)
(time2run = proc.time() - ptm)

# Acceptance Rate
cat("The acceptance rate is", u_nd$acp, "%.")



################  Results   #######################
### Pretty pictures
parms = cbind(u_nd$phi, sqrt(u_nd$q), u_nd$beta[,1], u_nd$beta[,2], u_nd$beta[,3])
names = c(expression(phi), expression(sigma), expression(beta[1]), expression(beta[2]), expression(beta[3]))


culer = c(rgb(.66,.12,.85), rgb(.12,.66,.85), rgb(.8*.52,.8*.87,.8*.08), rgb(.83,.31,.10), rgb(.63,.11, 0.40))
culerb = c(rgb(.66,.12,.85, .4), rgb(.12,.66,.85,.4), rgb(.8*.52,.8*.87,.8*.08,.4), rgb(.83,.31,.10,.4), rgb(.63,.11, 0.40,.4) )

###### parameters ##########
dev.new(height=7, width=9)  
IF = sprintf("%.2f", 
round(cbind(initseq(parms[,1])$var.pos/initseq(parms[,1])$gamma0, 
            initseq(parms[,2])$var.pos/initseq(parms[,2])$gamma0, 
			initseq(parms[,3])$var.pos/initseq(parms[,3])$gamma0, 
            initseq(parms[,4])$var.pos/initseq(parms[,4])$gamma0, 
			initseq(parms[,5])$var.pos/initseq(parms[,5])$gamma0), 3 ))
par(mfrow=c(3,5))
for (i in 1:5){
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
  acf1(parms[,4], 100, ylim=c(L,U), main='', col=culer[4])
    legend('topright', legend=paste('IF =', IF[4]),  bty='n', cex=1.25)
  acf1(parms[,5], 100, ylim=c(L,U), main='', col=culer[5])
    legend('topright', legend=paste('IF =', IF[5]),  bty='n', cex=1.25)	
 for (i in 1:5){
   hist(parms[,i], prob=TRUE, breaks=15, main='', xlab='')
   Grid(nx=0, ny=NULL)
   hist(parms[,i], breaks=15, col=culerb[i], border=culer[i], prob=TRUE, add=TRUE)
  lines(density(parms[,i], adjust=3))
  abline(v=mean(parms[,i]), col=rgb(0,0,0,.8), lwd=2, lty=1)
} 


#####  data and states ############

tspar =  tsp(y)
mX = ts(apply(u_nd$X, 2, mean), start=tspar[1], frequency=tspar[3])
lX = ts(apply(u_nd$X, 2, quantile, 0.005), start=tspar[1], frequency=tspar[3])
uX = ts(apply(u_nd$X, 2, quantile, 0.995), start=tspar[1], frequency=tspar[3])

#X = window(X, start=2005.008) 
dev.new(width=10, height=8) 
par(mfrow=c(4,1), cex.lab=1.25)
tsplot(y[,1], ylab='BOA', col=culer[1], lwd=2)
tsplot(y[,2], ylab='Citi', col=culer[2], lwd=2)
tsplot(y[,3], ylab='JPM', col=culer[3], lwd=2)
tsplot(mX, col=culer[4], ylim=c(-6,6), lwd=1, ylab='Log Volatiltiy (%)')
  xx=c(time(mX), rev(time(mX)))
  yy=c(lX, rev(uX))
polygon(xx, yy, border=NA, col=culerb[4]) 
lines(lowess(mX, f=.04, delta=.001), col=rgb(0,0,0,.25))	  






















