##############################################
## SP500 SV Model using PGAS_adaptive       ##
##    2 parameters model                    ##
##  this is using the new method            ##
##############################################


library(astsa)
source("R/pgasSV_binorm.R")
library(plyr)
library(MASS) 
library(mcmc)
load('data/sp500.gr.rda')

y = 100*window(sp500.gr, start=2005) 


#-------------------------------------------------------------------
npart   = 20       # Number of particles used in pgas 
nmcmc   = 5000     # Number of iterations in the mcmc samplers after burnin       
burnin  = 100      # Number of iterations to burn                               
mcmseed = 90210 

#  Bi-normal prior and initial values
init        =  c(.9, 0.35, 0) 
hyper       =  c(0.9, 0.35, 0.125, 0.25, -0.6)
sigma_MH    = .1*matrix(c(1,-0.25,-0.25,1),nrow=2,ncol=2)


# Run it
 ptm <- proc.time()
uu = pgasSV_binorm(nmcmc, burnin, y, init, hyper, sigma_MH, npart,  parms_2 = T, mcmseed)
 (time2run = proc.time() - ptm)

# Acceptance Rate
cat("The acceptance rate is", uu$acp, "%.")



############### Results  ########################
culerb = c(rgb(.85,.30,.12,.3), rgb(.12,.65,.85,.3))
culera = c(rgb(.85,.30,.12,.6), rgb(.12,.65,.85,.7))
culer  = c(rgb(.85,.30,.12),    rgb(.12,.65,.85))

### Pretty pictures
parms = cbind(uu$phi, sqrt(uu$q))
names = c(expression(phi), expression(sigma))
nobs  = length(y)

tspar =  tsp(sp500.gr)
mX = ts(apply(uu$X, 2, mean), start=tspar[1], frequency=tspar[3])
lX = ts(apply(uu$X, 2, quantile, 0.025), start=tspar[1], frequency=tspar[3])
uX = ts(apply(uu$X, 2, quantile, 0.975), start=tspar[1], frequency=tspar[3])

# state and data
dev.new(width=9, height=6) 
par(mfrow=c(2,1))
tsplot(y, ylab='S&P500 Returns (%)', col='dodgerblue3', lwd=2)
tsplot(mX, col=culera[2], ylim=c(min(lX)-.1, max(uX)+.1), lwd=2, ylab='Log Volatiltiy', )
  xx=c(time(mX), rev(time(mX)))
  yy=c(lX, rev(uX))
polygon(xx, yy, border=NA, col=culerb[2]) 	   


# parameters
dev.new(width=9, height=6)  
IF = sprintf("%.2f", 
     round(cbind(initseq(parms[,1])$var.pos/initseq(parms[,1])$gamma0, 
	             initseq(parms[,2])$var.pos/initseq(parms[,2])$gamma0), 2))
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
# legend('topright', legend=c(IF[1], IF[2]), lwd=2, col=culer, title.col=1, text.col=culer,  title="Inefficiency", bg='white', text.font=2)
 legend('topright', legend=c(as.expression(bquote(~~.(IF[1])~~"["~phi~"]")), 
		                     as.expression(bquote(.(IF[2])~~"["~sigma~"]"))), 
        lwd=2, col=culer, title.col=1, lty=c(6,1), text.col=culer, title="Inefficiency", bg='white', text.font=2)
plot(parms[,1], parms[,2], pch=20, col=rgb(24,116,205,max=255,alpha=150), xlab=expression(phi), 
      ylab=expression(sigma), cex=1.25, cex.lab=1.25, panel.first=Grid() )
abline(h=mean(parms[,2]), col=gray(.5))
abline(v=mean(parms[,1]), col=gray(.5))

# > cor(parms)
#            [,1]       [,2]
# [1,]  1.0000000 -0.5451042
# [2,] -0.5451042  1.0000000
# > colMeans(parms)
# [1] 0.8483645 0.4533468



############### data and states from new and old
tspar =  tsp(sp500.gr) 
# old stuff
mXo = ts(apply(u$X, 2, mean), start=tspar[1], frequency=tspar[3])
lXo = ts(apply(u$X, 2, quantile, 0.025), start=tspar[1], frequency=tspar[3])
uXo = ts(apply(u$X, 2, quantile, 0.975), start=tspar[1], frequency=tspar[3])
#####

# pdf(file="sp500_states_data.pdf", width=9, height=6) 
dev.new(width=9, height=6) 
culer = c(rgb(.85,.30,.12), rgb(.12,.65,.85))
culerb = c(rgb(.85,.30,.12,.3), rgb(.12,.65,.85,.3))
par(mfrow=c(3,1))
tsplot(y, ylab='S&P500 Returns (%)', col='dodgerblue3', lwd=2)
tsplot(mXo, col=culer[1], ylim=c(min(lXo)-.1, max(uXo)+.1), lwd=2, ylab='Log Volatiltiy', )
  xx=c(time(mXo), rev(time(mXo)))
  yy=c(lX, rev(uXo))
polygon(xx, yy, border=NA, col=culerb[1]) 	
tsplot(mX, col=culer[2], ylim=c(min(lX)-.1, max(uX)+.1), lwd=2, ylab='Log Volatiltiy', )
  xx=c(time(mX), rev(time(mX)))
  yy=c(lX, rev(uX))
polygon(xx, yy, border=NA, col=culerb[2]) 	



################## remap parms to causal #################
# parme = parms
# nn = nrow(parms)
# for (i in 1:nn){ 
#   if (parms[i,1]>1) {parme[i,2]= parms[i,2]/parms[i,1]; parme[i,1]=1/parms[i,1] }
#   }



