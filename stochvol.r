library(stochvol)
library(astsa)
library(MASS)

# as in stochvol example
data(exrates)
dat <- logret(exrates$CHF, demean = TRUE)
rest <- svsample(dat, priormu = c(-12, 1), priorphi = c(20, 1.1), priorsigma = 0.1, priornu = c(2, 100), burnin = 2000)
#

parms = as.matrix(rest$para[,c('phi', 'sigma')])
names = c(expression(phi), expression(sigma))
culer =  c(rgb(.85,.30,.12),    rgb(.12,.65,.85))
nobs  = length(dat)

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
 abline(h=c(0,-2/sqrt(nobs),2/sqrt(nobs)), col=gray(.5), lty=c(1,2,2))
z  = kde2d(parms[,1], parms[,2], h=c(.008,.06))
plot(parms[,1], parms[,2], pch=20, col=rgb(24,116,205,max=255,alpha=100), 
    xlab=expression(phi), ylab=expression(sigma), cex=1.25, cex.lab=1.25, panel.first=Grid())
 contour(z, drawlabels=FALSE, nlevels=20,  add=TRUE, col=gray(.5))
 abline(h=mean(parms[,2]), col=gray(.5))
 abline(v=mean(parms[,1]), col=gray(.5))
 