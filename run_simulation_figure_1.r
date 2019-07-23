library(astsa)
#
set.seed(1989)
num  = 1000
lev  =  0
beta =  1
Vsd  =  1
phi1 =  0.99  # Model I
sig1 =  0.15  # Model I
phi2 =  0.95  # Model II
sig2 =  0.35  # Model II
x1   =  arima.sim(list(order=c(1,0,0), ar=phi1), sd=sig1, n=num) + lev   
y1   =  beta*exp(x1/2)*rnorm(num,0,Vsd)  
x2   =  arima.sim(list(order=c(1,0,0), ar=phi2), sd=sig2, n=num) + lev   
y2   =  beta*exp(x2/2)*rnorm(num,0,Vsd) 
y1   =  as.vector(y1)  
y2   =  as.vector(y2)
L    =  min(y1,y2)  
U    =  max(y1,y2)

SVacf = function(phi,sig,h){
            # returns the acf of the squared returns
			sigx2 = sig^2/(1-phi^2)                
            numer = exp(sigx2*phi^h)-1
			denom = 3*exp(sigx2)-1
			return(numer/denom)
}

# ACFs
u1 = SVacf(phi1, sig1, 1:30)
u2 = SVacf(phi2, sig2, 1:30)


##
culera = c(rgb(.85,.30,.12,.6), rgb(.12,.65,.85,.6))
culer  = c(rgb(.85,.30,.12), rgb(.12,.65,.85))
 
dev.new(width=9, height=6)
layout(matrix(c(1,2,3, 1,2,4), nrow=3), heights=c(1,1,1.25))
tsplot(y1, col=culera[1], lwd=2,   ylab='%', main='Series A', col.main=culer[1])
tsplot(y2, col=culera[2], lwd=2,   ylab='%', main='Series B', col.main=culer[2])
acf1(y1^2, 30, main='', ylim=c(0,.6), col=culer[1], lwd=4)
 title('Squared Series A', col.main=culer[1])
 lines(u1, lwd=2)
 lines(u2, lwd=2, lty=2)
 op <- par(family = "serif")
 legend('topright', lty = 1:2, legend=c('Model I', 'Model II'), lwd=2, cex=1.1, bg='white')
 par(op)
acf1(y2^2, 30, main='', ylim=c(0,.6), col=culer[2], lwd=4)
 title('Squared Series B', col.main=culer[2])
 lines(u1, lwd=2)
 lines(u2, lwd=2, lty=2)
 op <- par(family = "serif")
 legend('topright', lty = 1:2, legend=c('Model I', 'Model II'), cex=1.1, bg='white')
 par(op)

