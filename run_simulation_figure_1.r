library(astsa)
#
set.seed(9753)
num  = 1000
lev  =  0
beta = .1
Vsd  =  1
sig1 =  1
sig2 = 1.5
x1   =  arima.sim(list(order=c(1,0,0), ar=.97), sd=sig1, n=num) + lev   
y1   =  beta*exp(x1/2)*rnorm(num,0,Vsd)  
x2   =  arima.sim(list(order=c(1,0,0), ar=.92), sd=sig2, n=num) + lev   
y2   =  beta*exp(x2/2)*rnorm(num,0,Vsd) 
y1   =  as.vector(y1)  
y2   =  as.vector(y2)
L    =  min(y1,y2)  
U    =  max(y1,y2)

numer = function(sig,phi,h) { exp(sig*phi^h)-1 }
denom = function(sig) { 3*exp(sig)-1 }

# ACFs
phi =.97
u1  = numer(sig1,phi,1:30)/denom(sig1)
phi =.92
u2  = numer(sig2,phi,1:30)/denom(sig2)

##
culera = c(rgb(.85,.30,.12,.6), rgb(.12,.65,.85,.6))
culer  = c(rgb(.85,.30,.12), rgb(.12,.65,.85))
 
layout(matrix(c(1,2, 1,3), 2), heights=c(1,.75))
tsplot(y1, col=culera[1], lwd=2, ylim=c(L,U), ylab='%')
 lines(y2, col=culera[2], lwd=2)
 legend('bottomleft', col=culer, lwd=c(2,2), legend=c('Series A', 'Series B'), cex=1.1, bg='white')
acf1(y1^2, 30, main='', ylim=c(0,.6), col=culer[1], lwd=4)
 title('Squared Series A', col.main=culer[1])
 lines(u1, lwd=2)
 lines(u2, lwd=2)
acf1(y2^2, 30, main='', ylim=c(0,.6), col=culer[2], lwd=4)
 title('Squared Series B', col.main=culer[2])
 lines(u1, lwd=2)
 lines(u2, lwd=2)

