################################################################################################# 
#   PGAS for simple SV model.  Input (nmcmc, burnin, y, prior, phi_init, q_init, npart, mcmseed)    
#        source this file then use "Run PGAS for SV model" code       
#    2 parameters: X(t)=phi X(t-1) + W(t) ~ iid N(0,q) ;  Y(t)=exp{X(t)/2)}V(t) ~ iid N(0,1)
##################################################################################################
#
#  This code was provided by Fredrik Lindsten https://liu.se/en/employee/freli29
#
##################################################################################################

pgasSV = function(nmcmc, burnin, y, prior, phi_init, q_init, npart, mcmseed){

numMCMC = nmcmc+burnin
N = npart
set.seed(mcmseed)

pr <- progress_text()         # displays progress (from plyr)
pr$init(numMCMC) 

T = length(y)
q = rep(0,numMCMC)             # q = state variance
X = matrix(0,numMCMC,T)
phi = rep(0,numMCMC)           #  phi-values

# Initialize the parameters
q[1] = q_init 
phi[1] = phi_init

# Initialize the state by running a PF
u = cpf_as_sv(y, phi[1], q[1], N,  X[1,])   #  changed from X to X[1,]
   particles = u$x           # returned particles
   w = u$w                   # returned weights
# Draw J
 J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
 X[1,] = particles[J,]      
 
   
# Run MCMC loop
 for(k in 2:numMCMC){
    # Sample the parameters (inverse gamma posteriors)
     Z = cbind(X[k-1,2:T], X[k-1,1:(T-1)])
     V = prior*diag(1,2) +   t(Z)%*%Z      
     m = V[2,1]/V[2,2]
     L = V[2,2]
     S = V[1,1] - m*V[1,2]
    # Sample from NIG posterior
     q[k] = 1/rgamma(1, (T-1)/2, S/2)
     phi[k] = m + sqrt(q[k]/L)*rnorm(1)
    # Run CPF-AS
      u = cpf_as_sv(y, phi[k], q[k], N, X[k-1,])
      particles = u$x   # returned particles
      w = u$w                   # returned weight
    # Draw J (extract a particle trajectory)
        J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
        X[k,] = particles[J,]     # center X
  pr$step()
}#end
bi = 1:burnin 
list(phi=phi[-bi], q=q[-bi],  X=X[-bi,])
}#end





#--------------------------------------------------------------------------
cpf_as_sv = function(y, phi, q, N, X){
# Conditional particle filter with ancestor sampling
# Input:
#   y - measurements
#   phi - transition parameter
#   q - state noise variance
#   N - number of particles
#   X - conditioned particles 

T = length(y)
x = matrix(0, N, T); # Particles
a = matrix(0, N, T); # Ancestor indices
w = matrix(0, N, T); # Weights
x[,1] = 0; # Deterministic initial condition
x[N,1] = X[1] 

for (t in 1:T){
    if(t != 1){
        ind = resamplew(w[,t-1]);             
        ind = ind[sample.int(N)];
          t1 = t-1
          xpred = phi*x[, t1] 
          x[,t] = xpred[ind] + sqrt(q)*rnorm(N);
          x[N,t] = X[t]; 
            # Ancestor sampling                                              
            m = exp(-1/(2*q)*(X[t]-xpred)^2);
            w_as = w[,t-1]*m
            w_as = w_as/sum(w_as);
            ind[N] =  which( (runif(1)-cumsum(w_as)) < 0 )[1]       
        # Store the ancestor indices
        a[,t] = ind;
    }#end
    # Compute importance weights
	var_now = exp(x[,t])
    logweights = -y[t]^2/(2*var_now) - 1/2*log(var_now)  # (up to an additive constant)
    const = max(logweights); # Subtract the maximum value for numerical stability
    weights = exp(logweights-const)
    w[,t] = weights/sum(weights) # Save the normalized weights
}#end

# Generate the trajectories from ancestor indices
ind = a[,T];
for(t in (T-1):1){
    x[,t] = x[ind,t];
    ind = a[ind,t];
}#end
list(x=x, w=w)
}#end
#-------------------------------------------------------------------

resamplew = function(w){
# multinomial resampling
N = length(w)
u = runif(N)
cw = cumsum(w)
cw = cw/cw[N]
 ucw = c(u,cw)
 ind1 = sort(ucw, index.return=TRUE)$ix
 ind2 = which(ind1<=N)
i = ind2-(0:(N-1))
return(i)
}
