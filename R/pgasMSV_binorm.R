################################################################################################# 
#     PGAS for Multivarite SV model.    
#        source this file for MSV model      
#    2 parameters:   X(t) = phi * X(t-1) + W(t);         W(t) ~ iid N(0,q)  
#                  Y_i(t) = beta_i*exp{X(t)/2}V_i(t);  V_i(t) ~ iid N(0,1)
#             see line 108 to change priors on betas        
##################################################################################################

pgasMSV_binorm = function(nmcmc, burnin, y, init, beta_init, hyper, sigma_MH, npart, mcmseed){
  # Input:
  #   nmcmc - number of MCMC
  #   burnin - number of burnin
  #   y - measurements
  #   init - initial value of (phi, q)
  #   beta_init - initial value of beta_i, e.g. user inputs (beta_1, beta_2, beta_3)
  #   npart - number of particles
  #   hyper - hyperparameters for bivariate normal (phi, q), user inputs (mu_phi, mu_q, sigma_phi, sigma_q, rho)
  #   sigma_MH - covariance matrix for Random Walk Metropolis Hastings
  #   mcmseed - seed for mcmc
  
  # Output:
  #   phi - sampled phi
  #   q - sampled phi
  #   beta - sampled beta_i
  #   X - sampled hidden states
  #   time - running time
  #   acp - acceptence rate of Random Walk Metropolis Hastings
  
  ptm <- proc.time()
  numMCMC = nmcmc+burnin
  N = npart
  set.seed(mcmseed)
  
  pr <- progress_text()         # displays progress (from plyr)
  pr$init(numMCMC) 
  
  T = nrow(y)
  p = ncol(y)
  X = matrix(0,numMCMC,T)
  q = rep(0,numMCMC)             # q = state variance ## W(t) ~ iid N(0,q)
  phi = rep(0,numMCMC)           # phi-values
  mu = rep(0,numMCMC)            # leverage term
  beta = matrix(0,numMCMC,p)
  

  # Initialize the parameters
  phi[1] = init[1]
  q[1] = init[2] 
  mu[1] = 0
  for(i in 1:p){
    beta[1,i] = beta_init[i]
  }
  
  # hyperparameters
  mu_phi = hyper[1]
  mu_q = hyper[2]
  sigma_phi = hyper[3]
  sigma_q = hyper[4]
  rho = hyper[5]  # negative correlation
  mu_MH = c(0,0)  # tuning parameter
  sigma_MH = sigma_MH # tuning parameter
  
  
  # Initialize the state by running a PF
  u = cpf_as_sv(y, phi[1], q[1], 0, N,  X[1,], beta[1,])   #  changed from X to X[1,]
  particles = u$x           # returned particles
  w = u$w                   # returned weights
  # Draw J
  J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
  X[1,] = scale(particles[J,], center=TRUE, scale=FALSE)        

  
  # Run MCMC loop
  for(k in 2:numMCMC){
    # Sample the parameters (phi, sqrt_q) ~ bivariate normal Random Walk Metropolis Hastings
    parms = cbind(phi[k-1], sqrt(q[k-1]))
    parms_star = parms + mvrnorm(1, mu_MH, sigma_MH)
    
    while(parms_star[2]^2 < 2e-2|parms_star[1]>2){
      parms_star = parms + mvrnorm(1, mu_MH, sigma_MH)
    }
    
    g_star = log_g_func(parms_star, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], mu[k-1])
    g_old = log_g_func(parms, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], mu[k-1])
    g_diff = g_star - g_old
    if(is.na(g_diff)){
      phi[k] = parms[1]
      q[k] = parms[2]^2
    }else {
      if(log(runif(1)) < g_diff){
        phi[k] = parms_star[1]
        q[k] = parms_star[2]^2
      }else{
        phi[k] = parms[1]
        q[k] = parms[2]^2
      }
    }
    
    mu[k] = 0
    
    
    # Sample beta, prior unif, N(0,b)=N(0,inf)
    #for (i in 1:p){
    #  atemp = T/2 - 1
    #  btemp = sum(y[,i]^2/exp(X[k-1,]))/2
    #  beta[k,i] = sqrt(1/rgamma(1, atemp, btemp))
    #}
    
     # Sample beta, prior IG(a,b)
     a = .001
     b = .001
     for (i in 1:p){
       atemp = T/2 + a
       btemp = sum(y[,i]^2/exp(X[k-1,]))/2 + b
       beta[k,i] = sqrt(1/rgamma(1, atemp, btemp))
     }
    
    # Run CPF-AS
    u = cpf_as_sv(y, phi[k], q[k], 0, N, X[k-1,], beta[k,])
    particles = u$x           # returned particles
    w = u$w                   # returned weight
    # Draw J (extract a particle trajectory)
    J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
    X[k,] = scale(particles[J,], center=TRUE, scale=FALSE)        
    pr$step()
  }#end
  
  bi = 1:burnin  
  time2run = proc.time() - ptm
  acp = dim(table(phi[-bi]))/nmcmc*100 
  list(phi=phi[-bi], q=q[-bi], beta=beta[-bi,], X=X[-bi,], time=time2run, acp=acp)
}#end


#--------------------------------------------------------------------------
cpf_as_sv = function(y, phi, q, mu, N, X, beta){
  # Conditional particle filter with ancestor sampling
  # Input:
  #   y - measurements
  #   phi - transition parameter
  #   q - state noise variance
  #   mu - leverage term
  #   N - number of particles
  #   X - conditioned particles
  
  T = nrow(y)
  pp= ncol(y)
  x = matrix(0, N, T); # Particles
  a = matrix(0, N, T); # Ancestor indices
  w = matrix(0, N, T); # Weights
  x[,1] = 0; # Deterministic initial condition
  x[N,1] = X[1] 
  
  for (t in 1:T){
    if(t != 1){
      ind = resamplew(w[,t-1]);             
      ind = ind[sample.int(N)]; # default: no replacement
      t1 = t-1
      xpred = phi*x[, t1] # length N Vector
      x[,t] = xpred[ind] + sqrt(q)*rnorm(N);
      x[N,t] = X[t]; #Line 6
      
      # Ancestor sampling ---# Eq(3) 
      
      logm = -1/(2*q)*(X[t]-xpred)^2;
      maxlogm = max(log(w[,t-1])+logm)
      w_as = exp(log(w[,t-1])+logm - maxlogm)
      w_as = w_as/sum(w_as);
      
      ind[N] =  which( (runif(1)-cumsum(w_as)) < 0 )[1]   
      
      # Store the ancestor indices
      a[,t] = ind;
    }#end IF
    # Compute importance weights
    exp_now = exp(x[,t])
    temp = (y[t,]/beta)^2
    logweights = -sum(log(beta)) - (pp/2)*x[,t] - 0.5*(sum(temp)/exp_now) # (up to an additive constant)
    const = max(logweights); # Subtract the maximum value for numerical stability
    weights = exp(logweights-const)
    w[,t] = weights/sum(weights) # Save the normalized weights
  }#end FOR
  
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
  u = runif(N)    # uniform random number
  cw = cumsum(w)  # Cumulative Sum
  cw = cw/cw[N]
  ucw = c(u,cw)
  ind1 = sort(ucw, index.return=TRUE)$ix   # $ix will give the index
  ind2 = which(ind1<=N)
  i = ind2-(0:(N-1))
  return(i)
}

#-------------------------------------------------------------------
log_g_func = function(parms, mu_phi, sigma_phi, mu_q, sigma_q, rho, x, mu){
  # Calculate acceptance probability for RWMH
  phi = parms[1]
  q = parms[2]
  Z = cbind(x[2:T]-mu, x[1:(T-1)]-mu)  # (T-1)*2 matrix
  V = t(Z)%*%Z
  term1 = -((phi-mu_phi)^2/(sigma_phi^2) + (q-mu_q)^2/(sigma_q^2) - 2*rho*(phi-mu_phi)*(q-mu_q)/(sigma_q*sigma_phi))/(2*(1-rho^2))
  term2 = -(1-phi^2)*(x[1]-mu)^2/(2*q^2) - V[1,1]/(2*q^2) - V[2,2]*(phi^2)/(2*q^2) + V[1,2]*phi/(q^2)
  g = .25*(log((1-phi^2)^2)) - .5*T*log((q)^2) + term1 + term2
  return(g)
}

