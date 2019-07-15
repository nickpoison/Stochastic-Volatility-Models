####################################################################################################################### 
#   PGAS for SV model.   
#        source this file for SV model      
#   3 parameters: X(t) = mu + phi (X(t-1)-mu) + W(t) ~ iid N(0,q) ;  Y(t) = beta*exp{X(t)/2)}V(t) ~ iid N(0,1); beta = exp(mu/2)
#   2 parameters: X(t) = phi*X(t-1) + W(t) ~ iid N(0,q) ;  Y(t) = beta*exp{X(t)/2)}V(t) ~ iid N(0,1); beta = exp(mu/2)
#######################################################################################################################

pgasSV_binorm = function(nmcmc, burnin, y, init, hyper, sigma_MH, npart, parms_2 = TRUE, mcmseed){
  # Input:
  #   nmcmc - number of MCMC
  #   burnin - number of burnin
  #   y - measurements
  #   init - initial value of (phi, q, mu)
  #   npart - number of particles
  #   hyper - hyperparameters for bivariate normal (phi, q), user inputs (mu_phi, mu_q, sigma_phi, sigma_q, rho)
  #   sigma_MH - covariance matrix for Random Walk Metropolis Hastings
  #   parms_2 - IF (parms_2=TRUE) THEN (2 parms SV model) ELSE (3 parms SV model)
  #   mcmseed - seed for mcmc
  # Output:
  #   phi - sampled phi
  #   q - sampled phi
  #   mu - sampled mu
  #   X - sampled hidden states
  #   time - running time
  #   acp - acceptence rate of Random Walk Metropolis Hastings
  
  ptm <- proc.time()
  numMCMC = nmcmc+burnin
  N = npart
  set.seed(mcmseed)
  
  pr <- progress_text()         # displays progress (from plyr)
  pr$init(numMCMC) 
  
  T = length(y)
  X = matrix(0,numMCMC,T)
  q = rep(0,numMCMC)             # q = state variance ## W(t) ~ iid N(0,q)
  phi = rep(0,numMCMC)           # phi-values
  mu = rep(0,numMCMC)            # leverage term
  
  # Initialize the parameters
  phi[1] = init[1]
  q[1] = init[2] 
  mu[1] = init[3]

  
  # hyperparameters
  mu_phi = hyper[1]
  mu_q = hyper[2]
  sigma_phi = hyper[3]
  sigma_q = hyper[4]
  rho = hyper[5]  # negative correlation
  mu_MH = c(0,0)  # tuning parameter
  sigma_MH = sigma_MH # tuning parameter
  
  # Initialize the state by running a PF
  if (parms_2){
    
    u = cpf_as_sv(y, phi[1], q[1], 0, N,  X[1,])   #  changed from X to X[1,]
    particles = u$x           # returned particles
    w = u$w                   # returned weights
    # Draw J
    J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
    X[1,] = particles[J,]      
    
  } else {
    
    u = cpf_as_sv(y, phi[1], q[1], mu[1], N,  X[1,])   #  changed from X to X[1,]
    particles = u$x           # returned particles
    w = u$w                   # returned weights
    # Draw J
    J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
    X[1,] = particles[J,]  
    
  }
  
  
  # Run MCMC loop
  if (parms_2){
    
    for(k in 2:numMCMC){
      # Sample the parameters (phi, sqrt_q) ~ bivariate normal Random Walk Metropolis Hastings
      parms = cbind(phi[k-1], sqrt(q[k-1]))
      (parms_star = parms + mvrnorm(1, mu_MH, sigma_MH))
      
      while(parms_star[2]^2 < 1e-2|parms_star[1]>2){
        (parms_star = parms + mvrnorm(1, mu_MH, sigma_MH))
      }
      
      (g_star = log_g_func(parms_star, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], 0))
      (g_old = log_g_func(parms, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], 0))
      (g_diff = g_star - g_old)
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
      
      # Run CPF-AS
      u1 = cpf_as_sv(y, phi[k], q[k], 0, N, X[k-1,])
      particles = u1$x           # returned particles
      w = u1$w                   # returned weight
      # Draw J (extract a particle trajectory)
      J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
      X[k,] = particles[J,]      # center X
      
      pr$step()
    }
    
  } else {
    
    for(k in 2:numMCMC){
      # Sample the parameters (phi, sqrt_q) ~ bivariate normal, Random Walk Metropolis Hastings
      parms = cbind(phi[k-1], sqrt(q[k-1]))
      (parms_star = parms + mvrnorm(1, mu_MH, sigma_MH))
      
      while(parms_star[2]^2 < 1e-2|parms_star[1]>2){
        (parms_star = parms + mvrnorm(1, mu_MH, sigma_MH))
      }
      
      (g_star = log_g_func(parms_star, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], mu[k-1]))
      (g_old = log_g_func(parms, mu_phi, sigma_phi, mu_q, sigma_q, rho, X[(k-1),], mu[k-1]))
      (g_diff = g_star - g_old)
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
      
      # Sample mu
      sigmu = q[k]/((T-1)*(1-phi[k])^2+(1-phi[k]^2))
      muhat = sigmu*((1-phi[k]^2)*X[k-1,1]/q[k]+(1-phi[k])*(sum(X[k-1,2:T])-phi[k]*sum(X[k-1,1:(T-1)]))/q[k])
      mu[k] = muhat + sqrt(abs(sigmu))*rnorm(1)
      
      while(abs(mu[k]) > 2){
        mu[k] = muhat + sqrt(abs(sigmu))*rnorm(1)
      }
      
      
      # Run CPF-AS
      u = cpf_as_sv(y, phi[k], q[k], mu[k], N, X[k-1,])
      particles = u$x           # returned particles
      w = u$w                   # returned weight
      # Draw J (extract a particle trajectory)
      J = which( (runif(1)-cumsum(w[,T])) < 0 )[1]
      X[k,] = particles[J,]     # center X
      pr$step()
    }#end
    
  }
  
  bi = 1:burnin 
  time2run = proc.time() - ptm
  acp = dim(table(phi[-bi]))/nmcmc*100
  list(phi=phi[-bi], q=q[-bi], mu=mu[-bi], X=X[-bi,], time=time2run, acp=acp)
}#end


#--------------------------------------------------------------------------
cpf_as_sv = function(y, phi, q, mu, N, X){
  # Conditional particle filter with ancestor sampling
  # Input:
  #   y - measurements
  #   phi - transition parameter
  #   q - state noise variance
  #   mu - leverage term
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
      ind = ind[sample.int(N)]; # default: no replacement
      t1 = t-1
      xpred = phi*x[, t1] + (1-phi)*mu # length N Vector
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
    var_now = exp(x[,t]+mu)
    logweights = -y[t]^2/(2*var_now) - 1/2*log(var_now)  # (up to an additive constant)
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