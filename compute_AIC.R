library(netmeta)
library(nmadb)
library(dplyr)

#  Compute AIC of datasets in nmadb given indices using
#   multiplicative and random/fixed effect models
compute_AIC <- function(dat = dat_nmadb, ind) {
  len <- length(ind)
  AIC_add <- numeric(len)
  AIC_fixed <- numeric(len)
  AIC_mul <- numeric(len)
  taus <- numeric(len)
  phis <- numeric(len)
  Q <- numeric(len)
  Q_pval <- numeric(len)
  
  for (j in seq_along(ind)) {
    i <- ind[j]
    
    net <- runnetmeta(dat_nmadb$recid[i])
    #observed TE
    theta <- net$TE
    m <- net$m
    n <- net$n
    V <- diag(net$seTE^2)
    
    # Multiplicative effect
    theta_me <- net$TE.nma.fixed # fitted TE
    phi <- as.numeric(t(theta-theta_me) %*% 
                        solve(V) %*% (theta-theta_me) / (m-n+1))
    if (phi < 1) {
      phi <- 1
    }
    
    logL_me <- -0.5*(m*log(2*pi)+
                        log(det(phi*V))+ 
                        t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
    AIC_mul[j] <- 2*n - 2*logL_me
    
    # Additive effects
    tau_hat <- net$tau
    sig <- V+tau_hat^2*diag(m)
    theta_ae <- net$TE.nma.random
    logL_ae <- -0.5*(m*log(2*pi)+
                       log(det(sig))+ 
                       t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
    AIC_add[j] <- 2*n-2*logL_ae
    
    # Fixed effect
    sig <- V
    theta_fe <- net$TE.nma.fixed
    logL_fe <- -0.5*(m*log(2*pi)+
                       log(det(V))+ 
                       t(theta-theta_fe) %*% solve(V) %*% (theta-theta_fe))
    AIC_fixed[j] <- 2*n-2*logL_fe
    
    taus[j] <- tau_hat
    phis[j] <- phi
    Q[j] <- net$Q.heterogeneity
    Q_pval[j] <- net$pval.Q.heterogeneity
    
  }
  return(list(
    ind = ind,
    recid =  dat$recid[ind],
    tau = taus,
    phi = phis,
    add  = AIC_add,
    fixed = AIC_fixed,
    mul  = AIC_mul,
    mul_add = AIC_mul - AIC_add,
    mul_fixed = AIC_mul - AIC_fixed,
    Q = Q,
    Q_pval = Q_pval
  ))
}


