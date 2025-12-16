library(netmeta)
library(nmadb)
library(dplyr)
library(ggplot2)
# install.packages("remotes")
# remotes::install_github("caitlin-h-daly/gflnma")
library(gflnma)

# Compute AIC of datasets in nmadb given indices using
#   multiplicative and random/fixed effect models

# Build reference-invariant B = bdiag( B_i ) with B_i = C_i P_i C_i^T
# dat must have columns: treat1, treat2, studlab (study id)
build_B_from_pairs <- function(dat) {
  # split row indices by study, preserving row order
  idx_by_study <- split(seq_len(nrow(dat)), factor(dat$studlab, levels = unique(dat$studlab)))
  
  # for each study, make its contrast matrix C_i using study-local arm labels
  B_blocks <- lapply(idx_by_study, function(idx) {
    di <- dat[idx, , drop = FALSE]
    arms <- sort(unique(c(di$treat1, di$treat2)))        # study-local arm labels
    T_i  <- length(arms)
    k_i  <- nrow(di)                                     # number of contrasts in this study
    
    # map global arm labels to study-local column indices 1..T_i
    arm_pos <- match
    
    # build C_i: each row is e_a - e_b with respect to the study-local arm indices
    Ci <- matrix(0, nrow = k_i, ncol = T_i)
    for (r in seq_len(k_i)) {
      a <- arm_pos(di$treat1[r], arms)
      b <- arm_pos(di$treat2[r], arms)
      Ci[r, a] <-  1
      Ci[r, b] <- -1
    }
    
    # centering projector P_i
    Pi <- diag(T_i) - matrix(1, T_i, T_i) / T_i
    
    # B_i = C_i P_i C_i^T
    Bi <- Ci %*% Pi %*% t(Ci)
    Bi
  })
  
  # block-diagonal combine (your bdiag helper is fine)
  bdiag <- function(mats) {
    nr <- sum(sapply(mats, nrow)); nc <- sum(sapply(mats, ncol))
    out <- matrix(0, nr, nc); r0 <- 0; c0 <- 0
    for (A in mats) {
      idxr <- (r0 + 1):(r0 + nrow(A)); idxc <- (c0 + 1):(c0 + ncol(A))
      out[idxr, idxc] <- A; r0 <- r0 + nrow(A); c0 <- c0 + ncol(A)
    }
    out
  }
  
  B <- bdiag(B_blocks)
  return(B)
}

compute_AIC_multiarm <- function(dat = dat_nmadb, ind) {
  len <- length(ind)
  AIC_re <- numeric(len)
  AIC_fe <- numeric(len)
  AIC_me <- numeric(len)
  taus <- numeric(len)
  phis <- numeric(len)
  Q <- numeric(len)
  Q_pval <- numeric(len)
  
  #computes log(det(A)) using cholesky decomposition
  logdet <- function(A) { R <- chol(A); 2*sum(log(diag(R))) }
  
  for (j in seq_along(ind)) {
    i <- ind[j]
    net <- tryCatch({runnetmeta(dat_nmadb$recid[i])}, 
                    error = function(e) {return(NULL)})
    # Skip if net is not a list
    if (!is.list(net)) next

    # ----------Multiplicative effect model------------------
    dat_raw <- data.frame(
      TE      = net$TE,
      seTE    = net$seTE,
      treat1  = net$treat1,
      treat2  = net$treat2,
      studlab = net$studlab
    )
    
    dat_list <- prep_gfl_data(dat = dat_raw, 
                              TE = "TE",
                              seTE = "seTE", 
                              treatment1 = "treat1", 
                              treatment2 = "treat2",
                              studlab = "studlab",
                              modtype = "ME")
    
    V <- dat_list$var_cov
    y <- dat_list$dat$TE
    ref <- net$reference.group
    X <- get_design_matrix(treatment1 = dat_list$dat$treat1,
                           treatment2 = dat_list$dat$treat2,
                           ref = ref)
    r <- ncol(X)
    M <- length(y)
    
    U  <- chol(solve(V))   # upper-triangular, U'U = V^{-1}
    Xw <- U %*% X
    yw <- U %*% y
    # equivalent to gls with weight matrix V^-1
    bhat <- coef(lm(yw ~ -1 + Xw))
    
    res   <- y - as.vector(X %*% bhat)
    SSE   <- as.numeric(crossprod(res, solve(V, res)))
    phi <- SSE / (M-r)
    phi <- max(1,phi)
    
    logL_me <- -0.5*(M*log(2*pi)+M*log(phi) + logdet(V) + SSE/phi)
    AIC_me[j] <- 2*(r+1) - 2*logL_me
    
    #--------------fixed-effect--------------------------------------------
    logL_fe <- -0.5*(M*log(2*pi) + logdet(V) + SSE)
    AIC_fe[j] <- 2*(r+1) - 2*logL_fe

    #------------Additive random-effects model------------------------------
    net_data <- as.data.frame(cbind(net$treat1, net$treat2))
    colnames(net_data) <- c("treat1", "treat2")
    net_data$TE <- net$TE
    net_data$seTE <- net$seTE
    
    net_2 <- netmeta(TE, seTE, treat1, treat2, data = net_data, method.tau = "REML")
    
    tau2 <- net_2$tau2
    # Using your dat_list (same ordering as y and V):
    B <- build_B_from_pairs(dat_list$dat)
    
    V_RE <- V + tau2 * B * 0.5
    R  <- chol(solve(V_RE))
    Xw <- R %*% X
    yw <- R %*% y
    bhat_re <- coef(lm(yw ~ -1 + Xw))
    res_re  <- y - as.vector(X %*% bhat_re)
    SSER    <- as.numeric(crossprod(res_re, solve(V_RE, res_re)))
    logL_re <- -0.5*(M*log(2*pi) + logdet(V_RE) + SSER)
    AIC_re[j] <- 2*(r+1) - 2*logL_re

    
    taus[j] <- net$tau
    phis[j] <- phi
    Q[j] <- net$Q.heterogeneity
    Q_pval[j] <- net$pval.Q.heterogeneity
    
  }
  return(list(
    ind = ind,
    recid =  dat$recid[ind],
    tau = taus,
    phi = phis,
    AIC_re  = AIC_re,
    AIC_fe = AIC_fe,
    AIC_me  = AIC_me,
    me_re = AIC_me - AIC_re,
    me_fe = AIC_me - AIC_fe,
    Q = Q,
    Q_pval = Q_pval
  ))
}


