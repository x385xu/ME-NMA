## -------------------------------------------------------------------------
## This dataset required manual correction of the standard errors before
## fitting the network meta-analysis model. 
## See Appendix A.4 Rheumatoid Arthritis Treatments Versus Placebo 
## for Pain Improvement
## -------------------------------------------------------------------------
# install.packages("netmeta") 
library(netmeta)
source("functions/fit_netmeta.R")

## Load network
twoarm_data_list <- readRDS("data/nmadb_twoarm_data_all.rds")
recid <-  as.character(480039)
dat_raw <- twoarm_data_list[[recid]]
net_raw <- fit_netmeta(
  indata = dat_raw,
  model = "random"
)

## Recompute standard errors for each treatment contrast
sd_arm <- dat_raw$data[, "sd"]

idx_1 <- seq(1, length(sd_arm) - 1, by = 2)
idx_2 <- seq(2, length(sd_arm),     by = 2)

se_corrected <- sqrt(sd_arm[idx_1]^2 + sd_arm[idx_2]^2)

## Construct corrected contrast-level dataset
dat_net <- data.frame(
  treat1 = net_raw$treat1,
  treat2 = net_raw$treat2,
  TE     = net_raw$TE,
  seTE   = se_corrected
)

## Fit network meta-analysis using corrected standard errors
net_corr <- netmeta(
  TE     = TE,
  seTE   = seTE,
  treat1 = treat1,
  treat2 = treat2,
  data   = dat_net
)

## Extract quantities used in AIC calculations
theta <- net_corr$TE
m     <- net_corr$m
n     <- net_corr$n
V     <- diag(net_corr$seTE^2)


## Multiplicative-effects model

theta_me <- net_corr$TE.nma.fixed

phi_hat <- as.numeric(net_corr$Q / (m - n + 1))

Sigma_me <- phi_hat * V

loglik_me <- -0.5 * (
  m * log(2 * pi) +
    log(det(Sigma_me)) +
    t(theta - theta_me) %*% solve(Sigma_me) %*% (theta - theta_me)
)

AIC_me <- 2 * n - 2 * as.numeric(loglik_me)


## Additive random-effects model

theta_re <- net_corr$TE.nma.random
tau_hat  <- net_corr$tau

Sigma_re <- V + tau_hat^2 * diag(m)

loglik_re <- -0.5 * (
  m * log(2 * pi) +
    log(det(Sigma_re)) +
    t(theta - theta_re) %*% solve(Sigma_re) %*% (theta - theta_re)
)

AIC_re <- 2 * n - 2 * as.numeric(loglik_re)

## Difference in AIC
## Negative values favour the multiplicative-effects model;
## positive values favour the additive random-effects model.
delta_AIC <- AIC_me - AIC_re
delta_AIC
