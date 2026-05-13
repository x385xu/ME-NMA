# Example code: fitting the additive random-effects (RE) model and
# multiplicative heterogeneity (ME) model for one NMA dataset

library(meta)
library(netmeta)

# Load one case-study dataset
twoarm_data_list <- readRDS("nmadb_twoarm_data_all.rds")
recid <-  as.character(501212)
dat <- twoarm_data_list[[recid]]

dat$type 

# The outcome is binary, so the arm-based data must first be converted
# to contrast-based data using pairwise().

D <- dat$data

pw <- meta::pairwise(
  treat = D$t,
  event = D$r,
  n = D$n,
  studlab = D$id,
  data = D,
  sm = dat$effect,
  allstudies = TRUE
)

# Replace numeric treatment codes by descriptive treatment names.
treatment_labels <- c(
  "1" = "Control",
  "2" = "Apixaban",
  "3" = "Dabigatran 110 mg",
  "4" = "Dabigatran 150 mg",
  "5" = "Rivaroxaban"
)
pw$treat1 <- treatment_labels[as.character(pw$treat1)]
pw$treat2 <- treatment_labels[as.character(pw$treat2)]

# -------------------------------------------------------------------------
# 1. Additive random-effects model
# -------------------------------------------------------------------------

fit <- netmeta::netmeta(
  TE = pw$TE,
  seTE = pw$seTE,
  treat1 = pw$treat1,
  treat2 = pw$treat2,
  studlab = pw$studlab,
  data = pw,
  sm = dat$effect,
  random = TRUE,
  method.tau = "DL" # change to "REML" for REML estimation
)

# Estimated additive heterogeneity variance
tau_hat <- fit$tau


Q_total <- fit$Q
df <- fit$df.Q
# Estimate multiplicative heterogeneity
phi_hat <- max(1, Q_total / df)
# The ME model uses the common-effect point estimates, but inflates the
# within-study variances by a multiplicative factor phi.
TE_ME <- fit$TE.common
seTE_ME <- sqrt(phi_hat) * fit$seTE.common


