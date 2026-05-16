# Example code: fitting the additive random-effects (RE) model and
# multiplicative heterogeneity (ME) model for one NMA dataset

library(meta)
library(netmeta)

# Load case study 2 dataset
twoarm_data_list <- readRDS(file.path("data/nmadb_twoarm_data_all.rds"))
recid <-  as.character(501235)
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
  "1" = "Usual Care",
  "2" = "Education",
  "3" = "Education+LCFE",
  "4" = "Education+LCFE+HSI",
  "5" = "Education+LCFE+Fitting",
  "6" = "Education+HSI",
  "7" = "Education+LCFE+Fitting+HSI"
)

pw$treat1 <- treatment_labels[as.character(pw$treat1)]
pw$treat2 <- treatment_labels[as.character(pw$treat2)]

# Fit the model with netmeta()
fit <- netmeta::netmeta(
  TE = pw$TE,
  seTE = pw$seTE,
  treat1 = pw$treat1,
  treat2 = pw$treat2,
  studlab = pw$studlab,
  data = pw,
  sm = dat$effect,
  random = TRUE,    # Set random = TRUE to fit RE model. The output still contains CE(ME) estimates
  method.tau = "DL" # change to "REML" for REML estimation
)

# Create netgraph
# Note: This network graph may differ visually from the version shown in the
# original paper because the treatment ordering and node layout are determined
# automatically by the plotting algorithm. The underlying network structure and
# study connections remain unchanged.
netgraph(
  fit,
  plastic = FALSE,
  iterate = FALSE,
  number.of.studies = TRUE,
  col = "grey",
  cex = 1,
  points = TRUE,
  col.points = "orange",
  cex.points = 3,
  cex.number.of.studies = 1,
  pos.number.of.studies = 0.4,
  lwd = 2,
  offset = 0.06
)

# -------------------------------------------------------------------------
# Additive random-effects (RE) model
# -------------------------------------------------------------------------

tau_hat <- fit$tau

# RE model estimates
TE_RE <- fit$TE.random
seTE_RE <- fit$seTE.random

# -------------------------------------------------------------------------
# Multiplicative heterogeneity (ME) model
# -------------------------------------------------------------------------

Q_total <- fit$Q
df <- fit$df.Q

# Estimate multiplicative heterogeneity parameter
phi_hat <- max(1, Q_total / df)

# The ME model uses the common-effect point estimates, while inflating
# the within-study variances by the multiplicative factor phi.
TE_ME <- fit$TE.common
seTE_ME <- sqrt(phi_hat) * fit$seTE.common






