# install.packages("netmeta")
library(netmeta)
# install.packages("dplyr") 
library(dplyr)
# install.packages("tibble") 
library(tibble)
# install.packages("ggplot2") 
library(ggplot2)
# install.packages("ggpattern") 
library(ggpattern)
# install.packages("patchwork") 
library(patchwork)
# install.packages("xtable") 
library(xtable)

source("functions/fit_netmeta.R")
source("functions/compute_AIC.R")
source("functions/make_AIC_histogram.R")
source("functions/make_quantile_table.R")

# Load data
load("data/idx_nmadb.RData")
twoarm_data_list <- readRDS("data/nmadb_twoarm_data_all.rds")

# Create results folder if it does not already exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Check the total number of available NMA datasets. (266)
nrow(idx_nmadb)

# DL estimation ----------------------------------------------------------------

# Compute AIC and related statistics for two-arm NMAs reporting OR, RR, or MD.
# The between-study heterogeneity parameter tau is estimated using the
# DerSimonian--Laird (DL) method.
# NOTE: This step may take several minutes to run
AIC_DL <- compute_AIC(twoarm_data_list,
                      idx_nmadb,
                      measure = c("odds ratio", "risk ratio", "mean difference"),
                      twoarm_only = TRUE,
                      method.tau = "DL")

# Figure 1: Plot the distribution of Delta AIC values under DL estimation.
p_DL <- make_AIC_histogram(AIC_DL)
ggsave(
  filename = file.path("results", "AIC_DL_stack_interval.png"),
  plot = p_DL,
  width = 20,
  height = 8,
  dpi = 200
)

# Appendix A.3 Table 1: DL quantile summaries ----------------------------------
df_AIC_DL <- tibble::as_tibble(AIC_DL)
# Create quantile summaries under different inclusion criteria.
quantile_tables_DL <- list(
  all = df_AIC_DL %>%
    make_quantile_table(),
  
  nonzero_AICs = df_AIC_DL %>%
    filter(me_re != 0) %>%
    make_quantile_table(),
  
  heterogeneous = df_AIC_DL %>%
    filter(Q_pval < 0.1) %>%
    make_quantile_table(),
  
  non_heterogeneous = df_AIC_DL %>%
    filter(Q_pval >= 0.1 | is.na(Q_pval)) %>%
    make_quantile_table()
)
# Print full tibble output in the console.
options(tibble.width = Inf)
quantile_tables_DL
# Combine the list of summary tables into one table.
# The .id column records which inclusion criterion was used.
quantile_tables_DL <- dplyr::bind_rows(
  quantile_tables_DL,
  .id = "Inclusion"
) %>%
  dplyr::relocate(Inclsion, .after = effect_measure)
# Generate LaTeX code for Appendix Table 1.
xtable::print.xtable(
  xtable::xtable(quantile_tables_DL),
  include.rownames = FALSE
)

# REML estimation --------------------------------------------------------------

# Compute AIC and related statistics for the same set of two-arm NMAs,
# now estimating tau using REML.
# REML estimation failed to converge for network recid = 479770.
# Since the DL method estimated tau = 0, tau is recorded
# as 0 for this network.
# NOTE: This step may take several minutes to run
AIC_REML <- compute_AIC(twoarm_data_list,
                        idx_nmadb,
                        measure = c("odds ratio", "risk ratio", "mean difference"),
                        twoarm_only = TRUE,
                        method.tau = "REML")
# Appendix A.4 Figure 10: Plot the distribution of Delta AIC values under REML estimation.
p_REML <- make_AIC_histogram(AIC_REML)
ggsave(
  filename = file.path("results", "AIC_REML_stack_interval.png"),
  plot = p_DL,
  width = 20,
  height = 8,
  dpi = 200
)

# Appendix A.4 Table 2: REML quantile summaries --------------------------------
df_AIC_REML <- tibble::as_tibble(AIC_REML)
quantile_tables_REML <- list(
  all = df_AIC_REML %>%
    make_quantile_table(),
  
  without_zero_AICs = df_AIC_REML %>%
    filter(me_re != 0) %>%
    make_quantile_table(),
  
  heterogeneous_only = df_AIC_REML %>%
    filter(Q_pval < 0.1) %>%
    make_quantile_table(),
  
  non_heterogeneous_only = df_AIC_REML %>%
    filter(Q_pval >= 0.1 | is.na(Q_pval)) %>%
    make_quantile_table()
)
# Print full tibble output in the console.
options(tibble.width = Inf)
quantile_tables_REML
# Combine the REML summary tables into one table.
quantile_tables_REML <- dplyr::bind_rows(
  quantile_tables_REML,
  .id = "Inclusion"
) %>%
  dplyr::relocate(Inclsion, .after = effect_measure)
# Generate LaTeX code for Table 2.
xtable::print.xtable(
  xtable::xtable(quantile_tables_REML),
  include.rownames = FALSE
)
