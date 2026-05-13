## -------------------------------------------------------------------------
## Case study figures and Q-decomposition tables
##
## This script generates the network plots, forest plots, and study-level
## Q_het decomposition tables for case studies 1--4. The helper function
## analyse_case_study() returns plotting functions rather than saving
## figures internally, so that figure dimensions can be adjusted separately
## for journal submission.
## -------------------------------------------------------------------------
# install.packages("netmeta") 
library(netmeta)
# install.packages("meta") 
library(meta)
# install.packages("dplyr") 
library(dplyr)
# install.packages("metafor") 
library(metafor)
source("functions/fit_netmeta.R")
source("functions/analyse_case_study.R")

twoarm_data_list <- readRDS("data/nmadb_twoarm_data_all.rds")

# Create results folder if it does not already exist
if (!dir.exists("results")) {
  dir.create("results")
}

## Case study 1: Topical NSAIDs Versus Placebo for Pain Relief
case1 <- analyse_case_study(
  recid = 501330,
  treatment = c("Placebo", "Ketoprofen",
                "Ibuprofen", "Felbinac",
                "Piroxicam", "Indomethacin",
                "Other NSAID"))

# Figure 2
png(filename = file.path("results", "netgraph1.png"),
    width = 1200, height = 600, res = 150)
case1$netgraph()
dev.off()
# Figure 3
png(filename = file.path("results", "forest1.png"), 
    width = 3200, height = 3500, res = 300)
case1$forest()
dev.off()
# Appendix A.5 Table 3
options(tibble.print_max = Inf)
case1$df_table

## Case study 2: Interventions to Increase Household Possession of a Functioning Smoke Alarm
case2 <- analyse_case_study(
  recid = 501235,
  treatment = c("Usual Care", "Education",
                "Education+LCFE","Education+LCFE+HSI",
                "Education+LCFE+Fitting","Education+HSI",
                "Education+LCFE+Fitting+HSI"))


# Figure 4
png(filename = file.path("results", "netgraph2.png"), 
    width = 1200, height = 600, res = 150)
case2$netgraph()
dev.off()
# Figure 5
png(filename = file.path("results", "forest2.png"), 
    width = 3200, height = 2500, res = 300)
case2$forest()
dev.off()
# Appendix A.5 Table 4
case2$df_table

## Case study 3: Biological Therapies Versus Placebo or Standard Care for ACR70 Improvement
case3 <- analyse_case_study(
  recid = 473552,
  treatment = c("Placebo/Standard Care", "Abatacept",
                "Adalimumab", "Certolizumab",
                "Etanercept", "Golimumab",
                "Infliximab", "Rituximab", "Tocilizumab"))


# Figure 6
png(filename = file.path("results", "netgraph3.png"), 
    width = 1200, height = 600, res = 150)
case3$netgraph()
dev.off()
# Figure 7
png(filename = file.path("results", "forest3.png"), 
    width = 3200, height = 3500, res = 300)
case3$forest()
dev.off()
# Appendix A.5 Table 5
case3$df_table

## Case study 4: Novel Oral Anticoagulants for Atrial Fibrillation
case4 <- analyse_case_study(
  recid = 501212,
  treatment = c("Control", "Apixaban",
                "Dabigatran 110 mg", "Dabigatran 150 mg",
                "Rivaroxaban"))


# Figure 8
png(filename = file.path("results", "netgraph4.png"), 
    width = 1200, height = 600, res = 150)
case4$netgraph()
dev.off()
# Figure 9
png(filename = file.path("results", "forest4.png"), 
    width = 3200, height = 2000, res = 300)
case4$forest()
dev.off()
# Appendix A.5 Table 6
case4$df_table

