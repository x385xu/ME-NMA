## -------------------------------------------------------------------------
## Case study figures and Q-decomposition tables
##
## This script generates the network plots, forest plots, and study-level
## Q_het decomposition tables for case studies 1--4. The helper function
## make_case_study_outputs() returns plotting functions rather than saving
## figures internally, so that figure dimensions can be adjusted separately
## for journal submission.
## -------------------------------------------------------------------------

library(netmeta)
library(nmadb)
library(meta)
library(dplyr)
library(metafor)
source("make_case_study_outputs.R")

## Case study 1: Topical NSAIDs Versus Placebo for Pain Relief
case1 <- make_case_study_outputs(
  recid = 501330,
  treatment = c("Placebo", "Ketoprofen",
                "Ibuprofen", "Felbinac",
                "Piroxicam", "Indomethacin",
                "Other NSAID"))
options(tibble.print_max = Inf)
case1$df_table
# Figure 2
png("netgraph1.png", width = 1200, height = 600, res = 150)
case1$netgraph()
dev.off()
# Figure 3
png("forest1.png", width = 3200, height = 3500, res = 300)
case1$forest()
dev.off()


## Case study 2: Interventions to Increase Household Possession of a Functioning Smoke Alarm
case2 <- make_case_study_outputs(
  recid = 501235,
  treatment = c("Usual Care", "Education",
                "Education+LCFE","Education+LCFE+HSI",
                "Education+LCFE+Fitting","Education+HSI",
                "Education+LCFE+Fitting+HSI"))

case2$df_table
# Figure 4
png("netgraph2.png", width = 1200, height = 600, res = 150)
case2$netgraph()
dev.off()
# Figure 5
png("forest2.png", width = 3200, height = 2500, res = 300)
case2$forest()
dev.off()

## Case study 3: Biological Therapies Versus Placebo or Standard Care for ACR70 Improvement
case3 <- make_case_study_outputs(
  recid = 473552,
  treatment = c("Placebo/Standard Care", "Abatacept",
                "Adalimumab", "Certolizumab",
                "Etanercept", "Golimumab",
                "Infliximab", "Rituximab", "Tocilizumab"))

case3$df_table
# Figure 6
png("netgraph3.png", width = 1200, height = 600, res = 150)
case3$netgraph()
dev.off()
# Figure 7
png("forest3.png", width = 3200, height = 3500, res = 300)
case3$forest()
dev.off()


## Case study 4: Novel Oral Anticoagulants for Atrial Fibrillation
case4 <- make_case_study_outputs(
  recid = 501212,
  treatment = c("Control", "Apixaban",
                "Dabigatran 110 mg", "Dabigatran 150 mg",
                "Rivaroxaban"))

case4$df_table
# Figure 8
png("netgraph4.png", width = 1200, height = 600, res = 150)
case4$netgraph()
dev.off()
# Figure 9
png("forest5.png", width = 3200, height = 2000, res = 300)
case4$forest()
dev.off()