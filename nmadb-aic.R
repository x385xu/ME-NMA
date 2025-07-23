# newest version of netmeta does not have the function pairwise() anymore
#remove.packages("netmeta")
#install.packages("netmeta_2.9-0.tar.gz", repos = NULL, type = "source")
library(netmeta)

# nmadb is removed from cran, can find it in the archive: https://cran.r-project.org/src/contrib/Archive/nmadb/?C=D;O=A
# https://github.com/cran/nmadb
#remove.packages("nmadb")
#install.packages("nmadb_1.2.0.tar.gz", repos = NULL, type = "source")
library(nmadb)
library(dplyr)
library(ggplot2)
# Load NMA database catalog
dat_nmadb <- getNMADB()

# Extract relevant info
dat_nmadb <- dat_nmadb %>% 
  select(Record.ID, Title, First.Author, Year, 
         Number.of.Studies.., Number.of.Treatments, 
         Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
         Primary.Outcome, Description.of.the.outcome.s., 
         Harmful.Beneficial.Outcome, dataset) %>%
  rename(recid = Record.ID) %>%
  mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), 
         .keep = "unused") 

#nma_test <- runnetmeta(dat_nmadb$recid[2])

source("get_index_nmadb.R")
source("compute_AIC.R")

#===============================================================================
#AIC_mul_add = AIC_mul - AIC_add,
#AIC_mul_fixed = AIC_mul - AIC_fixed
#\delta AIC=AIC_add-AIC_fixed
#\delta  AIC < -3 -----additive-1
# -3 <\delta  AIC < 3 ---fixed-2
# \delta AIC > 3 --------fixed-3


#---------or-------------------------------------------------------------------
ind_or <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "odds ratio",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_or <- compute_AIC(dat = dat_nmadb,ind_or)

# diffs_or <- AIC_or$add - AIC_or$fixed
# AICdiff_or_1 <- AIC_or$mul_add[intersect(which(AIC_or$Q_pval < 0.05), 
#                                          which(diffs_or < -3))]
# AICdiff_or_2 <- AIC_or$mul_fixed[intersect(which(AIC_or$Q_pval < 0.05), 
#                                            which(abs(diffs_or) <= 3))]
# AICdiff_or_3 <- AIC_or$mul_fixed[intersect(which(AIC_or$Q_pval < 0.05), 
#                                            which(diffs_or > 3))]

#-------------rr---------------------------------------------------------------
ind_rr <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "risk ratio",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_rr <- compute_AIC(dat = dat_nmadb,ind_rr)

# diffs_rr <- AIC_rr$add - AIC_rr$fixed
# AICdiff_rr_1 <- AIC_rr$mul_add[intersect(which(diffs_rr < -3), 
#                                          which(AIC_rr$Q_pval < 0.05))]
# AICdiff_rr_2 <- AIC_rr$mul_fixed[intersect(which(abs(diffs_rr) <= 3), 
#                                            which(AIC_rr$Q_pval < 0.05))]
# AICdiff_rr_3 <- AIC_rr$mul_fixed[intersect(which(diffs_rr > 3), 
#                                            which(AIC_rr$Q_pval < 0.05))]

#------------md----------------------------------------------------------------
ind_md <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "mean difference",
                           multiarm = FALSE,
                           fixed = TRUE)
AIC_md <- compute_AIC(dat = dat_nmadb,ind_md)

# diffs_md <- AIC_md$add - AIC_md$fixed
# AICdiff_md_1 <- AIC_md$mul_add[intersect(which(diffs_md < -3),
#                                          which(AIC_md$Q_pval < 0.05))]
# AICdiff_md_2 <- AIC_md$mul_fixed[intersect(which(abs(diffs_md) <= 3),
#                                                  which(AIC_md$Q_pval < 0.05))]
# AICdiff_md_3 <- AIC_md$mul_fixed[intersect(which(diffs_md > 3),
#                                            which(AIC_md$Q_pval < 0.05))]

#===========AIC plot===============================================================
AIC_diff_or <- AIC_or$mul_add[which(AIC_or$Q_pval < 0.05)]
AIC_diff_rr <- AIC_rr$mul_add[which(AIC_rr$Q_pval < 0.05)]
AIC_diff_md <- AIC_md$mul_add[which(AIC_md$Q_pval < 0.05)]

df_AICplot <- bind_rows(
  tibble(value = AIC_diff_or, measure = "OR"),
  tibble(value = AIC_diff_rr, measure = "RR"),
  tibble(value = AIC_diff_md, measure = "MD")) %>%
  mutate(measure = factor(measure, levels = c("OR", "RR", "MD")))

ggplot(df_AICplot, aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(.~measure, scales = "free") +
  labs(y = "Count", x="") +
  theme_minimal(base_size = 12)

median(AIC_diff_or)
median(AIC_diff_rr)
median(AIC_diff_md)

length(which(AIC_diff_or < 0))
length(which(AIC_diff_rr < 0))
length(which(AIC_diff_md < 0))

length(AIC_diff_or)
length(AIC_diff_rr)
length(AIC_diff_md)

# # build a data-frame with one row per observation
# df <- bind_rows(
#   tibble(value = AICdiff_or_1, category = "Additive better", measure = "OR"),
#   tibble(value = AICdiff_or_2, category = "Similar",         measure = "OR"),
#   tibble(value = AICdiff_or_3, category = "Fixed better",    measure = "OR"),
#   tibble(value = AICdiff_rr_1, category = "Additive better", measure = "RR"),
#   tibble(value = AICdiff_rr_2, category = "Similar",         measure = "RR"),
#   tibble(value = AICdiff_rr_3, category = "Fixed better",    measure = "RR"),
#   tibble(value = AICdiff_md_1, category = "Additive better", measure = "MD"),
#   tibble(value = AICdiff_md_2, category = "Similar",         measure = "MD"),
#   tibble(value = AICdiff_md_3, category = "Fixed better",    measure = "MD")
# ) %>%
#   mutate(
#     # specify the column and row order
#     category = factor(category,
#                       levels = c("Additive better",
#                                  "Similar",
#                                  "Fixed better")),
#     measure  = factor(measure,  levels = c("OR", "RR", "MD")))
# 
# 
# ggplot(df, aes(x = value)) + 
#   geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
#   geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
#   facet_grid(measure ~ category, scales = "free") +
#   labs(x = expression(Delta~"AIC (Multiplicative - Additive or Fixed)"), 
#        y = "Count") +
#   theme_minimal(base_size = 12)

#======================================================
# idx_or <- intersect(which(AIC_or$Q_pval < 0.05), which(diffs_or > -3))
# 
# table_or <- data.frame(
#   Q_pval = AIC_or$Q_pval[idx_or],
#   Q_heterogeneity = AIC_or$Q[idx_or],
#   tau = AIC_or$tau[idx_or],
#   phi = AIC_or$phi[idx_or],
#   RE_FE = diffs_or[idx_or],
#   ME_FE   = AIC_or$mul_fixed[idx_or]
# )
# 
# print(
#   format(table_or, digits = 2, nsmall = 2),
#   row.names = FALSE
# )
# 
# write.table(table_or, sep = "\t", row.names = FALSE, quote = FALSE)
# 
# 
# 
# idx_rr <- intersect(which(AIC_rr$Q_pval < 0.05), which(diffs_rr > -3))
# 
# table_rr <- data.frame(
#   Q_pval = AIC_rr$Q_pval[idx_rr],
#   Q_heterogeneity = AIC_rr$Q[idx_rr],
#   tau = AIC_rr$tau[idx_rr],
#   phi = AIC_rr$phi[idx_rr],
#   RE_FE = diffs_rr[idx_rr],
#   ME_FE   = AIC_rr$mul_fixed[idx_rr]
# )
# 
# print(
#   format(table_rr, digits = 2, nsmall = 2),
#   row.names = FALSE
# )
# 
# 
# 
# 
# idx_md <- intersect(which(AIC_md$Q_pval < 0.05), which(diffs_md > -3))
# 
# table_md <- data.frame(
#   Q_pval = AIC_md$Q_pval[idx_md],
#   Q_heterogeneity = AIC_md$Q[idx_md],
#   tau = AIC_md$tau[idx_md],
#   phi = AIC_md$phi[idx_md],
#   RE_FE = diffs_md[idx_md],
#   ME_FE   = AIC_md$mul_fixed[idx_md]
# )
# 
# print(
#   format(table_md, digits = 2, nsmall = 2),
#   row.names = FALSE
# )