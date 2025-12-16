library(netmeta)
library(nmadb)
library(dplyr)
library(ggplot2)
library(gflnma)

source("get_index_nmadb.R")
source("compute_AIC_multiarm.R")
#source("compute_AIC.R")

#==========load nmadb================
dat_nmadb <- getNMADB()

dat_nmadb <- dat_nmadb %>% 
  select(Record.ID, Title, First.Author, Year, 
         Number.of.Studies.., Number.of.Treatments, 
         Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
         Primary.Outcome, Description.of.the.outcome.s., 
         Harmful.Beneficial.Outcome, dataset) %>%
  rename(recid = Record.ID) %>%
  mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), 
         .keep = "unused") 

#==========get indices of multiarm studies============
idx_or <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "odds ratio")

idx_rr <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "risk ratio")

idx_md <-  get_index_nmadb(dat = dat_nmadb,
                           measure = "mean difference")

#=================compute AIC of two-arm studies===================
idx_or_twoarm <- na.omit(idx_or$twoarm)
idx_rr_twoarm <- na.omit(idx_rr$twoarm)
idx_md_twoarm <- na.omit(idx_md$twoarm)
AIC_rr_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_rr_twoarm)
AIC_or_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_or_twoarm)
AIC_md_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_md_twoarm)

AIC_diff_or_twoarm <- AIC_or_twoarm$me_re[which(AIC_or_twoarm$Q_pval < 0.05)]
AIC_diff_rr_twoarm <- AIC_rr_twoarm$me_re[which(AIC_rr_twoarm$Q_pval < 0.05)]
AIC_diff_md_twoarm <- AIC_md_twoarm$me_re[which(AIC_md_twoarm$Q_pval < 0.05)]

AIC_diff_or_twoarm_DL <- AIC_diff_or_twoarm
AIC_diff_rr_twoarm_DL <- AIC_diff_rr_twoarm
AIC_diff_md_twoarm_DL <- AIC_diff_md_twoarm

df_AICplot_twoarm_DL <- bind_rows(
  tibble(value = AIC_diff_or_twoarm_DL, measure = "OR"),
  tibble(value = AIC_diff_rr_twoarm_DL, measure = "RR"),
  tibble(value = AIC_diff_md_twoarm_DL, measure = "MD")) %>%
  mutate(measure = factor(measure, levels = c("OR", "RR", "MD")))

ggplot(df_AICplot_twoarm_DL, aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(. ~ measure, scales = "free") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1),
    labels = function(x) ifelse(x %% 3 == 0, x, ""),
  ) +
  labs(y = "", x = "") +
  theme_minimal(base_size = 12)

AIC_diff_or_twoarm_REML <- AIC_diff_or_twoarm
AIC_diff_rr_twoarm_REML <- AIC_diff_rr_twoarm
AIC_diff_md_twoarm_REML <- AIC_diff_md_twoarm

df_AICplot_twoarm_REML <- bind_rows(
  tibble(value = AIC_diff_or_twoarm_REML, measure = "OR"),
  tibble(value = AIC_diff_rr_twoarm_REML, measure = "RR"),
  tibble(value = AIC_diff_md_twoarm_REML, measure = "MD")) %>%
  mutate(measure = factor(measure, levels = c("OR", "RR", "MD")))

ggplot(df_AICplot_twoarm_REML, aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(. ~ measure, scales = "free") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1),
    labels = function(x) ifelse(x %% 3 == 0, x, ""),
  ) +
  labs(y = "", x = "") +
  theme_minimal(base_size = 12)
# check the cases
# or and md have the same outliers
which(AIC_diff_rr_twoarm_DL < -3)
which(AIC_diff_rr_twoarm_REML < -3)

median(AIC_diff_rr_twoarm_REML)

#==============multiarm=============================
AIC_or_multiarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_or$multiarm)
AIC_rr_multiarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_rr$multiarm)

AIC_diff_or_multiarm <- AIC_or_multiarm$me_re[which(AIC_or_multiarm$Q_pval < 0.05)]
AIC_diff_rr_multiarm <- AIC_rr_multiarm$me_re[which(AIC_rr_multiarm$Q_pval < 0.05)]
AIC_diff_md_multiarm <- corrected_AIC_md$me_re[which(corrected_AIC_md$Q_pval < 0.05)]

df_AICplot_multiarm <- bind_rows(
  tibble(value = AIC_diff_or_multiarm, measure = "OR"),
  tibble(value = AIC_diff_rr_multiarm, measure = "RR"),
  tibble(value = AIC_diff_md_multiarm, measure = "MD")) %>%
  mutate(measure = factor(measure, levels = c("OR", "RR", "MD")))

ggplot(df_AICplot_multiarm, aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(.~measure, scales = "free") +
  labs(y = "Count", x="") +
  theme_minimal(base_size = 12)

median(AIC_diff_or_multiarm)
median(AIC_diff_rr_multiarm)
median(AIC_diff_md_multiarm)

#============check for p-val >0.05=======
AIC_diff_or_multiarm_2 <- AIC_or_multiarm$me_re[which(AIC_or_multiarm$Q_pval > 0.05)]
AIC_diff_rr_multiarm_2 <- AIC_rr_multiarm$me_re[which(AIC_rr_multiarm$Q_pval > 0.05)]
AIC_diff_md_multiarm_2 <- corrected_AIC_md$me_re[which(corrected_AIC_md$Q_pval > 0.05)]

df_AICplot_multiarm_2 <- bind_rows(
  tibble(value = AIC_diff_or_multiarm_2, measure = "OR"),
  tibble(value = AIC_diff_rr_multiarm_2, measure = "RR"),
  tibble(value = AIC_diff_md_multiarm_2, measure = "MD")) %>%
  mutate(measure = factor(measure, levels = c("OR", "RR", "MD")))

ggplot(df_AICplot_multiarm_2, aes(x = value)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black", fill = "white") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(.~measure, scales = "free") +
  labs(y = "Count", x="") +
  theme_minimal(base_size = 12)
