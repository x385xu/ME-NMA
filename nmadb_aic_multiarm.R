library(netmeta)
library(nmadb)
library(dplyr)
library(ggplot2)
library(gflnma)
library(ggpattern)

source("get_index_nmadb.R")
source("compute_AIC_multiarm.R")
#source("compute_AIC.R")

#==========load nmadb================
# dat_nmadb <- getNMADB()
# 
# dat_nmadb <- dat_nmadb %>% 
#   select(Record.ID, Title, First.Author, Year, 
#          Number.of.Studies.., Number.of.Treatments, 
#          Type.of.Outcome., Effect.Measure, Fixed.effect.or.Random.effects.model,
#          Primary.Outcome, Description.of.the.outcome.s., 
#          Harmful.Beneficial.Outcome, dataset) %>%
#   rename(recid = Record.ID) %>%
#   mutate(Year = as.numeric(format(as.Date(Year, format="%Y-%m-%d"),"%Y")), 
#          .keep = "unused") 
# 
# saveRDS(dat_nmadb, "dat_nmadb.rds")

dat_nmadb <- readRDS("dat_nmadb.rds")
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

# change method.tau="REML" to get results under the REML method
AIC_rr_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_rr_twoarm, method.tau = "DL")
AIC_or_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_or_twoarm, method.tau = "DL")
AIC_md_twoarm <- compute_AIC_multiarm(dat = dat_nmadb,idx_md_twoarm, method.tau = "DL")

save(AIC_rr_twoarm, AIC_or_twoarm, AIC_md_twoarm,
     file = "AIC_twoarm_results_REML.RData")
load("AIC_twoarm_results_REML.RData")

save(AIC_rr_twoarm, AIC_or_twoarm, AIC_md_twoarm,
     file = "AIC_twoarm_results_DL.RData")
load("AIC_twoarm_results_DL.RData")

df_AICplot_twoarm <- bind_rows(
  tibble(value = AIC_or_twoarm$me_re, measure = "OR", pval = AIC_or_twoarm$Q_pval),
  tibble(value = AIC_rr_twoarm$me_re, measure = "RR", pval = AIC_rr_twoarm$Q_pval),
  tibble(value = AIC_md_twoarm$me_re, measure = "MD", pval = AIC_md_twoarm$Q_pval)
) %>%
  mutate(
    measure = factor(measure, levels = c("OR", "RR", "MD")),
    hetero_group = case_when(
      pval < 0.05 ~ "p1",
      pval >= 0.05 & pval < 0.1 ~ "p2",
      pval >= 0.1 ~ "p3"
    ),
    hetero_group = factor(hetero_group, levels = c("p1", "p2", "p3"))
  )

p <- ggplot(df_AICplot_twoarm, aes(x = value, fill = hetero_group)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    color = "black",
    position = "stack",
    alpha = 0.8
  ) +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  facet_grid(. ~ measure, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(
    values = c(
      "p1" = "grey90",  # light
      "p2" = "grey55",  # medium
      "p3" = "grey10"   # dark
    ),
    labels = c(
      expression(p %in% "[0, 0.05)"),
      expression(p %in% "[0.05, 0.1)"),
      expression(p %in% "[0.1, 1]")
    )
  ) +
  
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x) / 3) * 3,
                             ceiling(max(x) / 3) * 3,
                             by = 3)
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 8),
    minor_breaks = scales::pretty_breaks(n = 16)
  ) +
  labs(y = "", x = "", fill = "") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.35),
    panel.grid.minor.y = element_line(linewidth = 0.2)
  )

p

ggsave("AIC_REML_Signiat10p_stack.png", p, width = 10, height = 4, dpi = 300)


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
