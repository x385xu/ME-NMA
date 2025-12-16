library(netmeta)
library(nmadb)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)

# 357, recid=501330
dat$nmadb[357,]
net <- runnetmeta(501330)
read <- readByID(501330)

theta <- net$TE
m <- net$m
n <- net$n
V <- diag(net$seTE^2)

# Multiplicative effect
theta_me <- net$TE.nma.fixed # fitted TE
phi <- as.numeric(t(theta-theta_me) %*% 
                    solve(V) %*% (theta-theta_me) / (m-n+1))
logL_me <- -0.5*(m*log(2*pi)+
                   log(det(phi*V))+ 
                   t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
AIC_mul <- 2*n - 2*logL_me

# Additive effects
tau_hat <- net$tau
sig <- V+tau_hat^2*diag(m)
theta_ae <- net$TE.nma.random
logL_ae <- -0.5*(m*log(2*pi)+
                   log(det(sig))+ 
                   t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
AIC_add <- 2*n-2*logL_ae


phi
tau_hat

AIC_mul-AIC_add

net$Q.inconsistency
net$Q.heterogeneity
net$pval.Q.heterogeneity
net$Q.decomp

#=========netgraph=================
treatment <- c("Placebo", "ketoprofen",
               "ibuprofen", "felbinac",
               "piroxicam", "indomethacin",
               "other NSAID")

netgraph(net,
         labels = treatment,
         plastic = FALSE, 
         iterate=FALSE,
         number.of.studies = TRUE, 
         col="black", 
         cex=1, 
         multiarm=FALSE, 
         points=TRUE, 
         col.points="black", 
         cex.points=3, 
         cex.number.of.studies = 1,
         lwd = 3,
         offset = 0.04,
#         srt.labels = -20,
#         rotate = -45
)

#net$TE = -log risk ratio
# [1] -0.46753754 -1.28093385 -0.65392647 -0.25608771 -1.52696169 -0.60261315 -1.73460106
# [8] -0.06007812 -0.69314718 -1.13943428  0.00000000 -0.46752882 -0.57054486 -1.07967509
# [15] -0.49899117 -0.24783616 -0.07333127 -0.07283732 -0.55961579 -1.24653242 -0.17429641
# [22] -0.10227885 -0.22314355 -0.83290912 -0.17289642 -2.39789527 -1.94591015 -0.34174929
# [29] -0.44731222
#===========forest plot====================

d_mul <- net$TE.fixed[-1,1]
se_mul <- sqrt(phi)*net$seTE.fixed[-1,1]
ci_mul_upper <- d_mul + qnorm(0.975)*se_mul
ci_mul_lower <- d_mul + qnorm(0.025)*se_mul

d_add <- net$TE.random[-1,1]
ci_add_lower <- net$lower.random[-1,1]
ci_add_upper <- net$upper.random[-1,1]

treatment_full <- c("Placebo","ketoprofen","ibuprofen","felbinac",
                    "piroxicam","indomethacin","other NSAID")

# drop "Placebo"
treats <- treatment_full[-1]

# Rank treatments by RE (d_add), descending:
order_idx  <- order(d_add, decreasing = TRUE)
treats_ord <- treats[order_idx]
#d_add[order_idx]
#treats[order_idx]

# df_plot: df of the estimates from RE and ME
df_plot <- tibble(
  treatment = rep(treats, 2),
  model     = rep(c("RE","ME"), each = length(treats_ord)),
  est       = c(d_add, d_mul),
  lo        = c(ci_add_lower, ci_mul_lower),
  hi        = c(ci_add_upper, ci_mul_upper)
  ) %>%
  mutate(
    treatment = factor(treatment, levels = treats_ord),
    y_base    = as.numeric(treatment),
    # change location on the y-axis
    y_num     = y_base + if_else(model=="RE", +0.28, +0.1)
  ) 


obs_df <- tibble(
  trt_code = as.integer(net$treat2),
  est      = -net$TE,
  se       = net$seTE,
  Q_het_decomp = df_02$Q_het_decomp,
) %>%
  filter(!is.na(est), !is.na(se)) %>%
  mutate(
    treatment = treatment_full[trt_code],
    lo        = est + qnorm(0.025)*se,
    hi        = est + qnorm(0.975)*se,
    treatment = factor(treatment, levels = treats_ord),
    y_base    = as.numeric(treatment),
    Q_het_decomp = ifelse(Q_het_decomp >= net$Q.heterogeneity/29, 
                          Q_het_decomp, NA)
  ) %>%
  group_by(treatment) %>%
  arrange(se, .by_group = TRUE) %>%
  mutate(
    n_obs = n(),
    idx   = row_number(),
    # keep them within ±0.2 of the integer row
    y_num = y_base + (idx - (n_obs+1)/2)*(0.4/n_obs)-0.2
  ) %>%
  ungroup() %>%
  mutate(
    alpha_val = rescale(max(se, na.rm=TRUE) - se, to = c(0.5,1)),
    treatment = factor(treatment,
                       levels = rev(treats_ord)),
    prec = 1/se
  ) 

# Plot: obs first, then model; axis breaks at 1:length(treats_ord)
ggplot() +
  geom_errorbarh(
    data   = obs_df,
    aes(xmin = lo, xmax = hi, y = y_num),
    colour = "grey70", height = 0
  ) +
  geom_point(data = obs_df,
             aes(x = est, y = y_num,
                 colour = treatment,
                 size   = prec,
                 alpha = alpha_val)
  ) +
  scale_size_continuous(range = c(0.8, 2.5), name = "", guide = "none") +
  geom_text(
    data = obs_df %>% filter(!is.na(Q_het_decomp)),
    aes(
      x     = est,
      y     = y_num,
      label = sprintf("%.2f", Q_het_decomp)
    ),
    hjust = -0.2,
    size = 2.8
  ) +
  geom_pointrange(
    data   = df_plot,
    aes(x = est, xmin = lo, xmax = hi, y = y_num, shape = model),
    colour = "black", size = 0.5
  ) +
  scale_colour_brewer(palette = "Set2", name = "", guide = "none") +
  scale_shape_manual(
    values = c(RE = 15, ME = 17),
    breaks = c("RE","ME"),
    name   = ""
  ) +
  scale_alpha_continuous(range = c(0.5, 1), guide = "none") +
  scale_y_continuous(
    breaks = seq_along(treats_ord),
    labels = treats_ord,
    expand = expansion(add = 0.6)
  ) +
  labs(x = "Log Risk Ratio (vs. Placebo)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(1, "lines")
  )

#=========Q.het decomposition===========
df <- tibble(
  y    = -net$TE,
  se  = net$seTE,
  design = paste(treatment[as.numeric(net$treat1)],
                  treatment[as.numeric(net$treat2)], sep="-")
) 

# design-specific estimates
df <- df %>%
  group_by(design) %>%
  mutate(thetaD = sum(y/(se^2)) / sum(1/(se^2))) %>%
  ungroup()  %>%
  mutate(Q_het_decomp = (y - thetaD)^2 / (se^2)) %>%
  mutate(thetaFE = net$TE.nma.fixed)

df_02 <- df %>% 
  select(design, Q_het_decomp, y, se) %>% arrange(design)

print(df_02, n=29)

df_02 %>%
  group_by(design) %>%
  summarise(sum_Q_het = sum(Q_het_decomp, na.rm = TRUE))

# #==============remove study with Q_het=23 =======================

net_data <- as.data.frame(cbind(net$treat1, net$treat2))
colnames(net_data) <- c("treat1", "treat2")
net_data$TE <- net$TE
net_data$seTE <- net$seTE

net_1 <- netmeta(TE, seTE, treat1, treat2, data = net_data[-26, ])


theta <- net_1$TE
m <- net_1$m
n <- net_1$n
V <- diag(net_1$seTE^2)

# Multiplicative effect
theta_me <- net_1$TE.nma.fixed # fitted TE
phi <- as.numeric(t(theta-theta_me) %*%
                    solve(V) %*% (theta-theta_me) / (m-n+1))
logL_me <- -0.5*(m*log(2*pi)+
                   log(det(phi*V))+
                   t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
AIC_mul <- 2*n - 2*logL_me

# Additive effects
tau_hat <- net_1$tau
sig <- V+tau_hat^2*diag(m)
theta_ae <- net_1$TE.nma.random
logL_ae <- -0.5*(m*log(2*pi)+
                   log(det(sig))+
                   t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
AIC_add <- 2*n-2*logL_ae

phi
tau_hat

AIC_mul-AIC_add

net_1$Q.heterogeneity
net_1$pval.Q.heterogeneity