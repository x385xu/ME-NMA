library(netmeta)
library(nmadb)
library(dplyr)
library(tidyr)
library(ggplot2)
#3, recid = 473552
dat_nmadb[3,]
net <- runnetmeta(473552)
read <- readByID(473552)

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

dat_nmadb[ind, ]
summary(net)
net$TE

# [1] -2.0250540 -2.3225004 -2.1887301 -2.0034024 -0.7777473 -1.1392833
# [7] -1.3217558 -0.1673743 -3.1921727 -0.6947491 -1.2169390 -1.7002937
# [13] -2.4784566 -2.0118214 -1.5066884 -2.5041175 -0.3673645 -2.4401923
# [19] -0.8721769 -2.4086121  0.5042706 -1.5416250 -1.7811979 -2.9714973
# [25] -2.6021630 -1.5889394 -2.2983207  1.2171428 -1.1026296 -0.1115178
# [31] -1.3837912 -1.3362839
#=========netgraph=================
read$data[c(41,42, 55,56),]
treatment <- c("Placebo/Standard Care", "Abatacept",
               "Adalimumab", "Certolizumab",
               "Etanercept", "Golimumab",
               "Infliximab", "Rituximab", "Tocilizumab")

netgraph(net,
         labels = treatment,
         plastic = FALSE, 
         iterate = FALSE,
         number.of.studies = TRUE, 
         col="black", 
         cex=1, 
         multiarm=FALSE, 
         points=TRUE, 
         col.points="black", 
         cex.points=3, 
         cex.number.of.studies = 1,
         pos.number.of.studies = 0.3,
         lwd = 3,
         offset = 0.04,
         #         srt.labels = -20,
         #         rotate = -45
)

#===========forest plot====================


# net_data <- as.data.frame(cbind(net$treat1, net$treat2))
# colnames(net_data) <- c("treat1", "treat2")
# net_data$TE <- net$TE
# net_data$seTE <- net$seTE
# 
# net <- netmeta(TE, seTE, treat1, treat2, data = net_data)

d_mul <- net$TE.fixed[-1, 1]
se_mul <- sqrt(phi)*net$seTE.fixed[-1, 1]
ci_mul_upper <- d_mul + qnorm(0.975)*se_mul
ci_mul_lower <- d_mul + qnorm(0.025)*se_mul

d_add <- net$TE.random[-1, 1]
ci_add_lower <- net$lower.random[-1,1 ]
ci_add_upper <- net$upper.random[-1,1 ]

treatment_full <- treatment

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
  Q_het_decomp = df_04$Q_het_decomp,
) %>%
  filter(!is.na(est), !is.na(se)) %>%
  mutate(
    treatment = treatment_full[trt_code],
    lo        = est + qnorm(0.025)*se,
    hi        = est + qnorm(0.975)*se,
    treatment = factor(treatment, levels = treats_ord),
    y_base    = as.numeric(treatment),
    Q_het_decomp = ifelse(Q_het_decomp >= net$Q.heterogeneity/n, 
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
  scale_size_continuous(range = c(0.8, 2.3), name = "", guide = "none") +
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
  labs(x = "Log Odds Ratio (vs. Placebo/Standard Care)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(1, "lines")
  )
#========Q decomp===============

df <- tibble(
  y    = -net$TE,
  se  = net$seTE,
  design = paste(treatment[as.numeric(net$treat1)],
                 treatment[as.numeric(net$treat2)], sep="-"),
  treatment1 = paste(treatment[as.numeric(net$treat1)]),
  treatment2 = paste(treatment[as.numeric(net$treat2)])
) %>%
  group_by(design) %>%
  mutate(thetaD = sum(y/(se^2)) / sum(1/(se^2))) %>%
  ungroup() %>%
  mutate(Q_het_decomp = (y - thetaD)^2 / (se^2)) %>%
  mutate(thetaFE = net$TE.nma.fixed)



df_04 <- df %>% select(treatment1, treatment2, Q_het_decomp, y, se)
print(df_04, n=32)

Q_total <- sum((df$y - df$thetaFE)^2 / df$var)


#==============remove study with negative log OR =======================

net_data <- as.data.frame(cbind(net$treat1, net$treat2))
colnames(net_data) <- c("treat1", "treat2")
net_data$TE <- net$TE
net_data$seTE <- net$seTE

net_1 <- netmeta(TE, seTE, treat1, treat2, data = net_data[-c(21, 28), ])


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

net<-net_1
