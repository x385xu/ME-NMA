library(netmeta)
library(nmadb)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 276, recid = 501235
dat_nmadb[276,]
net <- runnetmeta(501235)
read <- readByID(501235)

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
# net$TE
# [1]  0.05064373 -0.57738796 -0.30818059 -0.05540604 -0.50962231 -1.83651737 -0.16653474
# [8] -2.99187870 -2.77258872 -0.10536052 -0.60549919 -1.97311984 -1.10454702 -0.83034830
# [15] -2.29253476  0.20067070  0.00000000  3.14595058 -0.01675702 -1.57304431
#=========netgraph=================
read <- readByID(dat_nmadb$recid[ind])

#LCFE = low-cost/free equipment
#HSI = home safety inspection
treatment <- c("Usual Care", "Education",
               "Education+LCFE",
               "Education+LCFE+HSI",
               "Education+LCFE+Fitting",
               "Education+HSI",
               "Education+LCFE+Fitting+HSI")

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
    y_num     = y_base + if_else(model=="RE", +0.28, +0.1),
    y_num = if_else(y_base == 5, y_num-0.19, y_num)
  ) 

trt1 <- (net$treat1 == 1)

obs_df <- tibble(
  trt_code = as.integer(net$treat2),
  est      = -net$TE*trt1,
  se       = net$seTE*trt1,
  Q_het_decomp = df_03$Q_het_decomp,
) %>%
  filter(!is.na(est), !is.na(se)) %>%
  mutate(
    treatment = treatment_full[trt_code],
    lo        = est + qnorm(0.025)*se,
    hi        = est + qnorm(0.975)*se,
    treatment = factor(treatment, levels = treats_ord),
    y_base    = as.numeric(treatment),
    Q_het_decomp = ifelse(trt1 == 1, Q_het_decomp, NA),
    Q_het_decomp = ifelse(Q_het_decomp >= 1.2, Q_het_decomp, NA)
  ) %>%
  group_by(treatment) %>%
  arrange(se, .by_group = TRUE) %>%
  mutate(
    n_obs = n(),
    idx   = row_number(),
    # keep them within ±0.2 of the integer row
    y_num = y_base + (idx - (n_obs+1)/2)*(0.4/n_obs) - 0.2
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
  labs(x = "Log Odds Ratio (vs. Usual Care)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(1, "lines")
  )
  
#===========Q statistic==========================

# example in R for design "1-2"
# i <- which(net$treat1==1 & net$treat2==2)
# wi <- 1/(net$seTE[i]^2)
# thetaD_12 <- sum(wi * net$TE[i]) / sum(wi)

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



df_03 <- df %>% select(treatment1, treatment2, Q_het_decomp, y, se)
print(df_03)

df %>%
  group_by(design) %>%
  summarise(sum_Q_het = sum(Q_het_decomp, na.rm = TRUE))
df %>%
  summarise(sum_Q_het = sum(Q_het_decomp, na.rm = TRUE))

Q_total <- sum((df$y - df$thetaFE)^2 / ((df$se)^2))
#> Q_total == 36.12

#==============remove outliers=======================

net_data <- as.data.frame(cbind(net$treat1, net$treat2))
colnames(net_data) <- c("treat1", "treat2")
net_data$TE <- net$TE
net_data$seTE <- net$seTE

net_1 <- netmeta(TE, seTE, treat1, treat2, data = net_data[-c(8,9,18), ])


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


#==============only include RCTs=======================
rct<-c(1,0,1,0,0,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1)


net_data <- as.data.frame(cbind(net$treat1, net$treat2))
colnames(net_data) <- c("treat1", "treat2")
net_data$TE <- net$TE
net_data$seTE <- net$seTE

net_1 <- netmeta(TE, seTE, treat1, treat2, data = net_data[which(rct==1), ])


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
