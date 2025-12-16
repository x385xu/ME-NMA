library(netmeta)
library(nmadb)
library(dplyr)
library(tidyr)
library(ggplot2)

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

source("get_index_nmadb.R")
source("compute_AIC.R")


#=====================================================================

#59, recid=480039
dat_nmadb[59, ]
net <- runnetmeta(480039)
read <- readByID(480039)


treatment <- c("Placebo", "MTX", "Abatacept+MTX", "Anakinra+MTX",
               "aTNF+MTX", "aTNF", "Tocilizaumab+MTX", 
               "Tocilizumab")

net$seTE

netgraph(net,
         labels = treatment,
         plastic = FALSE, 
         iterate=TRUE,
         number.of.studies = TRUE, 
         col="black", 
         cex=1, 
         multiarm=FALSE, 
         points=TRUE, 
         col.points="black", 
         cex.points=2, 
         cex.number.of.studies = 1,
         offset = 0.03,
         adj = matrix(c(0, -0.5, 0, 0.5, 1, -0.4, -0.1, -0.2, 
                        0.5, 0.3, 0, -2, 1, 0.3, 0.3, 0.3),
                      nrow = 8),
         lwd = 2,
#         srt.labels = -20,
         rotate = 40
)

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

# Fixed Effect
# theta_fe <- net$TE.nma.fixed
# logL_fe <- -0.5*(m*log(2*pi)+
#                    log(det(V))+ 
#                    t(theta-theta_fe) %*% solve(V) %*% (theta-theta_fe))
# AIC_fixed <- 2*n-2*logL_fe

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



#------------Compare Observed TE with fitted TE----------------------------------------------------------
# multiplicative d is the same as the fixed effects model

# Residual
net$TE.nma.fixed-net$TE
net$TE.nma.random-net$TE

#---------Residual Plot------------------
df_resid <- tibble(
  comparison = seq_along(net$TE),
  ObsTE = net$TE,
  RE = net$TE.nma.random - net$TE,
  ME = net$TE.nma.fixed - net$TE
) %>%
  pivot_longer(
    cols      = c(RE, ME),
    names_to  = "Model",
    values_to = "Residual")


ggplot(df_resid, aes(x = ObsTE, y = Residual, color = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = Model), size = 2, alpha = 0.8) +
  labs(x = "Observed Treatment Effect", y = "Residual") +
  theme_minimal(base_size = 12)


#----------TE & CI----------------------------------------------------
d_mul <- net$TE.fixed[1, -1]
se_mul <- sqrt(phi)*net$seTE.fixed[1, -1]
ci_mul_upper <- d_mul + qnorm(0.975)*se_mul
ci_mul_lower <- d_mul + qnorm(0.025)*se_mul

d_add <- net$TE.random[1,-1 ]
ci_add_lower <- net$lower.random[1,-1 ]
ci_add_upper <- net$upper.random[1,-1 ]



#==============correct seTE================================
# the correct seTE
sd  <- read1$data[, "sd"]

odd <- seq(from = 1, to = 25, by = 2)
even <- seq(from = 2, to = 26, by = 2)

se <- sqrt(sd[odd]^2+sd[even]^2)
net_data <- as.data.frame(cbind(net$treat1, net$treat2))
colnames(net_data) <- c("treat1", "treat2")
net_data$TE <- net$TE
net_data$seTE <- se

net_2 <- netmeta(TE, seTE, treat1, treat2, data = net_data)


theta <- net_2$TE
m <- net_2$m
n <- net_2$n
V <- diag(net_2$seTE^2)

# Multiplicative effect
theta_me <- net_2$TE.nma.fixed # fitted TE
phi <- as.numeric(t(theta-theta_me) %*% 
                    solve(V) %*% (theta-theta_me) / (m-n+1))
logL_me <- -0.5*(m*log(2*pi)+
                   log(det(phi*V))+ 
                   t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
AIC_mul <- 2*n - 2*logL_me

# Additive effects
tau_hat <- net_2$tau
sig <- V+tau_hat^2*diag(m)
theta_ae <- net_2$TE.nma.random
logL_ae <- -0.5*(m*log(2*pi)+
                   log(det(sig))+ 
                   t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
AIC_add <- 2*n-2*logL_ae

phi
tau_hat

AIC_mul-AIC_add

net_2$Q.heterogeneity
net_2$pval.Q.heterogeneity

#--------forest plot with corrected seTE-------------
d_mul <- net_2$TE.fixed[1, -1]
se_mul <- sqrt(phi)*net_2$seTE.fixed[1, -1]
ci_mul_upper <- d_mul + qnorm(0.975)*se_mul
ci_mul_lower <- d_mul + qnorm(0.025)*se_mul

d_add <- net_2$TE.random[1,-1 ]
ci_add_lower <- net_2$lower.random[1,-1 ]
ci_add_upper <- net_2$upper.random[1,-1 ]

treatment_full <- treatment
# drop "Placebo"
treats <- treatment_full[-1]
order_idx  <- order(d_add, decreasing = TRUE)
treats_ord <- treats[order_idx]

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
    y_num     = y_base + if_else(model=="RE", +0.09, -0.09),
    y_num = y_num +if_else(treatment == "aTNF", +0.2, 0)
  ) 

trt1 <- (net$treat1 == 1)
obs_df <- tibble(
  trt_code = as.integer(net_2$treat2),
  est      = net_2$TE*trt1,
  se       = net_2$seTE*trt1
) %>%
  filter(!is.na(est), !is.na(se)) %>%
  mutate(
    treatment = treatment_full[trt_code],
    lo        = est + qnorm(0.025)*se,
    hi        = est + qnorm(0.975)*se,
    treatment = factor(treatment, levels = treats_ord),
    y_base    = as.numeric(treatment)
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
  labs(x = "Mean Difference (vs. Placebo)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position   = "right",
    legend.key.height = unit(1, "lines")
  )

#============Q decomp====================
df <- tibble(
  y    = net$TE,
  var  = net$seTE^2,
  #  design = paste(net$treat1, net$treat2, sep="-"),
  design = paste(treatment[as.numeric(net$treat1)],
                treatment[as.numeric(net$treat2)], sep=" - ")
) %>%
  group_by(design) %>%
  mutate(thetaD = sum(y/var) / sum(1/var)) %>%
  ungroup() %>%
  mutate(Q_het_decomp = (df$y - df$thetaD)^2 / df$var) %>%
  mutate(thetaFE = net$TE.nma.fixed)

df_01 <- df %>% 
  select(design, Q_het_decomp, y, var)

# #==============remove small study=======================
# 
# net_data <- as.data.frame(cbind(net$treat1, net$treat2))
# colnames(net_data) <- c("treat1", "treat2")
# net_data$TE <- net$TE
# net_data$seTE <- net$seTE
# 
# net_1 <- netmeta(TE, seTE, treat1, treat2, data = net_data[-10, ])
# 
# 
# theta <- net_1$TE
# m <- net_1$m
# n <- net_1$n
# V <- diag(net_1$seTE^2)
# 
# # Multiplicative effect
# theta_me <- net_1$TE.nma.fixed # fitted TE
# phi <- as.numeric(t(theta-theta_me) %*% 
#                     solve(V) %*% (theta-theta_me) / (m-n+1))
# logL_me <- -0.5*(m*log(2*pi)+
#                    log(det(phi*V))+ 
#                    t(theta-theta_me) %*% solve(phi*V) %*% (theta-theta_me))
# AIC_mul <- 2*n - 2*logL_me
# 
# # Additive effects
# tau_hat <- net_1$tau
# sig <- V+tau_hat^2*diag(m)
# theta_ae <- net_1$TE.nma.random
# logL_ae <- -0.5*(m*log(2*pi)+
#                    log(det(sig))+ 
#                    t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
# AIC_add <- 2*n-2*logL_ae
# 
# phi
# tau_hat
# 
# AIC_mul-AIC_add
# 
# net_1$Q.heterogeneity
# net_1$pval.Q.heterogeneity
# 
# #======-32.60178 (MD) forest plot with estimates in original paper=======================
# readByID(480039)
# 
# d_mul <- net_2$TE.fixed[1, ]
# se_mul   <- sqrt(phi)*net_2$seTE.fixed[1, ]
# ci_mul_upper   <- d_mul + qnorm(0.975)*se_mul
# ci_mul_lower   <- d_mul + qnorm(0.025)*se_mul
# 
# d_add <- net_2$TE.random[1, ]
# ci_add_lower  <- net_2$lower.random[1, ]
# ci_add_upper  <- net_2$upper.random[1, ]
# 
# d_bayesian <- c(0, 14.71, 37.63, 22.00, 32.53, 20.17, 30.71, 31.28)
# ci_bayesian_lower <- c(0, -3.85, 6.71, 0.86, 13.46, 12.33, 15.14, 18.69)
# ci_bayesian_upper <- c(0, 33.43, 67.22, 42.52, 52.09, 29.73, 46.97, 45.21)
# 
# df_plot <- data.frame(
#   treatment = factor(rep(1:8, times = 3)),
#   model     = rep(c("Bayesian", "Multiplicative", "Random"), each = 8),
#   est  = c(d_bayesian, d_mul, d_add),
#   lo   = c(ci_bayesian_lower, ci_mul_lower, ci_add_lower),
#   hi   = c(ci_bayesian_upper, ci_mul_upper, ci_add_upper)
# ) 
# 
# ggplot(df_plot, aes(x = est, y = treatment, color = model)) +
#   geom_pointrange(aes(xmin = lo, xmax = hi),
#                   position = position_dodge(width = 0.6),
#                   size = 0.8) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(
#     x     = "MD",
#     y     = "Treatment",
#     color = "Model",
#     title = "Comparison of Treatment Effects Across Models"
#   ) +
#   theme_minimal(base_size = 14)
# 
