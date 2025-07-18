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
# where phi < 1
# or-2: 6,12,18,26
# ind <- AIC_or$ind[26]
# rr-2: 17,24,29
# ind <- AIC_rr$ind[29]

# min(AIC_md$mul_add)= -32.60178 (MD)
#ind <- AIC_md$ind[which(AIC_md$mul_add == min(AIC_md$mul_add))]

# max(AIC_or$mul_add)= 11.46659 (OR)
#ind <- AIC_or$ind[which(AIC_or$mul_add == max(AIC_or$mul_add))]

#ind <- AIC_rr$ind[which(AIC_rr$mul_fixed == max(AIC_rr$mul_fixed))]
#ind <- AIC_rr$ind[which(AIC_rr$mul_add == min(AIC_rr$mul_add))]

net <- runnetmeta(dat_nmadb$recid[ind])

netgraph(net, plastic = FALSE, iterate=FALSE,
         number.of.studies = TRUE, col="darkgray", cex=1.5, 
         multiarm=FALSE, points=TRUE, col.points="black", 
         cex.points=3, allfigures=TRUE)

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
theta_fe <- net$TE.nma.fixed
logL_fe <- -0.5*(m*log(2*pi)+
                   log(det(V))+ 
                   t(theta-theta_fe) %*% solve(V) %*% (theta-theta_fe))
AIC_fixed <- 2*n-2*logL_fe
# Additive effects
tau_hat <- net$tau
sig <- V+tau_hat^2*diag(m)
theta_ae <- net$TE.nma.random
logL_ae <- -0.5*(m*log(2*pi)+
                   log(det(sig))+ 
                   t(theta-theta_ae) %*% solve(sig) %*% (theta-theta_ae))
AIC_add <- 2*n-2*logL_ae

dat_nmadb[ind, ]
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

df_plot <- data.frame(
  treatment = factor(rep(2:(length(d_mul)+1), times = 2)),
  model     = rep(c( "ME", "RE"), each = length(d_mul)),
  est  = c(d_mul, d_add),
  lo   = c(ci_mul_lower, ci_add_lower),
  hi   = c(ci_mul_upper, ci_add_upper)
)

ggplot(df_plot, aes(y = treatment, x = est, color = model)) +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5),
                  size = 0.1) +
  geom_point(position = position_dodge(width = 0.5), 
             size = 2, shape = 15, alpha = 0.5) +
  labs(y = "Treatment", x = "Treatment Effect", color = "Model") +
  theme_minimal(base_size = 14)



#==========min(AIC_md$mul_add)= -32.60178 (MD)=======================
readByID(480039)

d_mul <- net$TE.fixed[1, ]
se_mul   <- sqrt(phi)*net$seTE.fixed[1, ]
ci_mul_upper   <- d_mul + qnorm(0.975)*se_mul
ci_mul_lower   <- d_mul + qnorm(0.025)*se_mul

d_add <- net$TE.random[1, ]
ci_add_lower  <- net$lower.random[1, ]
ci_add_upper  <- net$upper.random[1, ]

d_bayesian <- c(0, 14.71, 37.63, 22.00, 32.53, 20.17, 30.71, 31.28)
ci_bayesian_lower <- c(0, -3.85, 6.71, 0.86, 13.46, 12.33, 15.14, 18.69)
ci_bayesian_upper <- c(0, 33.43, 67.22, 42.52, 52.09, 29.73, 46.97, 45.21)

df_plot <- data.frame(
  treatment = factor(rep(1:8, times = 3)),
  model     = rep(c("Bayesian", "Multiplicative", "Random"), each = 8),
  est  = c(d_bayesian, d_mul, d_add),
  lo   = c(ci_bayesian_lower, ci_mul_lower, ci_add_lower),
  hi   = c(ci_bayesian_upper, ci_mul_upper, ci_add_upper)
)

ggplot(df_plot, aes(x = treatment, y = est, color = model)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.6),
                  size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x     = "Treatment",
    y     = "TE",
    color = "Model",
    title = "Comparison of Treatment Effects Across Models"
  ) +
  theme_minimal(base_size = 14)

