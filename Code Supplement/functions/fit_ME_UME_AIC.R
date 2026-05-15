#' @description
#' Fit a multiplicative unrelated mean effects model and compute its AIC.
#'
#' This function fits an unrelated mean effects (UME) model using the
#' study-level pairwise comparisons stored in a fitted \code{netmeta}
#' object. The UME model relaxes the consistency constraints of a standard
#' network meta-analysis by estimating a separate mean effect for each
#' observed direct comparison.
#'
#' @param net A fitted \code{netmeta} object containing study-level
#'   pairwise treatment contrasts. The object must include the components
#'   \code{TE}, \code{seTE}, \code{treat1}, \code{treat2}, and
#'   \code{studlab}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{estimates}: Data frame of UME estimates, with one
#'     estimated mean effect for each observed direct comparison.
#'   \item \code{loglik}: Maximized marginal log-likelihood under the
#'     ME-UME model.
#'   \item \code{k}: Number of estimated parameters, equal to the number
#'     of direct-comparison mean effects plus one multiplicative
#'     variance inflation parameter.
#'   \item \code{AIC}: Akaike information criterion under the ME-UME
#'     model.
#'   \item \code{phi}: Estimated multiplicative variance inflation
#'     factor.
#'   \item \code{Q}: Residual weighted sum of squares under the UME
#'     model.
#'   \item \code{df}: Residual degrees of freedom used to estimate
#'     \code{phi}.
#' }
#' 
fit_ME_UME_AIC <- function(net) {
  
  marginal_loglik <- function(y, se, mu, phi) {
    sum(dnorm(y, mean = mu, sd = sqrt(phi) * se, log = TRUE))
  }
  
  # ── 1. Build oriented observation-level data ───────────────────────────────
  dat <- data.frame(
    TE      = net$TE,
    seTE    = net$seTE,
    treat1  = as.character(net$treat1),
    treat2  = as.character(net$treat2),
    studlab = as.character(net$studlab)
  )
  
  dat <- dat[complete.cases(dat[, c("TE", "seTE", "treat1", "treat2")]), ]
  
  dat$treat_low  <- pmin(dat$treat1, dat$treat2)
  dat$treat_high <- pmax(dat$treat1, dat$treat2)
  dat$sign       <- ifelse(dat$treat1 == dat$treat_low, 1, -1)
  dat$y          <- dat$sign * dat$TE
  dat$comp       <- paste(dat$treat_low, dat$treat_high, sep = " vs ")
  
  # ── 2. UME fitted values under ME weighting ────────────────────────────────
  # Since phi is common to all observations, it cancels out of the WLS estimate.
  comps  <- unique(dat$comp)
  mu_ume <- numeric(nrow(dat))
  
  for (comp_i in comps) {
    idx         <- which(dat$comp == comp_i)
    w_c         <- 1 / dat$seTE[idx]^2
    mu_ume[idx] <- sum(w_c * dat$y[idx]) / sum(w_c)
  }
  
  # ── 3. Restrict to direct-evidence comparisons only ───────────────────────
  sp        <- netsplit(net)
  direct_df <- as.data.frame(sp$direct.random)
  direct_df <- direct_df[!is.na(direct_df$TE), ]
  
  split_comp     <- strsplit(as.character(direct_df$comparison), ":")
  direct_df$comp <- paste(
    sapply(split_comp, function(x) pmin(trimws(x[1]), trimws(x[2]))),
    sapply(split_comp, function(x) pmax(trimws(x[1]), trimws(x[2]))),
    sep = " vs "
  )
  
  direct_obs    <- dat$comp %in% direct_df$comp
  dat_direct    <- dat[direct_obs, ]
  mu_ume_direct <- mu_ume[direct_obs]
  
  # ── 4. Estimate multiplicative heterogeneity phi ──────────────────────────
  Q_ume <- sum((dat_direct$y - mu_ume_direct)^2 / dat_direct$seTE^2)
  df_ume <- nrow(dat_direct) - length(unique(dat_direct$comp))
  
  phi_hat <- max(1, Q_ume / df_ume)
  
  # ── 5. UME log-likelihood and AIC ─────────────────────────────────────────
  loglik_ume <- marginal_loglik(
    y   = dat_direct$y,
    se  = dat_direct$seTE,
    mu  = mu_ume_direct,
    phi = phi_hat
  )
  
  k_ume   <- length(unique(dat_direct$comp)) + 1  # means + phi
  AIC_ume <- -2 * loglik_ume + 2 * k_ume
  
  # ── 6. Per-comparison estimate table ──────────────────────────────────────
  estimates <- aggregate(
    mu_ume_direct,
    by = list(comp = dat_direct$comp),
    FUN = mean
  )
  names(estimates)[2] <- "mu_ume"
  
  list(
    estimates = estimates,
    loglik    = loglik_ume,
    k         = k_ume,
    AIC       = AIC_ume,
    phi       = phi_hat,
    Q         = Q_ume,
    df        = df_ume
  )
}
