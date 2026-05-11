#' @description
#' Compute AIC values for fixed-effect, additive random-effects, and multiplicative models
#'
#' This function fits network meta-analysis models to selected datasets and computes
#' AIC values under three variance structures: fixed-effect (FE), additive
#' random-effects (RE), and multiplicative-error (ME). It also extracts
#' heterogeneity, inconsistency, and I-squared statistics from the fitted
#' \code{netmeta} objects.
#'
#' @param dat Data frame containing the network meta-analysis database. Must include
#'   a \code{recid} column used by \code{runnetmeta()} to retrieve each network.
#'   Defaults to \code{dat_nmadb}.
#' @param idx_nmadb Data frame indexing the networks to be analyzed. Must include
#'   columns \code{index}, \code{effect_measure}, and \code{twoarm_only}.
#' @param measure Character vector specifying which effect measures to include.
#'   Defaults to \code{c("odds ratio", "risk ratio", "mean difference")}.
#' @param twoarm_only Logical. If \code{TRUE}, only networks consisting exclusively
#'   of two-arm studies are included. If \code{FALSE}, only networks containing
#'   multi-arm studies are included.
#' @param method.tau Character string specifying the estimator used for the
#'   between-study standard deviation \eqn{\tau} in \code{netmeta}. Defaults to
#'   \code{"DL"}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{index}: Indices of the analyzed networks.
#'   \item \code{recid}: Record IDs of the analyzed networks.
#'   \item \code{effect_measure}: Effect measure for each network.
#'   \item \code{tau}: Estimated between-study heterogeneity standard deviation.
#'   \item \code{phi}: Estimated multiplicative variance inflation factor.
#'   \item \code{AIC_re}: AIC under the additive random-effects model.
#'   \item \code{AIC_fe}: AIC under the fixed-effect model.
#'   \item \code{AIC_me}: AIC under the multiplicative-error model.
#'   \item \code{me_re}: Difference \code{AIC_me - AIC_re}.
#'   \item \code{me_fe}: Difference \code{AIC_me - AIC_fe}.
#'   \item \code{Q}, \code{Q_pval}: Total Q statistic and corresponding p-value.
#'   \item \code{Qh}, \code{Qh_pval}: Heterogeneity Q statistic and p-value.
#'   \item \code{Qi}, \code{Qi_pval}: Inconsistency Q statistic and p-value.
#'   \item \code{I2}, \code{I2_lower}, \code{I2_upper}: I-squared estimate and
#'     confidence interval limits.
#' }
#'

compute_AIC <- function(dat = dat_nmadb,
                        idx_nmadb,
                        measure = c("odds ratio", "risk ratio", "mean difference"),
                        twoarm_only = TRUE,
                        method.tau = "DL") {
  
  idx_use <- idx_nmadb %>%
    dplyr::filter(
      effect_measure %in% measure,
      twoarm_only == !!twoarm_only
    )
  
  idx <- idx_use$index
  len <- length(idx)
  
  # Preallocate
  AIC_re <- AIC_fe <- AIC_me <- numeric(len)
  taus <- phis <- numeric(len)
  Q <- Q_pval <- numeric(len)
  I2 <- I2_lower <- I2_upper <- numeric(len)
  Qh <- Qh_pval <- Qi <- Qi_pval <- numeric(len)
  
  # helper
  safe_num <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }
  
  for (j in seq_along(idx)) {
    
    i <- idx[j]
    
    net <- tryCatch(runnetmeta(dat$recid[i]), error = function(e) NULL)
    if (!is.list(net)) next
    
    net_data <- data.frame(
      treat1 = net$treat1,
      treat2 = net$treat2,
      TE = net$TE,
      seTE = net$seTE
    )
    
    net_fit <- tryCatch(
      netmeta(TE, seTE, treat1, treat2, data = net_data, method.tau = method.tau),
      error = function(e) {
        message("netmeta failed for recid ", dat$recid[i], ": ", e$message)
        NULL
      }
    )
    
    if (is.null(net_fit)) next
    
    theta <- net_fit$TE
    V <- diag(net_fit$seTE^2)
    m <- net_fit$m
    n <- net_fit$n
    theta_fe <- net_fit$TE.nma.fixed
    
    Q_total <- net_fit$Q
    phi_hat <- max(1, Q_total / (m - n + 1))
    
    logL_me <- -0.5 * (
      m * log(2 * pi) +
        log(det(phi_hat * V)) +
        t(theta - theta_fe) %*% solve(phi_hat * V) %*% (theta - theta_fe)
    )
    AIC_me[j] <- 2 * n - 2 * logL_me
    
    tau_hat <- net_fit$tau
    Sigma_re <- V + tau_hat^2 * diag(m)
    theta_re <- net_fit$TE.nma.random
    
    logL_re <- -0.5 * (
      m * log(2 * pi) +
        log(det(Sigma_re)) +
        t(theta - theta_re) %*% solve(Sigma_re) %*% (theta - theta_re)
    )
    AIC_re[j] <- 2 * n - 2 * logL_re
    
    logL_fe <- -0.5 * (
      m * log(2 * pi) +
        log(det(V)) +
        t(theta - theta_fe) %*% solve(V) %*% (theta - theta_fe)
    )
    AIC_fe[j] <- 2 * n - 2 * logL_fe
    
    
    taus[j] <- safe_num(tau_hat)
    phis[j] <- safe_num(phi_hat)
    
    Q[j] <- safe_num(net_fit$Q)
    Q_pval[j] <- safe_num(net_fit$pval.Q)
    
    Qh[j] <- safe_num(net_fit$Q.heterogeneity)
    Qh_pval[j] <- safe_num(net_fit$pval.Q.heterogeneity)
    
    Qi[j] <- safe_num(net_fit$Q.inconsistency)
    Qi_pval[j] <- safe_num(net_fit$pval.Q.inconsistency)
    
    I2[j] <- safe_num(net_fit$I2)
    I2_lower[j] <- safe_num(net_fit$lower.I2)
    I2_upper[j] <- safe_num(net_fit$upper.I2)
  }
  
  return(list(
    index = idx,
    recid = dat$recid[idx],
    effect_measure = idx_use$effect_measure,
    tau = taus,
    phi = phis,
    AIC_re = AIC_re,
    AIC_fe = AIC_fe,
    AIC_me = AIC_me,
    me_re = AIC_me - AIC_re,
    me_fe = AIC_me - AIC_fe,
    Q = Q,
    Q_pval = Q_pval,
    Qh = Qh,
    Qh_pval = Qh_pval,
    Qi = Qi,
    Qi_pval = Qi_pval,
    I2 = I2,
    I2_lower = I2_lower,
    I2_upper = I2_upper
  ))
}