#' @description
#' Compute AIC values for common-effect, additive random-effects, and
#' multiplicative models for network meta-analysis.
#'
#' This function fits network meta-analysis models using locally saved
#' datasets and computes AIC values under three variance structures:
#' common-effect (CE), additive random-effects (RE), and
#' multiplicative-effects (ME). It also extracts heterogeneity,
#' inconsistency, and I-squared statistics from the fitted
#' \code{netmeta} objects.
#'
#' @param twoarm_data_list Named list of locally saved NMA datasets.
#'   Each element should correspond to one network indexed by its
#'   \code{recid}.
#' @param idx_nmadb Data frame indexing the networks to be analyzed.
#'   Must include columns \code{index}, \code{recid},
#'   \code{effect_measure}, and \code{twoarm_only}.
#' @param measure Character vector specifying which effect measures
#'   to include. Defaults to
#'   \code{c("odds ratio", "risk ratio", "mean difference")}.
#' @param twoarm_only Logical. If \code{TRUE}, only networks
#'   consisting exclusively of two-arm studies are included.
#'   If \code{FALSE}, only networks containing multi-arm studies
#'   are included.
#' @param method.tau Character string specifying the estimator used
#'   for the between-study variance \eqn{\tau^2} in
#'   \code{netmeta}. Defaults to \code{"DL"}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{index}: Indices of the analyzed networks.
#'   \item \code{recid}: Record IDs of the analyzed networks.
#'   \item \code{effect_measure}: Effect measure for each network.
#'   \item \code{tau2}: Estimated between-study heterogeneity
#'     variance.
#'   \item \code{phi}: Estimated multiplicative variance inflation
#'     factor.
#'   \item \code{AIC_re}: AIC under the additive random-effects model.
#'   \item \code{AIC_ce}: AIC under the common-effect model.
#'   \item \code{AIC_me}: AIC under the multiplicative-error model.
#'   \item \code{me_re}: Difference \code{AIC_me - AIC_re}.
#'   \item \code{me_ce}: Difference \code{AIC_me - AIC_ce}.
#'   \item \code{Q}, \code{Q_pval}: Total Q statistic and
#'     corresponding p-value.
#'   \item \code{Qh}, \code{Qh_pval}: Heterogeneity Q statistic
#'     and p-value.
#'   \item \code{Qi}, \code{Qi_pval}: Inconsistency Q statistic
#'     and p-value.
#'   \item \code{I2}, \code{I2_lower}, \code{I2_upper}:
#'     I-squared estimate and confidence interval limits.
#' }

compute_AIC <- function(twoarm_data_list,
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
  AIC_re <- AIC_ce <- AIC_me <- numeric(len)
  taus <- phis <- numeric(len)
  Q <- Q_pval <- numeric(len)
  I2 <- I2_lower <- I2_upper <- numeric(len)
  Qh <- Qh_pval <- Qi <- Qi_pval <- numeric(len)
  
  for (j in seq_along(idx_use$recid)) {
    
    recid_j <- as.character(idx_use$recid[j])
    dat_i <- twoarm_data_list[[recid_j]]
    
    if (is.null(dat_i)) {
      message("Dataset not found for recid ", recid_j)
      next
    }
    
    net_fit <- tryCatch(
      fit_netmeta(
        dat_i,
        model = "random",
        method.tau = method.tau
      ),
      error = function(e) {
        message("netmeta failed for recid ", recid_j, ": ", e$message)
        NULL
      }
    )
    
    if (is.null(net_fit)) next
    
    # compute AIC for ME model
    theta <- net_fit$TE
    V <- diag(net_fit$seTE^2)
    m <- net_fit$m
    n <- net_fit$n
    theta_ce <- net_fit$TE.nma.common
    
    Q_total <- net_fit$Q
    phi_hat <- max(1, Q_total / (m - n + 1))
    
    logL_me <- -0.5 * (
      m * log(2 * pi) +
        log(det(phi_hat * V)) +
        t(theta - theta_ce) %*% solve(phi_hat * V) %*% (theta - theta_ce)
    )
    AIC_me[j] <- 2 * n - 2 * logL_me
    
    # Compute AIC for RE model
    tau_hat2 <- net_fit$tau2
    Sigma_re <- V + tau_hat2 * diag(m)
    theta_re <- net_fit$TE.nma.random
    
    logL_re <- -0.5 * (
      m * log(2 * pi) +
        log(det(Sigma_re)) +
        t(theta - theta_re) %*% solve(Sigma_re) %*% (theta - theta_re)
    )
    AIC_re[j] <- 2 * n - 2 * logL_re
    
    # Compute AIC for CE model
    logL_ce <- -0.5 * (
      m * log(2 * pi) +
        log(det(V)) +
        t(theta - theta_ce) %*% solve(V) %*% (theta - theta_ce)
    )
    AIC_ce[j] <- 2 * n - 2 * logL_ce
    
    
    taus[j] <- tau_hat2
    phis[j] <- phi_hat
    
    Q[j] <- net_fit$Q
    Q_pval[j] <- net_fit$pval.Q
    
    Qh[j] <- net_fit$Q.heterogeneity
    Qh_pval[j] <- net_fit$pval.Q.heterogeneity
    
    Qi[j] <- net_fit$Q.inconsistency
    Qi_pval[j] <- net_fit$pval.Q.inconsistency
    
    I2[j] <- net_fit$I2
    I2_lower[j] <- net_fit$lower.I2
    I2_upper[j] <- net_fit$upper.I2
  }
  
  return(list(
    index = idx,
    recid = idx_use$recid,
    effect_measure = idx_use$effect_measure,
    tau2 = taus,
    phi = phis,
    AIC_re = AIC_re,
    AIC_ce = AIC_ce,
    AIC_me = AIC_me,
    me_re = AIC_me - AIC_re,
    me_ce = AIC_me - AIC_ce,
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