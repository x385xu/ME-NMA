#' @description
#' Summarize quantiles of \eqn{\Delta}AIC by effect measure
#'
#' This function computes descriptive summaries of \eqn{\Delta}AIC (\code{me_re}) 
#' within each effect measure category. The effect measures are ordered as odds 
#' ratio, risk ratio, and mean difference. The returned table contains the minimum,
#' first quantile, median, third quantile, maximum, and sample size.
#'
#' @param dat Data frame. Input dataset containing at least \code{effect_measure}
#' and \code{me_re}.
#'
#' @return A tibble with one row per effect measure and the following columns:
#' \itemize{
#'   \item \code{effect_measure}: Effect measure category.
#'   \item \code{Min}: Minimum value of \code{me_re}.
#'   \item \code{q25}: First quantile of \code{me_re}.
#'   \item \code{Median}: Median of \code{me_re}.
#'   \item \code{q75}: Third quantile of \code{me_re}.
#'   \item \code{Max}: Maximum value of \code{me_re}.
#'   \item \code{n}: Number of observations.
#' }
#'
make_quantile_table <- function(dat) {
  dat %>%
    mutate(
      effect_measure = factor(
        effect_measure,
        levels = c("odds ratio", "risk ratio", "mean difference")
      )
    ) %>%
    group_by(effect_measure) %>%
    summarise(
      Min = min(me_re, na.rm = TRUE),
      q25 = quantile(me_re, 0.25, na.rm = TRUE),
      Median = median(me_re, na.rm = TRUE),
      q75 = quantile(me_re, 0.75, na.rm = TRUE),
      Max = max(me_re, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}