#' @description
#' Construct case study outputs for network meta-analysis
#' This function processes a network meta-analysis dataset (from \texttt{nmadb})
#' and generates study-level summaries, heterogeneity decomposition, and
#' visualization tools for case study analysis. It produces a netgraph and a
#' customized forest plot comparing treatments against a reference (treatment 1),
#' alongside a table of study-level contributions to heterogeneity.
#'
#' @param recid Integer. Record ID used by \texttt{nmadb} to identify the dataset.
#'
#' @param treatment Character vector. Ordered list of treatment labels. The first
#' element is treated as the reference treatment in all comparisons and determines
#' the orientation of treatment effects in the forest plot.
#'
#' @param xlab Character (optional). Label for the x-axis in the forest plot.
#' If \code{NULL}, defaults to \code{"Log-OR (vs. <reference treatment>)"}.
#'
#' @param refline Numeric. Location of the vertical reference line in the forest
#' plot (typically 0 for log effect measures, corresponding to no effect).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{net}: Fitted \texttt{netmeta} object.
#'   \item \code{dat}: Raw dataset retrieved via \texttt{readByID()}.
#'   \item \code{df_table}: Data frame of study-level estimates and heterogeneity
#'   contributions (\eqn{Q_{\text{het}}^i}), suitable for reporting.
#'   \item \code{df_plot}: Processed data used to construct the forest plot.
#'   \item \code{netgraph}: Function that generates the network graph when called.
#'   \item \code{forest}: Function that generates the customized forest plot when called.
#' }
#'
#' @note
#' The plotting functions (\code{netgraph} and \code{forest}) are returned as
#' callable closures to allow flexible control of graphical devices (e.g.,
#' \texttt{png()}, \texttt{pdf()}) outside the function.

analyse_case_study <- function(recid,
                               treatment,
                               xlab = NULL,
                               refline = 0) {
  
  recid <- as.character(recid)
  dat <- twoarm_data_list[[recid]]
  effect <- dat$effect
  
  net <- fit_netmeta(
    indata = dat,
    model = "random"
  )
  
  ## -----------------------------
  ## Netgraph function
  ## -----------------------------
  
  plot_netgraph <- function() {
    netgraph(
      net,
      labels = treatment,
      plastic = FALSE,
      iterate = FALSE,
      number.of.studies = TRUE,
      col = "grey",
      cex = 1,
      points = TRUE,
      col.points = "orange",
      cex.points = 3,
      cex.number.of.studies = 1,
      pos.number.of.studies = 0.4,
      lwd = 2,
      offset = 0.06
    )
  }
  
  ## -----------------------------
  ## Study-level estimates and Q decomposition
  ## -----------------------------
  
  ref_treatment <- treatment[1]
  
  df_study <- tibble(
    study  = net$studlab,
    treat1 = treatment[as.numeric(net$treat1)],
    treat2 = treatment[as.numeric(net$treat2)],
    t1_num = as.numeric(net$treat1),
    t2_num = as.numeric(net$treat2),
    TE     = -net$TE,
    seTE   = net$seTE
  ) %>%
    mutate(
      subgroup   = treat2,
      comparison = paste(treat2, "vs", treat1)
    ) %>%
    group_by(comparison) %>%
    mutate(
      theta_pair_FE = sum(TE / seTE^2) / sum(1 / seTE^2),
      Q_het_i = (TE - theta_pair_FE)^2 / seTE^2,
      studlab2 = sprintf("Q_het^i = %.2f", Q_het_i)
    ) %>%
    ungroup()
  
  df_table <- df_study %>%
    select(treat1, treat2, Q_het_i, TE, seTE)
  
  ## -----------------------------
  ## ME estimates
  ## -----------------------------
  
  m <- net$m
  n <- net$n
  
  theta    <- net$TE
  theta_me <- net$TE.nma.common
  V        <- diag(net$seTE^2)
  
  phi <- as.numeric(
    t(theta - theta_me) %*% solve(V) %*% (theta - theta_me) /
      (m - n + 1)
  )
  
  d_me  <- net$TE.common[-1, 1]
  se_me <- sqrt(phi) * net$seTE.common[-1, 1]
  
  ## -----------------------------
  ## RE estimates
  ## -----------------------------
  
  d_re  <- net$TE.random[-1, 1]
  se_re <- net$seTE.random[-1, 1]
  
  ## -----------------------------
  ## Plot data
  ## -----------------------------
  df_study_plot <- df_study %>%
    filter(t1_num == 1 | t2_num == 1) %>%
    mutate(
      active_treatment = ifelse(t1_num == 1, treat2, treat1),
      TE = ifelse(t1_num == 1, TE, -TE),  # orient vs reference
      subgroup = active_treatment,
      comparison = paste(active_treatment, "vs", ref_treatment),
      treat1 = ref_treatment,
      treat2 = active_treatment
    )
  
  direct_treatments <- df_study_plot %>%
    distinct(subgroup) %>%
    pull(subgroup)
  
  direct_idx <- match(direct_treatments, treatment)
  
  df_summary <- tibble(
    study = rep(c("RE estimate", "ME estimate"), each = length(direct_treatments)),
    subgroup = rep(direct_treatments, times = 2),
    TE = c(d_re[direct_idx - 1], d_me[direct_idx - 1]),
    seTE = c(se_re[direct_idx - 1], se_me[direct_idx - 1]),
    Q_het_i = NA_real_,
    studlab2 = study
  )
  
  df_plot <- bind_rows(df_summary, df_study_plot) %>%
    mutate(
      row_type = case_when(
        study == "RE estimate" ~ "RE",
        study == "ME estimate" ~ "ME",
        TRUE ~ "Study"
      ),
      row_order = case_when(
        row_type == "RE"    ~ 1,
        row_type == "ME"    ~ 2,
        row_type == "Study" ~ 3
      ),
      pch = case_when(
        row_type == "Study" ~ 16,
        row_type == "RE"    ~ 15,
        row_type == "ME"    ~ 17
      ),
      col_point = case_when(
        row_type == "Study" ~ "black",
        row_type == "RE"    ~ "#1b9e77",
        row_type == "ME"    ~ "#d95f02"
      ),
      col_ci = case_when(
        row_type == "Study" ~ "grey50",
        row_type == "RE"    ~ "#1b9e77",
        row_type == "ME"    ~ "#d95f02"
      ),
      ci.lb = TE - qnorm(0.975) * seTE,
      ci.ub = TE + qnorm(0.975) * seTE,
      label_right = sprintf("%.2f [%.2f, %.2f]", TE, ci.lb, ci.ub)
    ) %>%
    arrange(subgroup, row_order, study)
  
  df_plot <- df_plot %>%
    group_by(subgroup) %>%
    mutate(row_in_group = row_number()) %>%
    ungroup()
  
  group_sizes <- df_plot %>%
    count(subgroup, name = "n_rows") %>%
    mutate(
      gap = 1,
      start = cumsum(lag(n_rows + gap, default = 0))
    )
  
  df_plot <- df_plot %>%
    left_join(group_sizes, by = "subgroup") %>%
    mutate(rows = max(start + row_in_group) - (start + row_in_group) + 1)
  
  x_min <- min(df_plot$ci.lb, na.rm = TRUE)
  x_max <- max(df_plot$ci.ub, na.rm = TRUE)
  
  x_range <- x_max - x_min
  
  alim <- c(
    floor(x_min),
    ceiling(x_max)
  )
  
  left_pad  <- 0.35 * x_range
  right_pad <- 0.75 * x_range
  
  xlim <- c(
    x_min - left_pad,
    x_max + right_pad
  )
  
  ilab.xpos <- x_max + 0.35 * x_range
  
  left_text_x <- xlim[1] + 0.05 * x_range
  subgroup_x  <- xlim[1] + 0.02 * x_range
  
  subgroup_y <- df_plot %>%
    group_by(subgroup) %>%
    summarise(y = max(rows) + 1, .groups = "drop")
  
  ref_treatment <- treatment[1]
  if (is.null(xlab)) {
    xlab <- paste0("Log-",effect, " (vs. ", ref_treatment, ")")
  }
  
  ## -----------------------------
  ## Forest plot function
  ## -----------------------------
  
  plot_forest <- function() {
    
    op <- par(mar = c(4, 10, 2, 12))
    on.exit(par(op))
    
    forest(
      x = df_plot$TE,
      ci.lb = df_plot$ci.lb,
      ci.ub = df_plot$ci.ub,
      slab = rep("", nrow(df_plot)),
      rows = df_plot$rows,
      pch = df_plot$pch,
      col = df_plot$col_ci,
      colout = df_plot$col_point,
      refline = refline,
      xlab = xlab,
      xlim = xlim,
      alim = alim,
      at = pretty(c(alim[1], alim[2], 0)),
      ilab = df_plot$label_right,
      ilab.xpos = ilab.xpos,
      annotate = FALSE,
      header = FALSE,
      cex = 0.9,
      top = 2
    )
    
    # erase default top horizontal rule
    segments(
      x0 = xlim[1],
      x1 = xlim[2],
      y0 = max(df_plot$rows) + 1,
      y1 = max(df_plot$rows) + 1,
      col = "white",
      lwd = 3
    )
    
    y_header <- max(df_plot$rows) + 1.5
    
    labels_left <- ifelse(
      df_plot$row_type == "Study",
      sprintf("Q[het]^i == %.2f", df_plot$Q_het_i),
      paste0('"', df_plot$studlab2, '"')
    )
    
    text(
      x = ilab.xpos,
      y = y_header,
      labels = "Estimate [95% CI]",
      font = 2,
      cex = 0.9
    )
    
    text(
      x = subgroup_x,
      y = subgroup_y$y,
      labels = subgroup_y$subgroup,
      pos = 4,
      font = 2,
      cex = 0.95
    )
    
    text(
      x = left_text_x,
      y = df_plot$rows,
      labels = parse(text = labels_left),
      pos = 4,
      cex = ifelse(df_plot$row_type == "Study", 0.8, 0.9)
    )
  }
  
  return(list(
    net = net,
    dat = dat,
    df_table = df_table,
    df_plot = df_plot,
    netgraph = plot_netgraph,
    forest = plot_forest
  ))
}