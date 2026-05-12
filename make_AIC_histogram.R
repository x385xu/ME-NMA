#' Create stratified histograms of \eqn{\Delta AIC} values for NMA models
#'
#' This function generates side-by-side histograms comparing the distribution
#' of \eqn{\Delta AIC = AIC_{ME} - AIC_{RE}} across network meta-analyses (NMAs)
#' for different effect measures (OR, RR, MD). Histograms are additionally
#' stratified by heterogeneity test significance (\eqn{p < 0.1} vs
#' \eqn{p \ge 0.1}).
#'
#' Negative \eqn{\Delta AIC} values favour the multiplicative-effects (ME)
#' model, while positive values favour the additive random-effects (RE) model.
#' Dashed vertical reference lines are drawn at \eqn{\Delta AIC = \pm 3}.
#'
#' For mean difference (MD) networks, the far left tail is optionally compressed
#' to improve visualization of extreme negative values.
#'
#' @param AIC_out A data frame or tibble containing model comparison results.
#'   Required columns include:
#'   \describe{
#'     \item{effect_measure}{Effect measure name (e.g. `"odds ratio"`,
#'       `"risk ratio"`, `"mean difference"`).}
#'     \item{me_re}{Difference in AIC between ME and RE models
#'       (\eqn{AIC_{ME} - AIC_{RE}}).}
#'     \item{Q_pval}{P-value for the heterogeneity test.}
#'   }
#'
#' @return A `patchwork`/`ggplot2` object containing three aligned histograms
#'   for OR, RR, and MD networks.
#'

make_AIC_histogram <- function(AIC_out) {
  
  # ---------------------------------------------------------------------------
  # Prepare plotting dataset
  # ---------------------------------------------------------------------------
  # - Convert effect measure labels to short forms (OR/RR/MD)
  # - Extract ΔAIC values (ME - RE)
  # - Create grouping variable based on heterogeneity p-value
  # - Remove missing values needed for plotting
  df_AICplot <- tibble::as_tibble(AIC_out) %>%
    dplyr::mutate(
      measure = dplyr::case_when(
        effect_measure == "odds ratio" ~ "OR",
        effect_measure == "risk ratio" ~ "RR",
        effect_measure == "mean difference" ~ "MD",
        TRUE ~ effect_measure
      ),
      measure = factor(measure, levels = c("OR", "RR", "MD")),
      value = me_re,
      pval = Q_pval,
      hetero_group = dplyr::case_when(
        is.na(pval) ~ NA_character_,
        pval < 0.10 ~ "p_less_0.1",
        TRUE ~ "p_ge_0.1"
      ),
      hetero_group = factor(hetero_group, levels = c("p_less_0.1", "p_ge_0.1"))
    ) %>%
    dplyr::filter(!is.na(value), !is.na(hetero_group))
  
  # Maximum y-lab
  common_ymax <- 22
  
  # Helper function to create one histogram panel
  make_one <- function(dat, title,
                       compress_md = FALSE,
                       show_legend = TRUE,
                       ymax_common = common_ymax) {
    # Parameters used to compress the extreme left MD tail
    md_break_left <- -8
    md_tail_target <- -7
    # Find minimum x-position among extreme left-tail bins
    tail_min <- dat %>%
      dplyr::mutate(
        bin_left = floor(value),
        x_raw = bin_left + 0.5,
        x_plot_tmp = dplyr::if_else(value < 0, x_raw - 0.5, x_raw + 0.5)
      ) %>%
      dplyr::filter(x_plot_tmp < md_break_left) %>%
      dplyr::summarise(tail_min = min(x_plot_tmp, na.rm = TRUE)) %>%
      dplyr::pull(tail_min)
    # Amount by which extreme left bins are shifted
    shift_md <- if (compress_md && is.finite(tail_min)) {
      md_tail_target - tail_min
    } else {
      0
    }
    # Construct histogram bins
    bins <- dat %>%
      dplyr::mutate(
        bin_left = floor(value),
        bin_type = dplyr::case_when(
          value == 0 ~ "zero",
          value < 0  ~ "neg",
          value > 0  ~ "pos"
        ),
        bin_label = dplyr::case_when(
          value == 0 ~ "0",
          value < -1  ~ paste0("(", bin_left, ", ", bin_left + 1, "]"),
          value < 0  ~ paste0("(", bin_left, ", ", bin_left + 1, ")"),
          value > 0  ~ paste0("(", bin_left, ", ", bin_left + 1, "]")
        ),
        x_raw = dplyr::case_when(
          value == 0 ~ 0,
          value < 0  ~ bin_left + 0.5,
          value > 0  ~ bin_left + 0.5
        ),
        # Shift negative and positive bins away from zero
        x_plot = dplyr::case_when(
          bin_type == "neg"  ~ x_raw - 0.5,
          bin_type == "zero" ~ 0,
          bin_type == "pos"  ~ x_raw + 0.5
        ),
        # Compress extreme left MD tail if requested
        x_plot = if (compress_md) {
          dplyr::if_else(x_plot < md_break_left, x_plot + shift_md, x_plot)
        } else {
          x_plot
        }
      ) %>%
      # Count networks within each bin
      dplyr::count(x_plot, bin_label, hetero_group)
    
    # Vertical reference lines at x = ±3
    vlines <- tibble::tibble(
      x = c(-3.5, 3.5)
    ) %>%
      dplyr::mutate(
        x = if (compress_md) {
          dplyr::if_else(x < -8, x + 23, x)
        } else {
          x
        }
      )
    
    # Build histogram plot
    ggplot(bins, aes(x = x_plot, y = n, fill = hetero_group)) +
      # Red dashed reference lines at ±3
      geom_segment(
        data = vlines,
        aes(x = x, xend = x, y = 0, yend = ymax_common),
        inherit.aes = FALSE,
        linetype = "dashed",
        color = "red",
        linewidth = 0.7
      ) +
      # Grey vertical line at 0
      geom_segment(
        aes(x = 0, xend = 0, y = 0, yend = ymax_common),
        inherit.aes = FALSE,
        linewidth = 0.5,
        colour = "grey75"
      ) +
      # Histogram bars
      geom_col(
        width = 0.95,
        color = "black",
        alpha = 0.8,
        position = "stack",
        show.legend = show_legend
      ) +
      # Rotated interval labels below bars
      geom_text(
        data = dplyr::distinct(bins, x_plot, bin_label),
        aes(x = x_plot, y = -0.08 * ymax_common, label = bin_label),
        inherit.aes = FALSE,
        show.legend = FALSE,
        angle = 90,
        hjust = 0.7,
        size = 6
      ) +
      scale_fill_manual(
        values = c(
          "p_less_0.1" = "grey80",
          "p_ge_0.1" = "grey20"
        ),
        labels = c(
          expression(p < 0.1),
          expression(p >= 0.1)
        )
      ) +
      # Shared y-axis formatting
      scale_y_continuous(
        breaks = seq(0, ymax_common, by = 2),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(
        ylim = c(0, ymax_common),
        clip = "off"
      ) +
      labs(title = title, x = "", y = "", fill = "") +
      # Add visual break marker for compressed MD tail
      {
        if (compress_md) {
          list(
            annotate(
              "segment",
              x = md_tail_target + 1.2, xend = md_tail_target + 1.45,
              y = -0.35, yend = 0.2,
              linewidth = 0.8,
              colour = "grey40"
            ),
            annotate(
              "segment",
              x = md_tail_target + 1.45, xend = md_tail_target + 1.7,
              y = -0.35, yend = 0.2,
              linewidth = 0.8,
              colour = "grey40"
            )
          )
        }
      } +
      theme_minimal(base_size = 16) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),,
        panel.grid.major.y = element_line(linewidth = 0.5, colour = "grey75"),
        panel.grid.minor = element_blank(),
        plot.margin = if (title == "MD") {
          margin(20, 10, 80, 10)
        } else {
          margin(5.5, 10, 25, 10)
        }
      )
  }
  
  # Create separate panels for OR, RR, and MD
  p_or <- make_one(dplyr::filter(df_AICplot, measure == "OR"), "OR", show_legend = TRUE)
  p_rr <- make_one(dplyr::filter(df_AICplot, measure == "RR"), "RR", show_legend = FALSE)
  p_md <- make_one(dplyr::filter(df_AICplot, measure == "MD"), "MD", compress_md = TRUE, show_legend = FALSE)
  
  final_plot <-
    p_or + p_rr + p_md +
    patchwork::plot_layout(
      nrow = 1,
      guides = "collect",
      widths = c(1.5, 1, 0.85)
    ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.margin = margin(t = 25),
      legend.box.margin = margin(t = 25),
      
      legend.text = element_text(size = 18),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.spacing.x = grid::unit(0.8, "cm")
    )

  
  return(final_plot)
}