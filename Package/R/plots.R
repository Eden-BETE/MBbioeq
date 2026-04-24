# =============================================================================
# Plotting functions for MBbioeq
# =============================================================================

#' Spaghetti plot of individual PK concentration-time profiles
#'
#' @description
#' Plots individual concentration-time profiles coloured by treatment group
#' (Test or Reference), with an optional overlay of the arithmetic mean
#' profile per group.  Useful for checking data quality, identifying outliers,
#' and visually comparing the two formulations before model fitting.
#'
#' @param data_path Character. Path to the dataset file (same format as
#'   accepted by [run_MBTost()] and [run_MBbot()]; separator is detected
#'   automatically).
#' @param title     Character. Plot title.
#' @param col_T     Character. Colour for Test (T) lines and mean profile
#'   (default `"#E63946"`, red).
#' @param col_R     Character. Colour for Reference (R) lines and mean profile
#'   (default `"#457B9D"`, blue).
#' @param add_mean  Logical. If `TRUE` (default), overlay the group mean
#'   profile as a thicker line.
#' @param alpha_ind Numeric in (0, 1]. Transparency of individual lines
#'   (`1` = fully opaque, default `0.35`).
#' @param xlab      Character. X-axis label.
#' @param ylab      Character. Y-axis label.
#'
#' @return Invisibly returns the pre-processed `data.frame` (output of
#'   [prep_parallel_TR()], NA concentrations removed).
#'
#' @seealso [prep_parallel_TR()], [run_MBTost()], [run_MBbot()]
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics lines legend plot
#' @importFrom utils read.table
#'
#' @export
#' @examples
#' \dontrun{
#' dp <- system.file("extdata", "dataset_1.txt", package = "MBbioeq")
#' plot_spaghetti(dp)
#' # Custom colours:
#' plot_spaghetti(dp, col_T = "darkorange", col_R = "steelblue")
#' }
plot_spaghetti <- function(data_path,
                           title     = "Individual PK Profiles",
                           col_T     = "#E63946",
                           col_R     = "#457B9D",
                           add_mean  = TRUE,
                           alpha_ind = 0.35,
                           xlab      = "Time (hr)",
                           ylab      = "Concentration (mg/L)") {

  tab <- .read_pk_data(data_path)
  tab <- prep_parallel_TR(tab)
  tab <- tab[!is.na(tab$Concentration), ]

  ids       <- unique(tab$Id)
  treat_map <- setNames(tab$Treat[match(ids, tab$Id)], as.character(ids))

  col_T_tr <- grDevices::adjustcolor(col_T, alpha.f = alpha_ind)
  col_R_tr <- grDevices::adjustcolor(col_R, alpha.f = alpha_ind)

  xlim <- range(tab$Time, na.rm = TRUE)
  ylim <- c(0, max(tab$Concentration, na.rm = TRUE) * 1.10)

  graphics::plot(NA,
                 xlim = xlim, ylim = ylim,
                 xlab = xlab,  ylab = ylab,
                 main = title, las = 1, bty = "l")

  # Individual lines -------------------------------------------------------
  for (id in ids) {
    sub      <- tab[tab$Id == id, ]
    sub      <- sub[order(sub$Time), ]
    col_line <- if (treat_map[as.character(id)] == 1L) col_T_tr else col_R_tr
    graphics::lines(sub$Time, sub$Concentration, col = col_line, lwd = 1)
  }

  # Group mean profiles ----------------------------------------------------
  if (add_mean) {
    for (tr in c(0L, 1L)) {
      sub_tr <- tab[tab$Treat == tr, ]
      if (nrow(sub_tr) == 0L) next
      time_pts  <- sort(unique(sub_tr$Time))
      mean_conc <- vapply(time_pts, function(t) {
        mean(sub_tr$Concentration[sub_tr$Time == t], na.rm = TRUE)
      }, numeric(1L))
      col_mean <- if (tr == 1L) col_T else col_R
      graphics::lines(time_pts, mean_conc, col = col_mean, lwd = 2.5)
    }
  }

  graphics::legend(
    "topright",
    legend = c("Test (T) — mean",       "Reference (R) — mean",
               "Test (T) — individual", "Reference (R) — individual"),
    col    = c(col_T, col_R, col_T_tr, col_R_tr),
    lty = 1L, lwd = c(2.5, 2.5, 1, 1),
    bty = "n", cex = 0.82
  )

  invisible(tab)
}


#' Forest plot of GMR and 90 percent confidence intervals from TOST analysis
#'
#' @description
#' Displays a horizontal forest plot of the geometric mean ratio (GMR) and its
#' 90 % confidence interval(s) for each available endpoint (AUC, Cmax).  The
#' bioequivalence acceptance zone \[0.80, 1.25\] is highlighted in green.
#' Normal and/or Gallant-corrected confidence intervals can be shown.
#'
#' @param result  A list returned by [run_MBTost()].
#' @param method  Character. Which inference to display: `"both"` (default),
#'   `"normal"`, or `"gallant"`.
#' @param delta   Numeric. Equivalence margin on the log scale (default
#'   `log(1.25)`); determines the acceptance zone drawn on the plot.
#' @param col_inside  Character. Colour used when the CI lies entirely within
#'   the acceptance zone (default dark green `"#2D6A4F"`).
#' @param col_outside Character. Colour used when the CI extends beyond a
#'   boundary (default red `"#D62828"`).
#' @param col_zone    Character. Fill colour of the acceptance zone rectangle
#'   (default light green `"#B7E4C7"`).
#'
#' @return Invisibly returns the `data.frame` used for plotting, with columns
#'   `metric`, `method`, `gmr`, `lower`, `upper`, and `be`.
#'
#' @seealso [run_MBTost()] to obtain the `result` argument;
#'   [plot_bot()] for the BOT equivalent.
#'
#' @importFrom graphics abline legend lines mtext par plot points rect text
#'
#' @export
#' @examples
#' \dontrun{
#' dp  <- system.file("extdata", "dataset_1.txt", package = "MBbioeq")
#' res <- run_MBTost(dp, model = "1cpt")
#' plot_tost(res)
#' plot_tost(res, method = "gallant")   # Gallant CIs only
#' }
plot_tost <- function(result,
                      method      = c("both", "normal", "gallant"),
                      delta       = log(1.25),
                      col_inside  = "#2D6A4F",
                      col_outside = "#D62828",
                      col_zone    = "#B7E4C7") {

  method <- match.arg(method)

  lo_bound <- exp(-delta)   # 0.80
  hi_bound <- exp(delta)    # 1.25

  # Build plotting data frame ----------------------------------------------
  rows <- list()
  add_row <- function(metric, meth, gmr, lower, upper, be) {
    rows[[length(rows) + 1L]] <<- data.frame(
      metric = metric, method = meth,
      gmr = gmr, lower = lower, upper = upper, be = be,
      stringsAsFactors = FALSE
    )
  }

  has_cmax <- !is.na(result$cmax$estimate)

  if (method %in% c("both", "normal")) {
    add_row("AUC",  "Normal",
            result$auc$estimate,
            result$auc$normal$ci_lower, result$auc$normal$ci_upper,
            result$auc$normal$decision)
    if (has_cmax)
      add_row("Cmax", "Normal",
              result$cmax$estimate,
              result$cmax$normal$ci_lower, result$cmax$normal$ci_upper,
              result$cmax$normal$decision)
  }
  if (method %in% c("both", "gallant")) {
    add_row("AUC",  "Gallant",
            result$auc$estimate,
            result$auc$gallant$ci_lower, result$auc$gallant$ci_upper,
            result$auc$gallant$decision)
    if (has_cmax)
      add_row("Cmax", "Gallant",
              result$cmax$estimate,
              result$cmax$gallant$ci_lower, result$cmax$gallant$ci_upper,
              result$cmax$gallant$decision)
  }

  df     <- do.call(rbind, rows)
  n_rows <- nrow(df)

  x_lo <- min(lo_bound * 0.82, min(df$lower, na.rm = TRUE) * 0.92)
  x_hi <- max(hi_bound * 1.12, max(df$upper, na.rm = TRUE) * 1.08)

  old_par <- graphics::par(mar = c(5, 9, 4, 3))
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::plot(NA,
                 xlim = c(x_lo, x_hi),
                 ylim = c(0.3, n_rows + 0.7),
                 xlab = "Geometric Mean Ratio  T / R",
                 ylab = "", yaxt = "n",
                 main = "TOST Bioequivalence \u2014 90% CI",
                 las = 1, bty = "l")

  # Acceptance zone --------------------------------------------------------
  graphics::rect(lo_bound, 0.1, hi_bound, n_rows + 0.9,
                 col = col_zone, border = NA)

  graphics::abline(v = 1,        lty = 2L, col = "gray40", lwd = 1.2)
  graphics::abline(v = lo_bound, lty = 3L, col = "gray50", lwd = 1.0)
  graphics::abline(v = hi_bound, lty = 3L, col = "gray50", lwd = 1.0)

  graphics::mtext(sprintf("%.2f", lo_bound), side = 1, at = lo_bound,
                  line = 2.3, cex = 0.72, col = "gray40")
  graphics::mtext(sprintf("%.2f", hi_bound), side = 1, at = hi_bound,
                  line = 2.3, cex = 0.72, col = "gray40")

  # CI rows ----------------------------------------------------------------
  for (i in seq_len(n_rows)) {
    y   <- n_rows + 1L - i
    be  <- isTRUE(df$be[i])
    col <- if (be) col_inside else col_outside

    graphics::lines(c(df$lower[i], df$upper[i]), c(y, y), col = col, lwd = 2)
    graphics::lines(rep(df$lower[i], 2), c(y - 0.13, y + 0.13), col = col, lwd = 2)
    graphics::lines(rep(df$upper[i], 2), c(y - 0.13, y + 0.13), col = col, lwd = 2)
    graphics::points(df$gmr[i], y, pch = 18L, cex = 1.9, col = col)

    # Row label in the left margin
    graphics::mtext(paste0(df$metric[i], "  (", df$method[i], ")"),
                    side = 2, at = y, las = 1, cex = 0.82, line = 0.5)

    # Numeric annotation above the bar
    graphics::text(df$gmr[i], y + 0.30,
                   labels = sprintf("%.3f  [%.3f, %.3f]",
                                    df$gmr[i], df$lower[i], df$upper[i]),
                   cex = 0.64, col = col)
  }

  graphics::legend(
    "bottomright",
    legend = c(
      sprintf("BE zone [%.2f \u2013 %.2f]", lo_bound, hi_bound),
      "Bioequivalent",
      "Not bioequivalent"
    ),
    fill   = c(col_zone, col_inside, col_outside),
    border = NA, bty = "n", cex = 0.82
  )

  invisible(df)
}


#' Folded-normal decision plot for BOT bioequivalence analysis
#'
#' @description
#' Visualises the Bayesian-Optimal Test (BOT) decision by plotting, for each
#' available endpoint, the folded-normal density under the worst-case null
#' hypothesis alongside the observed test statistic and the critical value.
#' This makes the strength of the bioequivalence evidence immediately visible.
#'
#' @section Statistical background:
#' The BOT test statistic is \eqn{T = |\hat\beta|}.  Under the worst-case null
#' \eqn{H_0: |\beta| = \delta}, \eqn{T} follows a folded-normal distribution
#' \eqn{FN(\delta, \hat\sigma^2)}, whose density is
#' \deqn{f(x) = \phi\!\left(\frac{x - \delta}{\hat\sigma}\right) / \hat\sigma
#'              + \phi\!\left(\frac{x + \delta}{\hat\sigma}\right) / \hat\sigma,
#'              \quad x \ge 0,}
#' where \eqn{\phi} is the standard normal density.  The critical value
#' \eqn{c^* = q_{FN}(\alpha; \delta, \hat\sigma)} is the \eqn{\alpha}-quantile
#' of that distribution.  Bioequivalence is concluded when
#' \eqn{T \le c^*} (the blue line falls in the **green region**).
#'
#' @param result     A list returned by [run_MBbot()].  Must contain `se`
#'   elements inside `auc` and (if available) `cmax`.
#' @param delta      Numeric. Equivalence margin (default `log(1.25)`).
#' @param alpha      Numeric. One-sided significance level (default `0.05`).
#' @param col_accept Character. Fill colour of the acceptance region
#'   \eqn{[0, c^*]} (default light green `"#B7E4C7"`).
#' @param col_reject Character. Fill colour of the rejection region
#'   \eqn{(c^*, \infty)} (default light red `"#FFCCD5"`).
#' @param col_stat   Character. Colour of the observed test-statistic line
#'   (default dark blue `"#023E8A"`).
#' @param col_crit   Character. Colour of the critical-value line
#'   (default red `"#D62828"`).
#'
#' @return Invisibly `NULL`.
#'
#' @seealso [run_MBbot()] to obtain the `result` argument;
#'   [plot_tost()] for the TOST equivalent.
#'
#' @importFrom graphics abline legend lines mtext par plot polygon text
#' @importFrom stats dnorm
#'
#' @export
#' @examples
#' \dontrun{
#' dp  <- system.file("extdata", "dataset_1.txt", package = "MBbioeq")
#' res <- run_MBbot(dp, model = "1cpt")
#' plot_bot(res)
#' }
plot_bot <- function(result,
                     delta      = log(1.25),
                     alpha      = 0.05,
                     col_accept = "#B7E4C7",
                     col_reject = "#FFCCD5",
                     col_stat   = "#023E8A",
                     col_crit   = "#D62828") {

  has_auc  <- !is.na(result$auc$statistic)  && !is.null(result$auc$se)
  has_cmax <- !is.na(result$cmax$statistic) && !is.null(result$cmax$se)

  if (!has_auc && !has_cmax) {
    warning(paste0(
      "No SE values found in the BOT result.\n",
      "Please re-run run_MBbot() with the current package version."
    ))
    return(invisible(NULL))
  }

  metrics  <- c(if (has_auc) "AUC", if (has_cmax) "Cmax")
  n_panels <- length(metrics)

  old_par <- graphics::par(
    mfrow = c(1L, n_panels),
    mar   = c(5, 4.5, 5.5, 2),
    oma   = c(0, 0, 2.5, 0)
  )
  on.exit(graphics::par(old_par))

  for (metric in metrics) {

    dat      <- result[[tolower(metric)]]
    stat_val <- dat$statistic       # |beta_hat|
    crit_val <- dat$critical_value  # qfoldnorm(alpha, delta, se)
    se_val   <- dat$se
    decision <- dat$decision

    # Folded-normal density -----------------------------------------------
    x_max <- max(delta * 3.2, stat_val * 1.5, crit_val * 1.5)
    xv    <- seq(0, x_max, length.out = 600L)
    yv    <- stats::dnorm(xv,  delta, se_val) +
             stats::dnorm(xv, -delta, se_val)

    ylim_top <- max(yv) * 1.28

    graphics::plot(xv, yv, type = "n",
                   xlim = c(0, x_max),
                   ylim = c(0, ylim_top),
                   xlab = expression(group("|", hat(beta), "|")),
                   ylab = "Density",
                   main = paste0(metric, " \u2014 BOT decision"),
                   las = 1, bty = "l")

    # Shade acceptance region [0, c*] -------------------------------------
    idx_acc <- xv <= crit_val
    graphics::polygon(
      c(xv[idx_acc], rev(xv[idx_acc])),
      c(yv[idx_acc], rep(0, sum(idx_acc))),
      col = col_accept, border = NA
    )

    # Shade rejection region (c*, x_max] ----------------------------------
    idx_rej <- xv >= crit_val
    graphics::polygon(
      c(xv[idx_rej], rev(xv[idx_rej])),
      c(yv[idx_rej], rep(0, sum(idx_rej))),
      col = col_reject, border = NA
    )

    # Density curve on top ------------------------------------------------
    graphics::lines(xv, yv, lwd = 2, col = "gray25")

    # Critical value line and label ---------------------------------------
    graphics::abline(v = crit_val, lty = 2L, col = col_crit, lwd = 2)
    pos_crit <- if (crit_val < x_max * 0.75) 4L else 2L
    graphics::text(crit_val, ylim_top * 0.90,
                   labels = bquote(c^"*" == .(sprintf("%.3f", crit_val))),
                   col = col_crit, cex = 0.82, pos = pos_crit)

    # Observed statistic line and label -----------------------------------
    graphics::abline(v = stat_val, lty = 1L, col = col_stat, lwd = 2.5)
    pos_stat <- if (stat_val < crit_val) 4L else 2L
    graphics::text(stat_val, ylim_top * 0.75,
                   labels = bquote(group("|", hat(beta), "|") ==
                                     .(sprintf("%.3f", stat_val))),
                   col = col_stat, cex = 0.82, pos = pos_stat)

    # Decision banner at top of panel -------------------------------------
    dec_txt <- if (isTRUE(decision)) "BIOEQUIVALENT" else "NOT BIOEQUIVALENT"
    dec_col <- if (isTRUE(decision)) "#1B4332"       else "#9B2226"
    graphics::mtext(dec_txt, side = 3, line = 1.2,
                    col = dec_col, font = 2, cex = 1.05)

    graphics::mtext(
      bquote(alpha == .(alpha) ~~ delta == .(sprintf("%.3f", delta)) ~~
               "(" * e^delta == .(sprintf("%.2f", exp(delta))) * ")"),
      side = 3, line = 0.1, cex = 0.72, col = "gray45"
    )

    # Legend --------------------------------------------------------------
    graphics::legend(
      x = x_max * 0.45, y = ylim_top * 0.60,
      legend = c(
        bquote("Accept  [0," ~ c^"*" * "]"),
        bquote("Reject  (" * c^"*" * ", " * infinity * ")")
      ),
      fill   = c(col_accept, col_reject),
      border = "gray70",
      bty    = "n", cex = 0.78
    )
  }

  graphics::mtext("BOT \u2014 Folded-Normal Decision Plot",
                  outer = TRUE, cex = 1.15, font = 2)
  invisible(NULL)
}
