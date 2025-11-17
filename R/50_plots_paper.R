paper_colors <- list(
  proposed = "#D55E00",
  scm = "#0072B2",
  obs = "#000000",
  ribbon = "grey70"
)

theme_scspill_paper <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "grey88",
        linewidth = 0.25
      ),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

.make_year_breaks <- function(times, every = 2L) {
  times <- as.numeric(times)
  if (!length(times)) {
    return(list(breaks = seq_along(times), labels = times))
  }
  idx <- seq_along(times)
  keep <- if (every <= 1) {
    idx
  } else {
    which(((times - min(times, na.rm = TRUE)) %% every) == 0)
  }
  keep <- unique(c(1L, keep, length(times))) # 端は必ず残す
  list(breaks = idx[keep], labels = times[keep])
}

# ---------------------------------------------------------
# Panel (a) : Observed vs Counterfactual (Proposed & SCM)
# Plots observed, SCM counterfactual, and Proposed (SCSPILL)
# counterfactual with credible intervals.
# Includes logic to skip a year and connect the gap with a
# different linetype.
# ---------------------------------------------------------
plot_panel_outcomes <- function(
  fit,
  data,
  treated_unit,
  treatment_dummy,
  y,
  unit_col,
  time_col,
  cred = 0.95,
  skip_year = NULL, # Year (e.g., 2011) to skip (plots as a gap)
  gap_linetype = "dotted", # Linetype for the skipped gap
  title = "(a) Observed Outcomes and Synthetic Control Outcomes",
  x_title = "Year",
  y_title = "Per-capita Cigarette Sales (in packs)",
  year_every = 2,
  force_perfect_pre = TRUE, # Force SCM to match pre-period observed
  save_path = NULL
) {
  # 1. Calculate counterfactual for the Proposed method (SCSPILL)
  cf <- scspill_counterfactual(fit, cred = cred, time_col = time_col)
  T0_idx <- max(cf$t_idx[cf$period == "pre"])

  # 2. Calculate counterfactual for the standard SCM method
  scm <- scm_counterfactual_light(
    data,
    treated_unit,
    treatment_dummy,
    y,
    unit_col,
    time_col
  )
  scm <- merge(
    data.frame(time = cf$time, .idx = cf$.idx),
    scm[, c("time", "y_cf", "period")],
    by = "time",
    all.x = TRUE
  )
  if (force_perfect_pre) {
    idx_pre <- cf$period == "pre"
    scm$y_cf[idx_pre] <- cf$y_obs[idx_pre]
  }

  # 3. Create axis labels
  x_axis_map <- data.frame(idx = cf$.idx, time = cf$time)
  yrs <- x_axis_map$time
  idx <- x_axis_map$idx

  if (is.numeric(yrs)) {
    keep_idx <- which(
      (as.numeric(yrs) - min(as.numeric(yrs), na.rm = TRUE)) %% year_every == 0
    )
    keep_idx <- unique(c(1, keep_idx, length(yrs)))
    x_breaks <- idx[keep_idx]
    x_labels <- yrs[keep_idx]
  } else {
    x_breaks <- idx
    x_labels <- yrs
  }

  # 4. Prepare data.frames for plotting
  obs_line <- data.frame(
    idx = cf$.idx,
    y = cf$y_obs,
    series = "Observed",
    time = cf$time
  )
  pre_last_cf <- subset(cf, period == "pre")
  pre_last_cf <- pre_last_cf[nrow(pre_last_cf), ]

  rib_post <- subset(cf, period == "post")
  rib_post <- rib_post[
    is.finite(rib_post$y_cf_lo) & is.finite(rib_post$y_cf_hi),
  ]
  rib_post_continuous <- rbind(
    data.frame(
      .idx = pre_last_cf$.idx,
      y_cf_lo = pre_last_cf$y_cf_mean,
      y_cf_hi = pre_last_cf$y_cf_mean,
      time = pre_last_cf$time
    ),
    rib_post[, c(".idx", "y_cf_lo", "y_cf_hi", "time")]
  )

  prop_pre <- data.frame(
    idx = cf$.idx[cf$period == "pre"],
    y = cf$y_cf_mean[cf$period == "pre"],
    time = cf$time[cf$period == "pre"]
  )
  prop_post_line <- rbind(
    data.frame(
      idx = pre_last_cf$.idx,
      y = pre_last_cf$y_cf_mean,
      time = pre_last_cf$time
    ),
    data.frame(
      idx = cf$.idx[cf$period == "post"],
      y = cf$y_cf_mean[cf$period == "post"],
      time = cf$time[cf$period == "post"]
    )
  )

  scm_line <- data.frame(
    idx = scm$.idx,
    y = scm$y_cf,
    series = "SCM",
    time = scm$time
  )
  scm_line <- scm_line[is.finite(scm_line$y), ]

  # 5. Handle 'skip_year': Create separate data.frames for gap
  obs_gap <- data.frame()
  prop_gap <- data.frame()
  scm_gap <- data.frame()

  if (!is.null(skip_year)) {
    time_skip <- skip_year
    time_bridge <- c(time_skip - 1, time_skip + 1) # Assumes annual data

    idx_skip <- cf$.idx[which(cf$time == time_skip)]
    idx_bridge <- cf$.idx[which(cf$time %in% time_bridge)]

    if (length(idx_skip) > 0 && length(idx_bridge) == 2) {
      # 5a. Create data.frames for the gap lines (2 points each)
      obs_gap <- obs_line[obs_line$idx %in% idx_bridge, ]
      prop_gap <- rbind(prop_pre, prop_post_line)
      prop_gap <- prop_gap[prop_gap$idx %in% idx_bridge, ]
      scm_gap <- scm_line[scm_line$idx %in% idx_bridge, ]

      # 5b. Set the 'y' value of the skipped year to NA in the *main* data
      # This forces the main geom_line to break
      obs_line$y[obs_line$idx == idx_skip] <- NA
      prop_pre$y[prop_pre$idx == idx_skip] <- NA
      prop_post_line$y[prop_post_line$idx == idx_skip] <- NA
      scm_line$y[scm_line$idx == idx_skip] <- NA
      rib_post_continuous$y_cf_lo[rib_post_continuous$.idx == idx_skip] <- NA
      rib_post_continuous$y_cf_hi[rib_post_continuous$.idx == idx_skip] <- NA
    }
  }

  # 6. Construct the ggplot object
  gg <- ggplot() +
    geom_line(
      data = obs_line,
      aes(idx, y, color = "Observed", linetype = "Observed"),
      linewidth = 0.5
    ) +
    geom_ribbon(
      data = rib_post_continuous,
      aes(
        x = .idx,
        ymin = y_cf_lo,
        ymax = y_cf_hi,
        fill = "95% Credible Interval"
      ),
      alpha = 0.18
    ) +
    geom_line(
      data = prop_pre,
      aes(idx, y, color = "Proposed", linetype = "Proposed"),
      linewidth = 0.6
    ) +
    geom_line(
      data = prop_post_line,
      aes(idx, y, color = "Proposed", linetype = "Proposed"),
      linewidth = 0.6
    ) +
    geom_line(
      data = scm_line,
      aes(idx, y, color = "SCM", linetype = "SCM"),
      linewidth = 0.6
    ) +
    # Intervention line
    geom_vline(xintercept = T0_idx, linetype = "dotted", color = "grey20") +

    # Scales, theme, and legends
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Observed" = paper_colors$obs,
        "Proposed" = paper_colors$proposed,
        "SCM" = paper_colors$scm
      ),
      breaks = c("Observed", "Proposed", "SCM")
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Observed" = "solid",
        "Proposed" = "dashed",
        "SCM" = "dotted"
      ),
      breaks = c("Observed", "Proposed", "SCM")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("95% Credible Interval" = paper_colors$ribbon)
    ) +
    labs(x = x_title, y = y_title, title = title) +
    theme_scspill_paper() +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    )

  interpolation_layers <- list(
    if (nrow(obs_gap) > 0) {
      geom_line(
        data = obs_gap,
        aes(idx, y, color = "Observed"),
        linetype = gap_linetype,
        linewidth = 0.5
      )
    },

    if (nrow(prop_gap) > 0) {
      geom_line(
        data = prop_gap,
        aes(idx, y, color = "Proposed"),
        linetype = gap_linetype,
        linewidth = 0.6
      )
    },

    if (nrow(scm_gap) > 0) {
      geom_line(
        data = scm_gap,
        aes(idx, y, color = "SCM"),
        linetype = gap_linetype,
        linewidth = 0.6
      )
    }
  )

  gg <- gg + interpolation_layers

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path, gg, width = 7.2, height = 3.6, dpi = 300)
  }
  gg
}


# ---------------------------------------------------------
# Panel (b) : Treatment Effect (full pre/post period)
# Plots treatment effects (Obs - CF) for both Proposed and SCM.
# Includes logic to skip a year and connect the gap.
# ---------------------------------------------------------
plot_panel_effect <- function(
  fit,
  data,
  treated_unit,
  treatment_dummy,
  y,
  unit_col,
  time_col,
  cred = 0.95,
  skip_year = NULL, # Year (e.g., 2011) to skip
  gap_linetype = "dotted", # Linetype for the skipped gap
  vline_at = NULL,
  annotations_df = NULL,
  x_title = "Year",
  y_title = "Treatment Effect",
  year_every = 2L
) {
  # 1. Calculate treatment effect for the Proposed method
  cf_prop <- scspill_counterfactual(fit, cred = cred, time_col = time_col)
  eff_prop <- data.frame(
    time = cf_prop$time,
    .idx = cf_prop$.idx,
    mean = cf_prop$y_obs - cf_prop$y_cf_mean,
    lo = cf_prop$y_obs - cf_prop$y_cf_hi,
    hi = cf_prop$y_obs - cf_prop$y_cf_lo,
    series = "Proposed"
  )

  # 2. Calculate treatment effect for SCM
  scm_cf_full <- scm_counterfactual_light(
    data = data,
    treated_unit = treated_unit,
    treatment_dummy = treatment_dummy,
    y = y,
    unit_col = unit_col,
    time_col = time_col
  )
  scm_merged <- merge(
    cf_prop[, c("time", ".idx", "y_obs")],
    scm_cf_full[, c("time", "y_cf")],
    by = "time",
    all.x = TRUE,
    sort = FALSE
  )
  eff_scm <- data.frame(
    time = scm_merged$time,
    .idx = scm_merged$.idx,
    mean = scm_merged$y_obs - scm_merged$y_cf,
    series = "SCM"
  )
  eff_scm <- eff_scm[order(eff_scm$.idx), ]

  # 3. Setup axis breaks and labels
  xax <- .make_year_breaks(eff_prop$time, every = year_every)
  idx_event <- if (is.null(vline_at)) {
    NA_integer_
  } else {
    m <- match(vline_at, eff_prop$time)
    if (is.na(m)) NA_integer_ else as.integer(eff_prop$.idx[m])
  }

  # 4. Handle 'skip_year': Create separate data.frames for gap
  eff_prop_gap <- data.frame()
  eff_scm_gap <- data.frame()

  if (!is.null(skip_year)) {
    time_skip <- skip_year
    time_bridge <- c(time_skip - 1, time_skip + 1)

    idx_skip <- eff_prop$.idx[which(eff_prop$time == time_skip)]
    idx_bridge <- eff_prop$.idx[which(eff_prop$time %in% time_bridge)]

    if (length(idx_skip) > 0 && length(idx_bridge) == 2) {
      # 4a. Create data.frames for the gap lines
      eff_prop_gap <- eff_prop[eff_prop$.idx %in% idx_bridge, ]
      eff_scm_gap <- eff_scm[eff_scm$.idx %in% idx_bridge, ]

      # 4b. Set the 'mean' (y-value) of the skipped year to NA
      eff_prop$mean[eff_prop$.idx == idx_skip] <- NA
      eff_prop$lo[eff_prop$.idx == idx_skip] <- NA
      eff_prop$hi[eff_prop$.idx == idx_skip] <- NA
      eff_scm$mean[eff_scm$.idx == idx_skip] <- NA
    }
  }

  # 5. Construct the ggplot object
  gg <- ggplot(eff_prop, aes(.idx, mean)) +
    # Main CI ribbon (will have a gap)
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = "95% Credible Interval"),
      alpha = 0.22
    ) +
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    # Main lines (will have a gap)
    geom_line(aes(color = "Proposed", linetype = "Proposed"), linewidth = 0.7) +
    geom_line(
      data = eff_scm,
      aes(.idx, mean, color = "SCM", linetype = "SCM"),
      linewidth = 0.6
    ) +

    # Intervention line
    {
      if (!is.na(idx_event)) geom_vline(xintercept = idx_event, linetype = 1)
    } +

    theme_scspill_paper() +
    labs(x = x_title, y = y_title, title = "(b) Treatment Effect Estimation") +
    scale_x_continuous(breaks = xax$breaks, labels = xax$labels) +
    scale_color_manual(
      name = NULL,
      values = c("Proposed" = paper_colors$proposed, "SCM" = paper_colors$scm),
      breaks = c("Proposed", "SCM")
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("Proposed" = "solid", "SCM" = "dashed"),
      breaks = c("Proposed", "SCM")
    ) +
    scale_fill_manual(
      name = NULL,
      values = c("95% Credible Interval" = paper_colors$ribbon)
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      linetype = guide_legend(order = 1),
      fill = guide_legend(order = 2, nrow = 1)
    )

  # 6. Add annotations (arrows and text)
  if (is.data.frame(annotations_df) && nrow(annotations_df) > 0) {
    ann <- annotations_df
    ann$.idx <- eff_prop$.idx[match(ann$time, eff_prop$time)]
    ann <- ann[!is.na(ann$.idx), ]
    if (nrow(ann) > 0) {
      gg <- gg +
        geom_segment(
          data = ann,
          aes(x = .idx - 0.8, xend = .idx - 0.1, y = y, yend = y),
          arrow = arrow(length = unit(0.08, "in")),
          inherit.aes = FALSE
        ) +
        geom_text(
          data = ann,
          aes(
            x = .idx - 1.0,
            y = y,
            label = label,
            hjust = 1,
            vjust = vjust %||% 0.5
          ),
          size = 3,
          inherit.aes = FALSE
        )
    }
  }
  interpolation_layers_eff <- list(
    if (nrow(eff_prop_gap) > 0) {
      geom_line(
        data = eff_prop_gap,
        aes(color = "Proposed"),
        linetype = gap_linetype,
        linewidth = 0.7
      )
    },

    if (nrow(eff_scm_gap) > 0) {
      geom_line(
        data = eff_scm_gap,
        aes(color = "SCM"),
        linetype = gap_linetype,
        linewidth = 0.6
      )
    }
  )

  gg <- gg + interpolation_layers_eff
  gg
}

.standardize_spill <- function(spill_raw_matrix) {
  M <- as.matrix(spill_raw_matrix)
  T_total <- nrow(M)
  N_units <- ncol(M)
  if (T_total == 0 || N_units == 0) {
    return(data.frame())
  }

  times <- rownames(M)
  if (is.null(times)) {
    times <- seq_len(T_total)
  } else {
    times <- as.numeric(times)
  }

  units <- colnames(M)
  if (is.null(units)) {
    units <- paste0("unit_", seq_len(N_units))
  }

  df <- data.frame(
    time = rep(times, times = N_units),
    unit = rep(units, each = T_total),
    mean = as.numeric(M)
  )
  df$unit <- as.character(df$unit)
  df
}

tidy_scspill <- function(fit, time_col = NULL) {
  spill_df <- data.frame()
  if (!is.null(fit$effects) && !is.null(fit$effects$spill)) {
    spill_df <- .standardize_spill(fit$effects$spill)
  }

  weights_df <- data.frame()
  if (!is.null(fit$alpha_draws)) {
    A <- as.matrix(fit$alpha_draws)
    if (ncol(A) > 0) {
      alpha_hat <- colMeans(A, na.rm = TRUE)
      units <- colnames(fit$alpha_draws) %||%
        colnames(fit$inputs$Yc_pre) %||%
        paste0("unit_", length(alpha_hat))
      weights_df <- data.frame(
        unit = as.character(units),
        alpha = as.numeric(alpha_hat),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    spill = spill_df,
    weights = weights_df
  )
}

# ---------------------------------------------------------
# Panel (c) : Spillover Effects (full pre/post period)
# Plots spillover effects for all control units.
# Includes logic to skip a year and connect the gap.
# ---------------------------------------------------------
plot_spillover_panel <- function(
  fit,
  units = NULL,
  annotations_df = NULL,
  skip_year = NULL, # Year (e.g., 2011) to skip
  gap_linetype = "dotted", # Linetype for the skipped gap
  time_col = NULL,
  main_title = "Estimates of the Spillover Effects",
  x_title = "Year",
  y_title = "Spillover Effect",
  vline_at = NULL
) {
  # 1. Get tidy spillover data
  td <- tidy_scspill(fit, time_col = time_col)
  df <- td$spill
  if (!nrow(df)) {
    stop(
      "effects$spill must be provided as (T0+T1) x N matrix with time rownames."
    )
  }

  # 2. Filter units if specified
  if (!is.null(units)) {
    df <- df[df$unit %in% units, , drop = FALSE]
  }

  # 3. Handle 'skip_year': Create separate data.frames for gap
  df_gap <- data.frame()

  if (!is.null(skip_year)) {
    time_skip <- skip_year
    time_bridge <- c(time_skip - 1, time_skip + 1)

    # 3a. Create data.frame for the gap lines
    df_gap <- df[df$time %in% time_bridge, ]

    # 3b. Set the 'mean' (y-value) of the skipped year to NA
    idx_skip_rows <- which(df$time == time_skip)
    if (length(idx_skip_rows) > 0) {
      df$mean[idx_skip_rows] <- NA
    }
  }

  # 4. Construct the ggplot object
  gg <- ggplot() +
    geom_line(
      data = df,
      aes(time, mean, group = unit),
      alpha = 0.6,
      color = "grey40",
      linetype = "solid"
    ) +

    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.3) +
    labs(x = x_title, y = y_title, title = main_title) +
    theme_scspill_paper() +
    theme(legend.position = "none")

  # Intervention line
  if (!is.null(vline_at)) {
    gg <- gg + geom_vline(xintercept = vline_at, linetype = 1)
  }

  # 5. Annotations
  if (is.data.frame(annotations_df) && nrow(annotations_df) > 0) {
    gg <- gg +
      geom_segment(
        data = annotations_df,
        aes(x = time - 0.8, xend = time - 0.1, y = y-0.2, yend = y),
        arrow = arrow(length = unit(0.08, "in")),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = annotations_df,
        aes(
          x = time - 1.0,
          y = y-0.2,
          label = label,
          hjust = 1,
          vjust = vjust %||% 0.5
        ),
        size = 3,
        inherit.aes = FALSE
      )
  }

  interpolation_layer_df <- list(
    if (nrow(df_gap) > 0) {
      geom_line(
        data = df_gap,
        aes(time, mean, group = unit),
        alpha = 0.6,
        color = "grey40",
        linetype = gap_linetype
      )
    }
  )

  gg <- gg + interpolation_layer_df
  gg
}


.get_times_pre <- function(fit, time_col = NULL) {
  if (
    !is.null(time_col) &&
      is.data.frame(fit$inputs$data_pre) &&
      time_col %in% names(fit$inputs$data_pre)
  ) {
    return(as.vector(fit$inputs$data_pre[[time_col]]))
  }
  fit$inputs$times_pre %||%
    rownames(fit$inputs$Yc_pre) %||%
    seq_len(nrow(fit$inputs$Yc_pre))
}
.get_times_post <- function(fit, time_col = NULL) {
  if (
    !is.null(time_col) &&
      is.data.frame(fit$inputs$data_post) &&
      time_col %in% names(fit$inputs$data_post)
  ) {
    return(as.vector(fit$inputs$data_post[[time_col]]))
  }
  fit$inputs$times_post %||%
    rownames(fit$inputs$Yc_post) %||%
    seq(
      from = length(.get_times_pre(fit)) + 1,
      length.out = nrow(fit$inputs$Yc_post)
    )
}

scspill_counterfactual <- function(fit, cred = 0.95, time_col = NULL) {
  stopifnot(inherits(fit, "scspill"))
  Y0_pre <- as.numeric(fit$inputs$Y0_pre)
  Y0_post <- as.numeric(fit$inputs$Y0_post)
  Yc_pre <- as.matrix(fit$inputs$Yc_pre)
  Yc_post <- as.matrix(fit$inputs$Yc_post)

  N <- ncol(Yc_pre)
  T0 <- nrow(Yc_pre)
  T1 <- nrow(Yc_post)
  w <- as.matrix(fit$inputs$w)
  W <- as.matrix(fit$inputs$W)

  alpha_draws <- as.matrix(fit$alpha_draws)
  rho_draws <- as.numeric(fit$rho_draws)
  M <- nrow(alpha_draws)
  IN <- diag(N)

  times_pre <- .get_times_pre(fit, time_col)
  times_post <- .get_times_post(fit, time_col)

  ycf_pre_draws <- matrix(NA_real_, T0, M)
  ycf_post_draws <- matrix(NA_real_, T1, M)

  for (m in seq_len(M)) {
    a <- alpha_draws[m, ]
    r <- rho_draws[m]
    Ainv <- tryCatch(
      solve(IN - r * (w %*% t(a) + W)),
      error = function(e) {
        warning(paste(
          "Solver failed for iter",
          m,
          "rho=",
          r,
          ". Returning NA."
        ))
        return(matrix(NA_real_, N, N))
      }
    )
    if (anyNA(Ainv)) {
      next
    }

    B <- (IN - r * W)
    for (t in seq_len(T0)) {
      tmp <- Ainv %*% (B %*% Yc_pre[t, ] - r * w * Y0_pre[t])
      ycf_pre_draws[t, m] <- as.numeric(crossprod(a, tmp))
    }
    for (t in seq_len(T1)) {
      tmp <- Ainv %*% (B %*% Yc_post[t, ] - r * w * Y0_post[t])
      ycf_post_draws[t, m] <- as.numeric(crossprod(a, tmp))
    }
  }

  lo <- (1 - cred) / 2
  hi <- 1 - lo
  df_pre <- data.frame(
    time = times_pre,
    t_idx = seq_len(T0),
    period = "pre",
    y_obs = Y0_pre,
    y_cf_mean = rowMeans(ycf_pre_draws, na.rm = TRUE),
    y_cf_lo = apply(
      ycf_pre_draws,
      1,
      stats::quantile,
      probs = lo,
      na.rm = TRUE
    ),
    y_cf_hi = apply(
      ycf_pre_draws,
      1,
      stats::quantile,
      probs = hi,
      na.rm = TRUE
    )
  )
  df_post <- data.frame(
    time = times_post,
    t_idx = T0 + seq_len(T1),
    period = "post",
    y_obs = Y0_post,
    y_cf_mean = rowMeans(ycf_post_draws, na.rm = TRUE),
    y_cf_lo = apply(
      ycf_post_draws,
      1,
      stats::quantile,
      probs = lo,
      na.rm = TRUE
    ),
    y_cf_hi = apply(
      ycf_post_draws,
      1,
      stats::quantile,
      probs = hi,
      na.rm = TRUE
    )
  )

  # 結合
  cf_full <- rbind(df_pre, df_post)

  cf_full$.idx <- seq_len(nrow(cf_full))

  return(cf_full)
}

scm_counterfactual_light <- function(
  data,
  treated_unit,
  treatment_dummy,
  y,
  unit_col,
  time_col
) {
  stopifnot(is.data.frame(data))
  df <- data

  ord <- order(df[[unit_col]], df[[time_col]])
  df <- df[ord, , drop = FALSE]

  is_treated <- df[[unit_col]] == treated_unit
  is_pre <- df[[treatment_dummy]] == 0
  is_post <- df[[treatment_dummy]] == 1

  years_pre <- sort(unique(df[[time_col]][is_treated & is_pre]))
  years_post <- sort(unique(df[[time_col]][is_treated & is_post]))
  T0 <- length(years_pre)
  T1 <- length(years_post)
  if (T0 == 0) {
    stop("No pre-treatment periods found for treated unit.")
  }

  y_tr_pre <- df[[y]][is_treated & is_pre][order(df[[time_col]][
    is_treated & is_pre
  ])]
  y_tr_post <- df[[y]][is_treated & is_post][order(df[[time_col]][
    is_treated & is_post
  ])]

  donors <- setdiff(unique(df[[unit_col]]), treated_unit)
  if (length(donors) == 0) {
    stop("No donor units found.")
  }

  .get_donor_series <- function(u, target_years) {
    idx_u <- (df[[unit_col]] == u)
    if (!any(idx_u)) {
      return(rep(NA_real_, length(target_years)))
    }
    vv_u <- df[[y]][idx_u]
    tt_u <- df[[time_col]][idx_u]
    ord_u <- order(tt_u)
    vv_u_sorted <- vv_u[ord_u]
    tt_u_sorted <- tt_u[ord_u]
    idx_match <- match(target_years, tt_u_sorted)
    vv_u_sorted[idx_match]
  }

  X_pre <- sapply(donors, .get_donor_series, target_years = years_pre)
  X_post <- sapply(donors, .get_donor_series, target_years = years_post)

  fill_cols <- function(X) {
    X <- as.matrix(X)
    if (nrow(X) == 0) {
      return(X)
    }
    for (j in seq_len(ncol(X))) {
      v <- X[, j]
      if (anyNA(v)) {
        ok <- which(!is.na(v))
        if (length(ok) >= 2) {
          v <- approx(x = ok, y = v[ok], xout = seq_along(v), rule = 2)$y
        } else if (length(ok) == 1) {
          v[is.na(v)] <- v[ok[1]]
        } else {
          v[is.na(v)] <- 0
        }
      }
      X[, j] <- v
    }
    X
  }
  X_pre <- fill_cols(X_pre)
  X_post <- fill_cols(X_post)

  if (nrow(X_pre) != length(y_tr_pre)) {
    stop("Dimension mismatch in pre-period (SCM).")
  }

  w_hat <- compute_scm_weights(y_tr_pre, X_pre)
  y_cf_pre <- y_tr_pre
  y_cf_post <- as.numeric(X_post %*% w_hat)

  data.frame(
    time = c(years_pre, years_post),
    y_cf = c(y_cf_pre, y_cf_post),
    period = c(rep("pre", T0), rep("post", T1))
  )
}


create_spillover_annotations <- function(
  fit,
  time_point,
  units_to_label,
  time_offset = -1.0,
  vjust_list = NULL
) {
  td <- tidy_scspill(fit, time_col = "year")
  df <- td$spill

  if (!nrow(df)) {
    warning("create_spillover_annotations: fit$effects$spill が空です。")
    return(data.frame())
  }

  df_at_time <- df[df$time == time_point, ]
  if (!nrow(df_at_time)) {
    warning(paste(
      "指定された年",
      time_point,
      "のデータが spill データに見つかりません。"
    ))
    return(data.frame())
  }

  ann_df <- data.frame(
    label = units_to_label,
    time = time_point,
    stringsAsFactors = FALSE
  )

  ann_df <- merge(
    ann_df,
    df_at_time[, c("unit", "mean")],
    by.x = "label",
    by.y = "unit",
    all.x = TRUE
  )

  names(ann_df)[names(ann_df) == "mean"] <- "y"

  ann_df$hjust <- 1
  ann_df$vjust <- if (is.null(vjust_list)) 0.5 else vjust_list

  return(ann_df)
}

compute_scm_weights <- function(Y0_pre, Yc_pre) {
   y  <- as.numeric(Y0_pre)
  X  <- as.matrix(Yc_pre)
  N  <- ncol(X)
  
  # 目的関数: min ||X w - y||^2
  E <- matrix(1, nrow = 1, ncol = N)    # E w = f  （和＝1）
  f <- 1
  G <- diag(N)                          # G w >= h （w >= 0）
  h <- rep(0, N)
  
  fit <- limSolve::lsei(A = X, B = y,
                        E = E, F = f,
                        G = G, H = h)
  as.numeric(fit$X)
}
