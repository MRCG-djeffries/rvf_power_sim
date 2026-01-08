plots_out1=function(){
  library(ggplot2)
  library(dplyr)
  
  # -----------------------------
  # Result directories (FIXED)
  # -----------------------------
  dirout = c("~/Documents/ed_YF3/case_split/shiny_graph/www/out9","~/Documents/ed_YF3/case_split/shiny_graph/www/out11")
  XA=load_results(dirout[1])
  XB=load_results(dirout[2])
  plot1(0.1,XA)
  #B01=gettab(XB,1.0)
  # create matirx with cols seropos and rows the VE
  #avg_line <- add_sero_average(XB)
  #XB2 <- dplyr::bind_rows(XB, avg_line)   # XB plus an “Average” pseudo-seropos
  
  plot_seroneg_vs_ve_by_seropos_black(XB,XA,0.9)
  plot_seroneg_vs_ve_by_seropos_all(XB,0.9)
}

plot_seroneg_vs_ve_by_seropos_black <- function(df, df_base,
                                                span = 0.8,
                                                sero_lines = c(10, 20, 30, 40, 50),
                                                ve_levels = seq(0.75, 0.95, by = 0.02),
                                                ar_levels = c(0.001, 0.01)) {
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  # ---- helper: interpolate y over sero grid to desired sero_lines ----
  interp_sero <- function(sero_grid, y_grid, xout) {
    ok <- is.finite(sero_grid) & is.finite(y_grid)
    sero_grid <- sero_grid[ok]; y_grid <- y_grid[ok]
    o <- order(sero_grid)
    sero_grid <- sero_grid[o]; y_grid <- y_grid[o]
    if (length(sero_grid) < 2) return(rep(NA_real_, length(xout)))
    approx(x = sero_grid, y = y_grid, xout = xout, rule = 2)$y
  }
  
  # ---- base grid ----
  d0 <- df %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  d0_base <- df_base %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  stopifnot(nrow(d0) > 0)
  
  # ---- interpolate to only the chosen sero lines at each (AR, VE) ----
  d_interp <- d0 %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_pct = sero_lines,
      y = interp_sero(sero_grid, y_grid, sero_lines),
      .groups = "drop"
    ) %>%
    mutate(sero_lab = paste0(sero_pct, "%")) %>%
    filter(is.finite(y))
  
  # ---- average line (mean over ALL available sero_grid in df_base) ----
  d_avg <- d0_base %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_lab = "Independent of risk",
      y = mean(y_grid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(y))
  
  # combine
  d_all <- bind_rows(
    d_interp %>% select(ar, AR_lab, VE_pct, ve_true, sero_lab, y),
    d_avg    %>% select(ar, AR_lab, VE_pct, ve_true, sero_lab, y)
  )
  
  # Put legend in the order you want
  d_all$sero_lab <- factor(
    d_all$sero_lab,
    levels = c(paste0(sero_lines, "%"), "Independent of risk")
  )
  
  # ---- manual colour map: sero_lines use default hue palette; Independent forced to black ----
  pal <- scales::hue_pal()(length(sero_lines))
  cols <- setNames(
    c(pal, "black"),
    c(paste0(sero_lines, "%"), "Independent of risk")
  )
  
  p=ggplot(d_all, aes(x = VE_pct, y = y, group = sero_lab)) +
    geom_smooth(
      aes(colour = sero_lab, linetype = sero_lab),
      method = "loess",
      span = span,
      se = FALSE,
      linewidth = 1.0
    ) +
    facet_wrap(~ AR_lab, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(75, 95, by = 2), limits = c(75, 95)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_colour_manual(values = cols, drop = FALSE) +
    scale_linetype_manual(
      values = setNames(rep("solid", length(levels(d_all$sero_lab))), levels(d_all$sero_lab)),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(title = "Baseline seropositivity",
                            override.aes = list(linewidth = 1.2)),
      linetype = "none"
    ) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Sero-negative person-years per arm",
      #title = "Required seronegative enrollment per arm vs VE"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  print(p)
  invisible(p)
}

plot_seroneg_vs_ve_by_seropos_all <- function(df,
                                              span = 0.8,
                                              sero_lines = c(10, 20, 30, 40, 50),
                                              ve_levels = seq(0.75, 0.95, by = 0.02),
                                              ar_levels = c(0.001, 0.01)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---- helper: interpolate y over sero grid to desired sero_lines ----
  interp_sero <- function(sero_grid, y_grid, xout) {
    ok <- is.finite(sero_grid) & is.finite(y_grid)
    sero_grid <- sero_grid[ok]
    y_grid    <- y_grid[ok]
    o <- order(sero_grid)
    sero_grid <- sero_grid[o]
    y_grid    <- y_grid[o]
    
    if (length(sero_grid) < 2) {
      return(rep(ifelse(length(y_grid) == 0, NA_real_, y_grid[1]), length(xout)))
    }
    approx(x = sero_grid, y = y_grid, xout = xout, rule = 2)$y
  }
  
  # ---- base grid (your simulated points) ----
  d0 <- df %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),     # 75..95
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_total_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  stopifnot(nrow(d0) > 0)
  
  # ---- interpolate to ONLY 10/20/30/40/50 at each (AR, VE) ----
  d_interp <- d0 %>%
    group_by(ar, AR_lab, VE_pct) %>%
    summarise(
      sero_pct = sero_lines,
      y = interp_sero(sero_grid, y_grid, sero_lines),
      .groups = "drop"
    ) %>%
    mutate(
      sero_lab = paste0(sero_pct, "%")
    ) %>%
    arrange(ar, sero_pct, VE_pct)
  
  # label positions at right edge
  lab <- d_interp %>%
    group_by(AR_lab, sero_lab) %>%
    filter(VE_pct == max(VE_pct, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  p=ggplot(d_interp, aes(x = VE_pct, y = y, group = sero_lab, colour = sero_lab)) +
    #geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_smooth(method = "loess", span = span, se = FALSE, linewidth = 0.9) +
    geom_text(
      data = lab,
      aes(label = sero_lab),
      hjust = -0.1,
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~ AR_lab, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(75, 95, by = 2), limits = c(75, 96)) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Person years per arm",
      colour = NULL,  # <-- CHANGE 1: remove legend title (optional but tidy)
      #title = "Required enrollment (for sero -ve and +ve) per arm vs VE"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",  # <-- CHANGE 2: drop legend (keep line labels)
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(5.5, 30, 5.5, 5.5)
    )
  print(p)
  invisible(p)
}
plot_seroneg_vs_ve_by_seropos <- function(df, df_base,
                                          span = 0.8,
                                          sero_lines = c(10, 20, 30, 40, 50),
                                          ve_levels = seq(0.75, 0.95, by = 0.02),
                                          ar_levels = c(0.001, 0.01)) {
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  # ---- helper: interpolate y over sero grid to desired sero_lines ----
  interp_sero <- function(sero_grid, y_grid, xout) {
    ok <- is.finite(sero_grid) & is.finite(y_grid)
    sero_grid <- sero_grid[ok]; y_grid <- y_grid[ok]
    o <- order(sero_grid)
    sero_grid <- sero_grid[o]; y_grid <- y_grid[o]
    if (length(sero_grid) < 2) return(rep(NA_real_, length(xout)))
    approx(x = sero_grid, y = y_grid, xout = xout, rule = 2)$y
  }
  
  # ---- base grid ----
  d0 <- df %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  d0_base <- df_base %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  stopifnot(nrow(d0) > 0)
  
  # ---- interpolate to only the chosen sero lines at each (AR, VE) ----
  d_interp <- d0 %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_pct = sero_lines,
      y = interp_sero(sero_grid, y_grid, sero_lines),
      .groups = "drop"
    ) %>%
    mutate(sero_lab = paste0(sero_pct, "%")) %>%
    filter(is.finite(y))
  
  # ---- average line (mean over ALL available sero_grid in df_base) ----
  d_avg <- d0_base %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_lab = "Independent of risk",
      y = mean(y_grid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(y))
  
  # combine (keep Average as its own group)
  d_all <- bind_rows(
    d_interp %>% select(ar, AR_lab, VE_pct, ve_true, sero_lab, y),
    d_avg    %>% select(ar, AR_lab, VE_pct, ve_true, sero_lab, y)
  )
  
  # Put legend in the order you want
  d_all$sero_lab <- factor(
    d_all$sero_lab,
    levels = c(paste0(sero_lines, "%"), "Independent of risk")
  )
  
  ggplot(d_all, aes(x = VE_pct, y = y, group = sero_lab)) +
    # coloured smoothed lines for seropos curves + average included
    geom_smooth(
      aes(colour = sero_lab, linetype = sero_lab),
      method = "loess",
      span = span,
      se = FALSE,
      linewidth = 1.0
    ) +
    facet_wrap(~ AR_lab, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(75, 95, by = 2), limits = c(75, 95)) +
    scale_y_continuous(labels = scales::label_comma()) +
    # Make Average black, others default palette; also make Average solid
    scale_colour_discrete(drop = FALSE) +
    scale_linetype_manual(
      values = c(rep("solid", length(sero_lines)), "solid"),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(title = "Baseline seropositivity", override.aes = list(linewidth = 1.2)),
      linetype = "none"
    ) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Sero-negative person-years per arm",
      title = "Required seronegative enrollment per arm vs VE"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}


plot_seroneg_vs_ve_by_seropos_v0 <- function(df,df_base,
                                          span = 0.8,
                                          sero_lines = c(10, 20, 30, 40, 50),
                                          ve_levels = seq(0.75, 0.95, by = 0.02),
                                          ar_levels = c(0.001, 0.01)) {
  library(dplyr)
  library(ggplot2)
  
  # ---- helper: interpolate y over sero grid to desired sero_lines ----
  interp_sero <- function(sero_grid, y_grid, xout) {
    ok <- is.finite(sero_grid) & is.finite(y_grid)
    sero_grid <- sero_grid[ok]; y_grid <- y_grid[ok]
    o <- order(sero_grid)
    sero_grid <- sero_grid[o]; y_grid <- y_grid[o]
    if (length(sero_grid) < 2) return(rep(NA_real_, length(xout)))
    approx(x = sero_grid, y = y_grid, xout = xout, rule = 2)$y
  }
  
  # ---- base grid ----
  d0 <- df %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  d0_base <- df_base %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  stopifnot(nrow(d0) > 0)
  
  # ---- interpolate to only the chosen sero lines at each (AR, VE) ----
  d_interp <- d0 %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_pct = sero_lines,
      y = interp_sero(sero_grid, y_grid, sero_lines),
      .groups = "drop"
    ) %>%
    mutate(sero_lab = paste0(sero_pct, "%")) %>%
    filter(is.finite(y))
  
  # ---- average line (mean over ALL available sero_grid in the data) ----
  d_avg <- d0_base %>%
    group_by(ar, AR_lab, VE_pct, ve_true) %>%
    summarise(
      sero_pct = -999,
      sero_lab = "Average",
      y = mean(y_grid, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(y))
  d_all <- bind_rows(d_interp, d_avg)
  
  # ---- label positions at right edge (use last VE point per curve) ----
  lab <- d_all %>%
    group_by(AR_lab, sero_lab) %>%
    filter(VE_pct == max(VE_pct, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  # ---- build plot: ONLY smoothed lines ----
  ggplot(d_all, aes(x = VE_pct, y = y, group = sero_lab)) +
    # smoothed coloured lines for seropos curves
    geom_smooth(
      data = d_all %>% filter(sero_lab != "Average"),
      aes(colour = sero_lab),
      method = "loess", span = span, se = FALSE, linewidth = 0.9
    ) +
    # smoothed black average line
    geom_smooth(
      data = d_all %>% filter(sero_lab == "Average"),
      method = "loess", span = span, se = FALSE, linewidth = 1.1,
      colour = "black"
    ) +
    # line labels at right edge (no legend)
    geom_text(
      data = lab,
      aes(label = sero_lab, colour = ifelse(sero_lab == "Average", NA, sero_lab)),
      hjust = -0.1,
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~ AR_lab, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(75, 95, by = 2), limits = c(75, 96)) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Sero -ve person years per arm",
      title = "Required seronegative enrollment per arm vs VE"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(5.5, 30, 5.5, 5.5)
    )
}



plot_seroneg_vs_ve_by_seropos_old <- function(df,
                                          span = 0.8,
                                          sero_lines = c(10, 20, 30, 40, 50),
                                          ve_levels = seq(0.75, 0.95, by = 0.02),
                                          ar_levels = c(0.001, 0.01)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---- helper: interpolate y over sero grid to desired sero_lines ----
  interp_sero <- function(sero_grid, y_grid, xout) {
    ok <- is.finite(sero_grid) & is.finite(y_grid)
    sero_grid <- sero_grid[ok]
    y_grid    <- y_grid[ok]
    o <- order(sero_grid)
    sero_grid <- sero_grid[o]
    y_grid    <- y_grid[o]
    
    if (length(sero_grid) < 2) {
      return(rep(ifelse(length(y_grid) == 0, NA_real_, y_grid[1]), length(xout)))
    }
    approx(x = sero_grid, y = y_grid, xout = xout, rule = 2)$y
  }
  
  # ---- base grid (your simulated points) ----
  d0 <- df %>%
    mutate(
      AR_pct    = 100 * ar,
      VE_pct    = round(100 * ve_true),     # 75..95
      sero_grid = round(sero_pct),
      y_grid    = as.numeric(N_seroneg_per_arm)
    ) %>%
    filter(
      ar %in% ar_levels,
      round(ve_true, 2) %in% ve_levels,
      is.finite(y_grid)
    ) %>%
    mutate(
      AR_lab = factor(
        AR_pct,
        levels = sort(unique(AR_pct)),
        labels = paste0("AR = ", sort(unique(AR_pct)), "%")
      )
    ) %>%
    arrange(ar, VE_pct, sero_grid)
  
  stopifnot(nrow(d0) > 0)
  
  # ---- interpolate to ONLY 10/20/30/40/50 at each (AR, VE) ----
  d_interp <- d0 %>%
    group_by(ar, AR_lab, VE_pct) %>%
    summarise(
      sero_pct = sero_lines,
      y = interp_sero(sero_grid, y_grid, sero_lines),
      .groups = "drop"
    ) %>%
    mutate(
      sero_lab = paste0(sero_pct, "%")
    ) %>%
    arrange(ar, sero_pct, VE_pct)
  
  # label positions at right edge
  lab <- d_interp %>%
    group_by(AR_lab, sero_lab) %>%
    filter(VE_pct == max(VE_pct, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  ggplot(d_interp, aes(x = VE_pct, y = y, group = sero_lab, colour = sero_lab)) +
    geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_smooth(method = "loess", span = span, se = FALSE, linewidth = 0.9) +
    geom_text(
      data = lab,
      aes(label = sero_lab),
      hjust = -0.1,
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~ AR_lab, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq(75, 95, by = 2), limits = c(75, 96)) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Number of seronegative subjects per arm",
      colour = NULL,  # <-- CHANGE 1: remove legend title (optional but tidy)
      title = "Required seronegative enrollment per arm vs VE"
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",  # <-- CHANGE 2: drop legend (keep line labels)
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(5.5, 30, 5.5, 5.5)
    )
}



# Example usage:
# p <- plot_seroneg_vs_ve_by_seropos(df, span = 0.6)
# print(p)

gettab = function(X,arval){
  d <-  X[X$ar_pct==arval,]
  req(nrow(d) > 0)
  
  d$VE_pct   <- round(100 * d$ve_true)
  d$sero_pct <- round(d$sero_pct)
  
  d$value <- d$N_seroneg_per_arm
  
  tab <- xtabs(value ~ sero_pct + VE_pct, data = d)
  out <- as.data.frame.matrix(tab)
  out <- cbind(sero_pos_pct = as.integer(rownames(out)), out)
  rownames(out) <- NULL
  out
}
format_metric_vec <- function(metric, x) {
  if (metric %in% c("power","p_stop_1","p_stop_2","p_stop_any","casesT")) return(round(x, 3))
  if (is_integer_metric(metric)) return(as.integer(ceiling(x)))  # ALWAYS round UP for integer metrics
  round(x, 3)
}
plot1=function(arval,X){
  #arval is 0.1 or 1
  X=X[X$ar_pct==arval,]
  x=X[,"ve_true"]
  y1=X[,"nevents"]
  y2=X[,"av_events"]
  #plot_events_vs_ve_safe_smooth_grid(x,y1,y2,0.4)
  p=plot_events_vs_ve_safe_smooth_grid_gg(x,y1,y2)
  print(p)
  invisible(p)
  # thi splots events for VE of 75% to 95% againts total number of eevents
  # for AMx events and expected early stopping
}

plot_events_vs_ve_safe_smooth_grid <- function(x, y1, y2,
                                               spar_grid = seq(0.3, 1.0, by = 0.05)) {
  
  stopifnot(length(x) == length(y1), length(x) == length(y2))
  
  ok <- is.finite(x) & is.finite(y1) & is.finite(y2)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]
  
  o <- order(x)
  x <- x[o]; y1 <- y1[o]; y2 <- y2[o]
  
  # Right-continuous step truth
  f1 <- stepfun(x[-length(x)], y1)
  f2 <- stepfun(x[-length(x)], y2)
  
  # ---- choose optimal spar ----
  choose_spar <- function(y, label = "") {
    errs <- sapply(spar_grid, function(sp) {
      s  <- smooth.spline(x, y, spar = sp)
      yh <- predict(s, x)$y
      yh <- pmax(yh, y)      # non-anti-conservative
      mean((yh - y)^2)
    })
    
    cat("\nSpline error grid", if (nzchar(label)) paste0(" (", label, ")"), ":\n")
    print(data.frame(
      spar = spar_grid,
      mse  = round(errs, 6)
    ))
    
    best <- spar_grid[which.min(errs)]
    cat("→ Selected spar =", best, "\n")
    
    best
  }
  
  spar1 <- choose_spar(y1)
  spar2 <- choose_spar(y2)
  
  # Dense grid
  xg <- seq(min(x), max(x), length.out = 600)
  
  # Final splines
  s1 <- smooth.spline(x, y1, spar = spar1)
  s2 <- smooth.spline(x, y2, spar = spar2)
  
  y1g <- pmax(predict(s1, xg)$y, f1(xg))
  y2g <- pmax(predict(s2, xg)$y, f2(xg))
  
  # Enforce monotone decrease
  y1g <- ceiling(cummin(y1g))
  y2g <- ceiling(cummin(y2g))
  
  # ---- plot ----
  plot(xg, y1g, type = "l", lwd = 2, col = "black",
       xlab = "VE (%)",
       ylab = "Total events",
       ylim = range(c(y1g, y2g), finite = TRUE))
  
  lines(xg, y2g, lwd = 2, col = "red")
  grid()
  
  legend("topright",
         legend = c("Max events", "Expected events"),
         col    = c("black", "red"),
         lwd    = 2,
         bty    = "n")
  
  invisible(list(spar_max = spar1, spar_exp = spar2))
}

plot_events_vs_ve_safe_smooth_grid_gg <- function(x, y1, y2,
                                                  spar_grid = seq(0.3, 1.0, by = 0.05)) {
  stopifnot(length(x) == length(y1), length(x) == length(y2))
  library(ggplot2)
  
  ok <- is.finite(x) & is.finite(y1) & is.finite(y2)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]
  
  o <- order(x)
  x <- x[o]; y1 <- y1[o]; y2 <- y2[o]
  
  # Right-continuous step truth
  f1 <- stepfun(x[-length(x)], y1)
  f2 <- stepfun(x[-length(x)], y2)
  
  # ---- choose optimal spar ----
  choose_spar <- function(y, label = "") {
    errs <- sapply(spar_grid, function(sp) {
      s  <- smooth.spline(x, y, spar = sp)
      yh <- predict(s, x)$y
      yh <- pmax(yh, y)  # never anti-conservative at the observed grid
      mean((yh - y)^2)
    })
    
    cat("\nSpline error grid", if (nzchar(label)) paste0(" (", label, ")"), ":\n")
    print(data.frame(spar = spar_grid, mse = round(errs, 6)))
    best <- spar_grid[which.min(errs)]
    cat("→ Selected spar =", best, "\n")
    best
  }
  
  spar1 <- choose_spar(y1, "Max events")
  spar2 <- choose_spar(y2, "Expected events")
  
  # Dense grid
  xg <- seq(min(x), max(x), length.out = 600)
  
  # Final splines
  s1 <- smooth.spline(x, y1, spar = spar1)
  s2 <- smooth.spline(x, y2, spar = spar2)
  
  y1g <- pmax(predict(s1, xg)$y, f1(xg))
  y2g <- pmax(predict(s2, xg)$y, f2(xg))
  
  # Enforce monotone decrease & integer counts
  y1g <- ceiling(cummin(y1g))
  y2g <- ceiling(cummin(y2g))
  
  # Build tidy df for ggplot
  dfp <- rbind(
    data.frame(VE_pct = 100 * xg, events = y1g, series = "Max events"),
    data.frame(VE_pct = 100 * xg, events = y2g, series = "Expected events")
  )
  
  ggplot(dfp, aes(x = VE_pct, y = events, colour = series)) +
    geom_line(linewidth = 1.1) +
    scale_x_continuous(
      breaks = seq(75, 95, by = 2),
      limits = c(75, 95)
    ) +
    scale_colour_manual(
      values = c("Max events" = "black", "Expected events" = "red")
    ) +
    labs(
      x = "Vaccine efficacy (VE, %)",
      y = "Total events",
      colour = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}



plot_events_vs_ve_safe_smooth <- function(x, y1, y2, spar = 0.75) {
  stopifnot(length(x) == length(y1), length(x) == length(y2))
  
  ok <- is.finite(x) & is.finite(y1) & is.finite(y2)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]
  
  o <- order(x)
  x <- x[o]; y1 <- y1[o]; y2 <- y2[o]
  
  # Dense grid
  xg <- seq(min(x), max(x), length.out = 600)
  
  # Step truth (right-continuous)
  f1 <- stepfun(x[-length(x)], y1)
  f2 <- stepfun(x[-length(x)], y2)
  y1_step <- f1(xg)
  y2_step <- f2(xg)
  
  # Smooth proposal
  s1 <- smooth.spline(x, y1, spar = spar)
  s2 <- smooth.spline(x, y2, spar = spar)
  y1_s <- predict(s1, xg)$y
  y2_s <- predict(s2, xg)$y
  
  # Clamp to never be anti-conservative
  y1_c <- pmax(y1_s, y1_step)
  y2_c <- pmax(y2_s, y2_step)
  
  # Enforce monotone non-increasing (VE increases -> required events should not increase)
  # Do this by taking cumulative minima from left->right
  y1_m <- cummin(y1_c)
  y2_m <- cummin(y2_c)
  
  # Round UP (counts)
  y1_m <- ceiling(y1_m)
  y2_m <- ceiling(y2_m)
  
  # Plot ONLY the final safe-smoothed lines
  plot(xg, y1_m, type = "l", lwd = 2,
       xlab = "VE (%)",
       ylab = "Total events",
       ylim = range(c(y1_m, y2_m), finite = TRUE))
  
  lines(xg, y2_m, lwd = 2, lty = 1,col="red")
  grid()
  legend("topright",
         legend = c("Max events", "Expected events"),
         col    = c("black", "red"),
         lwd    = 2,
         lty    = c(1, 1),   # or c(1,2) if you want dashed
         bty    = "n")
  
}




plot_events_vs_ve <- function(x, y1, y2, spar = 0.8) {
  stopifnot(length(x) == length(y1), length(x) == length(y2))
  
  # keep finite rows only
  ok <- is.finite(x) & is.finite(y1) & is.finite(y2)
  x <- x[ok]; y1 <- y1[ok]; y2 <- y2[ok]
  
  # sort by x
  o <- order(x)
  x <- x[o]; y1 <- y1[o]; y2 <- y2[o]
  
  # if duplicate x values exist, average them (spline wants unique-ish x)
  if (any(duplicated(x))) {
    y1 <- tapply(y1, x, mean)
    y2 <- tapply(y2, x, mean)
    x  <- as.numeric(names(y1))
    y1 <- as.numeric(y1)
    y2 <- as.numeric(y2)
    o <- order(x); x <- x[o]; y1 <- y1[o]; y2 <- y2[o]
  }
  
  # spline smoothing (ONLY smoothed lines will be plotted)
  s1 <- smooth.spline(x = x, y = y1, spar = spar)
  s2 <- smooth.spline(x = x, y = y2, spar = spar)
  
  xg <- seq(min(x), max(x), length.out = 300)
  y1g <- predict(s1, xg)$y
  y2g <- predict(s2, xg)$y
  
  
  # plot
  plot(xg, y1g, type = "l", lwd = 2,
       xlab = "VE (%)",
       ylab = "Total events",
       ylim = range(c(y1g, y2g), finite = TRUE))
  
  lines(xg, y2g, lwd = 2, lty = 2)
  
  grid()
  legend("topright",
         legend = c("Max events", "Expected events"),
         lwd = 2, lty = c(1, 2), bty = "n")
}

parse_bres_filename <- function(fname) {
  base <- basename(fname)
  base <- sub("\\.rds$", "", base)
  base <- sub("^Bres_", "", base)
  
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  stopifnot(length(parts) == 6)
  
  ar   <- as.numeric(paste0(parts[1], ".", parts[2]))
  hr   <- as.numeric(paste0(parts[3], ".", parts[4]))
  sero <- as.numeric(paste0(parts[5], ".", parts[6]))
  
  list(
    ar_pct   = ar,
    ar       = ar / 100,
    hr_true  = hr,
    ve_true  = 1 - hr,
    sero_pct = sero,
    seroprev = sero / 100
  )
}

load_results <- function(dir) {
  if (!dir.exists(dir)) stop("Directory not found: ", dir)
  
  files <- list.files(dir, pattern = "^Bres_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No Bres_*.rds files found in: ", dir)
  
  rows <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    meta <- parse_bres_filename(files[i])
    x <- readRDS(files[i])
    
    N_seroneg_per_arm <- x$Nvacc_mean
    N_total_per_arm   <- x$Nvacc_mean / (1 - meta$seroprev)
    
    p_stop_1   <- if (!is.null(x$p_stop_1)) x$p_stop_1 else NA_real_
    p_stop_2   <- if (!is.null(x$p_stop_2)) x$p_stop_2 else NA_real_
    p_stop_any <- if (is.finite(p_stop_1) && is.finite(p_stop_2)) (p_stop_1 + p_stop_2) else NA_real_
    
    rows[[i]] <- data.frame(
      file      = basename(files[i]),
      ar        = meta$ar,
      ar_pct    = meta$ar_pct,
      hr_true   = meta$hr_true,
      ve_true   = meta$ve_true,
      seroprev  = meta$seroprev,
      sero_pct  = meta$sero_pct,
      
      found     = isTRUE(x$found),
      nevents   = x$nevents,
      av_events = if (!is.null(x$av_events)) x$av_events else NA_real_,
      power     = x$power,
      looks     = I(list(x$looks)),
      
      N_seroneg_per_arm = N_seroneg_per_arm,
      N_total_per_arm   = N_total_per_arm,
      
      p_stop_1   = p_stop_1,
      p_stop_2   = p_stop_2,
      p_stop_any = p_stop_any,
      
      casesT     = if (!is.null(x$casesT)) x$casesT else NA_real_,
      
      stringsAsFactors = FALSE
    )
  }
  return(  do.call(rbind, rows))
}

add_sero_average <- function(df,
                             sero_col = "sero_pct",
                             y_col    = "N_seroneg_per_arm",
                             group_cols = c("ar","ar_pct","hr_true","ve_true","VE_pct")) {
  
  stopifnot(all(c(sero_col, y_col) %in% names(df)))
  
  d0 <- df
  if (!("VE_pct" %in% names(d0))) d0$VE_pct <- round(100 * d0$ve_true)
  
  # mean across seropos at each (AR, VE, HR, ...)
  davg <- d0 |>
    dplyr::mutate(y_grid = as.numeric(.data[[y_col]])) |>
    dplyr::filter(is.finite(y_grid)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols[group_cols %in% names(d0)]))) |>
    dplyr::summarise(
      sero_pct = -999,                  # sentinel
      sero_lab = "Average",
      y = mean(y_grid, na.rm = TRUE),
      .groups = "drop"
    )
  
  davg
}

