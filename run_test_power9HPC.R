#!/usr/bin/env Rscript
run_test_power9HPC = function(){
  # added increased seroprev in high risk clusters
  
  # ------------------------------------------------------------
  # NOTE ON SAWTOOTH POWER (exact discrete group-sequential tests)
  #
  # For the v4/v5 exact BIN-SPLIT design, the stopping boundaries
  # are recalibrated for EACH candidate total event target (nevents).
  # Because the test statistic is discrete, the achievable rejection
  # cutoffs jump in integer steps as nevents changes, which can make
  # estimated power vs nevents non-monotone ("sawtooth").
  #
  # Therefore: DO NOT use bisection/uniroot to search nevents for
  # the exact BIN-SPLIT design. Instead, use a small cached SCAN over
  # integer nevents values (here, a window around a Poisson starting
  # value read from the corresponding Poisson result file).
  # ------------------------------------------------------------
  # 19 Dec, this just uses total so should be quicker
  # array_id sets the row of parameters
  # B is the 
  # Now outputs the number of cases where the venets in vacc group were zero
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    stop("Usage: run_test_power4HPC.R <array_id 1..242> <B>")
  }
  array_id <- as.integer(args[[1]]) # 1 to 242 of below
  B        <- as.integer(args[[2]]) # number of sim reps per power calc
  f_dir_out <- args[[3]] # output directory
  testof    <- args[[4]] # "binsplit" or "poisson"
  f_dir_pois <- args[[5]]  # directory with Poisson results (for binsplit start)
  
  
  
  window = 5L; ne_floor = 3L
  
  ## -------- define parameter grid (5×5×5 = 125) --------
  wave_pct <- c(0.1,1)           # %
  hr_true  <- rev(seq(0.05,0.25,length.out=11))#c(0.25, 0.20, 0.15, 0.10, 0.05)
  sero_pct <- seq(10,50,length.out=11)#c(10, 20, 30, 40, 50)          # %
  
  grid <- expand.grid(
    wave = wave_pct,
    hr   = hr_true,
    sero = sero_pct,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (array_id < 1L || array_id > nrow(grid)) {
    stop("array_id out of range 1..", nrow(grid))
  }
  
  row <- grid[array_id, ]
  wave_attack <- row$wave / 100
  hr_val      <- row$hr
  seroprev    <- row$sero / 100
  
  cat(sprintf("\n=== Task %d of %d: wave_attack=%.2f%%, hr_true=%.2f, seroprev=%.2f%% ===\n",
              array_id, nrow(grid), 100*wave_attack, hr_val, 100*seroprev))
  OBF=getOBF_cuts()
  ## -------- file keys (used for output AND to read Poisson starting value) --------
  hr_str   <- gsub("\\.", "_", sprintf("%.2f", hr_val))  # hazard ratio
  ar_str   <- gsub("\\.", "_", sprintf("%.2f", row$wave)) # attack rate (%)
  sp_str   <- gsub("\\.", "_", sprintf("%.2f", row$sero)) # seropos (%)
  
  ## -------- run the search for this combo --------
  if (identical(testof, "binsplit")) {
    if (is.null(f_dir_pois)) stop("For binsplit v5 you must provide f_dir_pois (directory with Poisson Zres_*.rds files).")
    f_pois <- sprintf(paste0(f_dir_pois, "/Zres_%s_%s_%s.rds"), ar_str, hr_str, sp_str)
    if (!file.exists(f_pois)) stop("Poisson start file not found: ", f_pois)
    pois <- readRDS(f_pois)
    NpoisE <- as.integer(pois$nevents)  # Poisson-derived final (100%) look event target
    cat(sprintf("Using Poisson start NpoisE=%d from %s\n", NpoisE, f_pois))
    
    # Exact BIN-SPLIT has sawtooth power; use a cached window scan around NpoisE.
    res <- find_nevents_for_power_window("binsplit", OBF, Ncenter = NpoisE,
                                         window = window, ne_floor = ne_floor,
                                         target_power = 0.80, B = B, seed = 123L,
                                         HR0 = 0.80,
                                         waves = 12,
                                         wave_attack = wave_attack,
                                         hr_true = hr_val,
                                         seroprev = seroprev,
                                         frailty_var = 0.5,
                                         rho_target = 0.1,
                                         max_enrolled = 1e6)
  } else {
    # Poisson / score-test branch (approximately monotone in nevents): bisection is OK.
    res <- find_nevents_for_power("poisson", OBF,
                                  target_power = 0.80,
                                  B = B,
                                  ne_min = 10L,
                                  ne_max = 40L,
                                  seed = 123L,
                                  HR0 = 0.80,
                                  waves = 12,
                                  wave_attack = wave_attack,
                                  hr_true = hr_val,
                                  seroprev = seroprev,
                                  frailty_var = 0.5,
                                  rho_target = 0.1,
                                  max_enrolled = 1e6)
  }
  ## -------- save: file_wave_attack_hr_true_seroprev (hr uses '_' for '.') --------
  fname    <- sprintf(paste0(f_dir_out,"/",toupper(substr(testof,1,1)),"res_%s_%s_%s.rds"), ar_str, hr_str, sp_str)
  saveRDS(res, file = fname)
  cat("Saved:", fname, "\n")
  
  
}

## -------- search: find minimal nevents giving target power --------
# Find minimal nevents for a (possibly non-monotone) power curve by scanning
# a small integer window around a center value.
#
# For v5 BIN-SPLIT exact: Ncenter should come from the Poisson design output
# (final/100% look nevents), and this scan respects early stopping because
# power_at_nevents("binsplit", ...) runs the sequential monitoring.
find_nevents_for_power_window <- function(
    testof,
    cuts,
    target_power = 0.80,
    B = 3000,
    seed = 123L,
    Ncenter,
    window = 5L,
    ne_floor = 3L,
    expand_if_fail = TRUE,
    ne_cap = 40L,
    expand_step = 1L,
    keep_grid = TRUE,     # if FALSE, only minimal outputs (no stored grid)
    ...,
    verbose = TRUE
) {
  stopifnot(is.numeric(target_power), target_power > 0, target_power < 1)
  stopifnot(B >= 1, Ncenter >= 3)
  stopifnot(window >= 0, ne_floor >= 3, ne_cap >= ne_floor, expand_step >= 1)
  
  testof <- as.character(testof)
  Ncenter <- as.integer(Ncenter)
  window  <- as.integer(window)
  ne_floor <- as.integer(ne_floor)
  ne_cap <- as.integer(ne_cap)
  expand_step <- as.integer(expand_step)
  
  lo <- max(ne_floor, Ncenter - window)
  hi <- min(max(lo, Ncenter + window), ne_cap)
  
  if (verbose) {
    cat(sprintf("\n=== Window search (%s) target power=%.3f ===\n", testof, target_power))
    cat(sprintf("Scan order: nevents=%d..%d (center=%d, window=%d), B=%d\n",
                lo, hi, Ncenter, window, B))
    flush.console()
  }
  
  # In-process cache (safe on HPC: per-task memory only)
  cache <- new.env(parent = emptyenv())
  power_cached <- function(ne) {
    key <- as.character(ne)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    dum <- power_at_nevents(testof, cuts, as.integer(ne), B = B, seed = seed, ...)
    assign(key, dum, envir = cache)
    dum
  }
  
  # Optional diagnostics storage
  grid <- integer(0)
  pow  <- numeric(0)
  
  # ---- 1) scan upward within the initial window; BREAK at first success ----
  found_ne <- NA_integer_
  found_dum <- NULL
  
  for (ne in lo:hi) {
    dum <- power_cached(ne)
    p_ne <- dum$pow
    
    if (keep_grid) {
      grid <- c(grid, ne)
      pow  <- c(pow,  p_ne)
    }
    
    if (verbose) {
      cat(sprintf("nevents=%3d => power=%.3f\n", ne, p_ne)); flush.console()
    }
    
    if (is.finite(p_ne) && p_ne >= target_power) {
      found_ne <- ne
      found_dum <- dum
      break
    }
  }
  
  if (!is.na(found_ne)) {
    return(list(
      found = TRUE,
      nevents = found_ne,
      av_events = found_dum$av_events,
      power = found_dum$pow,
      looks = looks_from_nevents(found_ne),
      Nvacc_mean = ceiling(mean(found_dum$vaccN,  na.rm = TRUE)),
      sero_pos_vaccN_mean = ceiling(mean(found_dum$sero_pos_vaccN,  na.rm = TRUE)),
      Nvacc_median = ceiling(median(found_dum$vaccN,  na.rm = TRUE)),
      sero_pos_vaccN_median = ceiling(median(found_dum$sero_pos_vaccN,  na.rm = TRUE)),      
      p_stop_1 = mean(found_dum$stop1, na.rm = TRUE),
      p_stop_2 = mean(found_dum$stop2, na.rm = TRUE),
      casesT = mean(found_dum$casesT, na.rm = TRUE),
      grid = if (keep_grid) grid else NULL,
      power_grid = if (keep_grid) pow else NULL
    ))
  }
  
  # ---- 2) not found in window ----
  if (!isTRUE(expand_if_fail)) {
    if (!keep_grid) {
      return(list(found = FALSE, message = "Target power not reached in window (no expansion)."))
    }
    best_i <- which.max(pow)
    ne_best <- grid[best_i]
    dum_best <- power_cached(ne_best)
    return(list(
      found = FALSE,
      nevents = ne_best,
      av_events = dum_best$av_events,      
      power = dum_best$pow,
      looks = looks_from_nevents(ne_best),
      Nvacc_mean = ceiling(mean(dum_best$vaccN,  na.rm = TRUE)),
      sero_pos_vaccN_mean = ceiling(mean(dum_best$sero_pos_vaccN,  na.rm = TRUE)),
      Nvacc_median = ceiling(median(dum_best$vaccN,  na.rm = TRUE)),
      sero_pos_vaccN_median = ceiling(median(dum_best$sero_pos_vaccN,  na.rm = TRUE)),      
      p_stop_1 = mean(dum_best$stop1, na.rm = TRUE),
      p_stop_2 = mean(dum_best$stop2, na.rm = TRUE),
      casesT = mean(dum_best$casesT, na.rm = TRUE),
      grid = grid,
      power_grid = pow,
      message = "Target power not reached inside initial window (no expansion requested)."
    ))
  }
  
  if (verbose) {
    cat("No value in initial window reached target; expanding upward...\n")
    flush.console()
  }
  
  # ---- 3) expand upward; BREAK at first success ----
  ne <- hi + expand_step
  while (ne <= ne_cap) {
    dum <- power_cached(ne)
    p_ne <- dum$pow
    
    if (keep_grid) {
      grid <- c(grid, ne)
      pow  <- c(pow,  p_ne)
    }
    
    if (verbose) {
      cat(sprintf("expand nevents=%3d => power=%.3f\n", ne, p_ne)); flush.console()
    }
    
    if (is.finite(p_ne) && p_ne >= target_power) {
      return(list(
        found = TRUE,
        nevents = ne,
        av_events = dum$av_events,	
        power = p_ne,
        looks = looks_from_nevents(ne),
        Nvacc_mean = ceiling(mean(dum$vaccN,  na.rm = TRUE)),
        sero_pos_vaccN_mean = ceiling(mean(dum$sero_pos_vaccN,  na.rm = TRUE)),
        Nvacc_median = ceiling(median(dum$vaccN,  na.rm = TRUE)),
        sero_pos_vaccN_median = ceiling(median(dum$sero_pos_vaccN,  na.rm = TRUE)),        
        p_stop_1 = mean(dum$stop1, na.rm = TRUE),
        p_stop_2 = mean(dum$stop2, na.rm = TRUE),
        casesT = mean(dum$casesT, na.rm = TRUE),
        grid = if (keep_grid) grid else NULL,
        power_grid = if (keep_grid) pow else NULL
      ))
    }
    
    ne <- ne + expand_step
  }
  
  # ---- 4) hit cap without success ----
  if (!keep_grid) {
    return(list(found = FALSE,
                message = sprintf("Target power not reached up to nevents cap=%d.", ne_cap)))
  }
  
  best_i <- which.max(pow)
  ne_best <- grid[best_i]
  dum_best <- power_cached(ne_best)
  
  list(
    found = FALSE,
    nevents = ne_best,
    av_events = dum_best$av_events,    
    power = dum_best$pow,
    looks = looks_from_nevents(ne_best),
    Nvacc_mean = ceiling(mean(dum_best$vaccN,  na.rm = TRUE)),
    sero_pos_vaccN_mean = ceiling(mean(dum_best$sero_pos_vaccN,  na.rm = TRUE)),
    Nvacc_median = ceiling(median(dum_best$vaccN,  na.rm = TRUE)),
    sero_pos_vaccN_median = ceiling(median(dum_best$sero_pos_vaccN,  na.rm = TRUE)),    
    p_stop_1 = mean(dum_best$stop1, na.rm = TRUE),
    p_stop_2 = mean(dum_best$stop2, na.rm = TRUE),
    casesT = mean(dum_best$casesT, na.rm = TRUE),
    grid = grid,
    power_grid = pow,
    message = sprintf("Target power not reached up to nevents cap=%d.", ne_cap)
  )
}



find_nevents_for_power <- function(testof,cuts, target_power = 0.80,
                                   B = 3000,
                                   ne_min = 10L,
                                   ne_max = 40L,
                                   seed = 123L,
                                   down_factor = 0.67,   # shrink by ~1/1.5 each step
                                   up_factor   = 1.5,
                                   ne_floor = 1L,
                                   ...) {
  cat("\n=== Search target power =", target_power, "===\n")
  cat("Initial bounds:", ne_min, "to", ne_max, "\n")
  
  # --- ensure upper bound has >= target power ---
  dum   <- power_at_nevents(testof,cuts, ne_max, B = B, seed = seed, ...)
  p_max <- dum$pow
  cat(sprintf("Upper bound check: nevents=%d => power=%.3f\n", ne_max, p_max))
  
  while (p_max < target_power && ne_max < 1e6) {
    ne_max <- as.integer(ceiling(ne_max * up_factor))
    dum    <- power_at_nevents(testof,cuts, ne_max, B = B, seed = seed, ...)
    p_max  <- dum$pow
    cat(sprintf("Expanded upper bound: nevents=%d => power=%.3f\n", ne_max, p_max))
    flush.console()
  }
  if (p_max < target_power) stop("Could not find an upper bound achieving target power.")
  
  # --- ensure lower bound has < target power (expand downward if needed) ---
  dum   <- power_at_nevents(testof,cuts, ne_min, B = B, seed = seed, ...)
  p_min <- dum$pow
  cat(sprintf("Lower bound check: nevents=%d => power=%.3f\n", ne_min, p_min))
  
  while (p_min >= target_power && ne_min > ne_floor) {
    ne_min <- as.integer(max(ne_floor, floor(ne_min * down_factor)))
    dum    <- power_at_nevents(testof,cuts, ne_min, B = B, seed = seed, ...)
    p_min  <- dum$pow
    cat(sprintf("Shrunk lower bound: nevents=%d => power=%.3f\n", ne_min, p_min))
    flush.console()
    
    # protect against getting stuck if ne_min doesn't change due to integer rounding
    if (ne_min == ne_floor) break
  }
  
  # If even at floor we still exceed target, the minimum is the floor.
  if (p_min >= target_power && ne_min == ne_floor) {
    cat(sprintf("\n=== Result: nevents=%d already gives power=%.3f (at floor) ===\n",
                ne_min, p_min))
    return(list(
      nevents = ne_min,
      power   = p_min,
      looks   = looks_from_nevents(ne_min),
      Nvacc   = ceiling(median(dum$vaccN, na.rm = TRUE)),
      p_stop_1 = mean(dum$stop1, na.rm = TRUE),
      p_stop_2 = mean(dum$stop2, na.rm = TRUE),
      casesT  = mean(dum$casesT, na.rm = TRUE)
    ))
  }
  
  # Now we have a bracket: p(ne_min) < target <= p(ne_max)
  lo <- ne_min
  hi <- ne_max
  
  cat("\n--- Bisection search ---\n")
  while ((hi - lo) > 1L) {
    mid    <- as.integer(floor((lo + hi) / 2))
    dum_m  <- power_at_nevents(testof,cuts, mid, B = B, seed = seed, ...)
    p_mid  <- dum_m$pow
    cat(sprintf("Testing mid=%3d => power=%.3f (range %d-%d)\n", mid, p_mid, lo, hi))
    flush.console()
    if (p_mid >= target_power) {
      hi  <- mid
      dum <- dum_m   # keep last successful-ish result around
    } else {
      lo <- mid
    }
  }
  
  # Evaluate at hi for final report
  dum_f <- power_at_nevents(testof,cuts, hi, B = B, seed = seed, ...)
  final_power <- dum_f$pow
  cat(sprintf("\n=== Result: nevents=%d gives power=%.3f ===\n", hi, final_power))
  
  list(
    nevents = hi,
    power   = final_power,
    looks   = looks_from_nevents(hi),
    Nvacc   = ceiling(median(dum_f$vaccN, na.rm = TRUE)),
    p_stop_1 = mean(dum_f$stop1, na.rm = TRUE),
    p_stop_2 = mean(dum_f$stop2, na.rm = TRUE),
    casesT  = mean(dum_f$casesT, na.rm = TRUE)
  )
}



power_at_nevents <- function(testof, cuts, nevents, B = 3000, seed = 1L, ...) {
  set.seed(seed)
  looks_at <- looks_from_nevents(nevents)
  
  # ---- v4 binsplit: build exact discrete cutoffs for THIS nevents ----
  dots <- list(...)
  HR0_here <- if (!is.null(dots$HR0)) dots$HR0 else 0.80
  VE0_here <- 1 - HR0_here
  
  cuts_use <- cuts
  if (identical(testof, "binsplit")) {
    # Exact discrete group-sequential boundaries (strong alpha control)
    cuts_use <- exact_binom_cutoffs_from_nevents(
      nevents      = nevents,
      VE0          = VE0_here,
      n_v          = 1, n_c = 1,        # planned 1:1; ratio-only matters
      alpha_total  = 0.025,
      verbose      = FALSE
    )
  }
  
  rejections <- logical(B)
  stop1      <- logical(B)
  stop2      <- logical(B)
  vaccineN   <- numeric(B)
  sero_pos_vaccineN <- numeric(B) # *** number of sero-pos subjects in vaccine arm
  casesT     <- logical(B)  # indicator: vaccine-arm seroneg cases == 0
  
  for (b in seq_len(B)) {
    s <- sim_mult_wave8_plus(
      testof, cuts_use,
      nevents  = nevents,
      looks_at = looks_at,
      seedy    = sample.int(.Machine$integer.max, 1),
      ...
    )
    
    rejections[b] <- isTRUE(s$reject_one_sided)
    vaccineN[b]   <- s$total_vaccinated
    sero_pos_vaccineN[b] <- s$sero_pos_subs_vaccine # *** number of sero-pos subjects in vaccine arm  
    casesT[b]     <- (s$sero_neg_cases_vaccine == 0)
    
    if (isTRUE(s$stopped_early)) {
      if (s$stop_look == 1L) stop1[b] <- TRUE
      if (s$stop_look == 2L) stop2[b] <- TRUE
    }
    
    if (b %% 100 == 0) {
      cat(sprintf("  Rep %4d/%4d for nevents=%d for Nsero - vaccinations=%d\n", b, B, nevents,s$total_vaccinated)); flush.console()
    }
  }
  # Calculate expected number of events
  stop_look <- ifelse(stop1, 1L, ifelse(stop2, 2L, 3L))
  av_events = ceiling(mean(looks_at[stop_look]))
  
  list(
    av_events = av_events,
    pow   = mean(rejections, na.rm = TRUE),
    vaccN = vaccineN,
    sero_pos_vaccN = sero_pos_vaccineN, # *** number of sero-pos subjects in vaccine arm    
    stop1 = stop1,
    stop2 = stop2,
    casesT = casesT
  )
}


sim_mult_wave8_plus = function(
    testof, #"binsplit" or "poisson" 
    cuts,
    nevents = 20,
    looks_at = c(10, 15, 20),
    HR0 = 0.7,
    waves = rep(3,3),
    wave_attack = c(0.06, 0.06, 0.06),
    hr_true = 0.15,
    seroprev = 0.10,
    ## original per-subject frailty variance
    frailty_var = 0.5,
    ## OPTIONAL cluster frailty controls (default OFF)
    cluster_size = 5,
    rho_target = 0.3,
    cluster_frailty_var = 0,      # if 0 calcs else theta enables cluster frailty
    mapping = "exact", # second OK
    alpha_final = 0.025,          # one-sided alpha
    seedy = 7563685,
    max_enrolled = 1000000
    #testof =  "binsplit" #"poisson" 
) {
  set.seed(seedy)
  stopifnot(tail(looks_at, 1) == nevents)
  # 19 dec new version just uses totals for speed
  
  ## ===== Cluster frailty path (ONLY when requested) =====
  mapping = match.arg(mapping, choices = c("engineering","first","second","exact"))
  info_frac = looks_at / nevents
  # v4: cuts is either numeric z-bounds (poisson) or a list of exact cutoffs (binsplit)
  z_bound <- if (identical(testof, "poisson")) cuts else NULL
  Lambda_total = sum(-log(1 - wave_attack))
  
  theta_cluster = if (!is.null(cluster_frailty_var) && cluster_frailty_var > 0) {
    cluster_frailty_var
  } else {
    theta_from_rho(Lambda_total, rho_target, mapping)
    #0
  }
  
  #When events are extremely rare, most clusters have ~0 events, so the only way to still have a 
  # noticeable within-cluster correlation is if risk is wildly heterogeneous across clusters 
  #(some clusters effectively “hot”, most “cold”). 
  #In the Gamma-frailty model, that heterogeneity is controlled by Var(U)=θ. Tiny Λ ⇒ you need big θ to create enough clustering.
  
  theta_ind = max(0, frailty_var)
  
  enrolled = 0L; total_cases = 0L # total_cases is case_placebo + cases_vaccine (see below)
  cases_placebo = 0L; cases_vaccine = 0L # These are cases in the sero-neg group, others assume immun so no cases
  sero_neg_subs_placebo = 0L; sero_neg_subs_vaccine = 0L # These are just the number of sero-neg subjects enrolled
  sero_pos_subs_placebo = 0L; sero_pos_subs_vaccine = 0L # ditto for sero-pos
  # just uses totals for speed
  # time_all = numeric(max_enrolled); status_all = numeric(max_enrolled); group_all = numeric(max_enrolled)
  # cluster_all = integer(max_enrolled); immune_all = integer(max_enrolled)
  
  L = length(looks_at)
  z_seen = rep(NA_real_, L); ve_lower_seen = rep(NA_real_, L); hr_seen = rep(NA_real_, L)
  
  current_cluster_id = 0L; current_cluster_f = 1.0; within_left = 0L
  tiny_rate <- 1 / .Machine$double.xmax
  # --- Baseline seropositivity enriched in hotspots ---
  # try 0 (current model), 0.5, 1, 2
  b_sero=0.1
  if (b_sero == 0) {
    a_sero <- qlogis(seroprev)
    mu_logU <- 0
    sd_logU <- 1
  } else {
    Utmp <- rgamma(40000, shape = 1/theta_cluster, scale = theta_cluster)
    # CRITICAL: clamp before log to avoid -Inf / NaN moments
    logU <- log(clamp01e(Utmp, eps = 1e-12, cap = 1e12))
    mu_logU <- mean(logU)
    sd_logU <- sd(logU)
    
    if (!is.finite(sd_logU) || sd_logU <= 0) {
      stop("Bad sd_logU in sero calibration: sd_logU=", sd_logU,
           " theta_cluster=", theta_cluster)
    }
    a_sero <- calibrate_a_for_seroprev_scaled(seroprev, b_sero, logU, mu_logU, sd_logU)
  }
  
  while (total_cases < nevents & enrolled < max_enrolled) {
    
    if (within_left == 0L) {
      current_cluster_id = current_cluster_id + 1L
      current_cluster_f  = draw_gamma(theta_cluster)
      within_left        = cluster_size
    }
    
    enrolled   = enrolled + 1L
    within_left = within_left - 1L
    
    group  = rbinom(1, 1, 0.5)
    #    immune = rbinom(1, 1, seroprev) == 1
    if (b_sero == 0) {
      immune <- (runif(1) < seroprev)
    } else {
      # Guard cluster frailty
      if (!is.finite(current_cluster_f) || current_cluster_f <= 0) {
        current_cluster_f <- 1e-12
        # stop("current_cluster_f not finite/positive: ", current_cluster_f,
        #      " theta_cluster=", theta_cluster)
      }
      log_cf <- log(clamp01e(current_cluster_f, eps = 1e-12, cap = 1e12))
      zU <- (log_cf - mu_logU) / sd_logU
      eta <- a_sero + b_sero * zU
      
      # Force finite eta
      if (!is.finite(eta)) {
        stop("Non-finite eta in immune model: ",
             "eta=", eta,
             " a_sero=", a_sero,
             " b_sero=", b_sero,
             " zU=", zU,
             " current_cluster_f=", current_cluster_f,
             " mu_logU=", mu_logU,
             " sd_logU=", sd_logU,
             " theta_cluster=", theta_cluster)
      }
      
      p_immune <- plogis(eta)
      # p_immune <- pmin(pmax(p_immune, 0.02), 0.98)
      if (p_immune<0.02) p_immune=0.02
      if (p_immune>0.98) p_immune=0.98
      # Final guard
      if (!is.finite(p_immune) || p_immune < 0 || p_immune > 1) {
        stop("Non-finite p_immune: ",
             "p_immune=", p_immune,
             " eta=", eta,
             " theta_cluster=", theta_cluster)
      }
      
      immune <- (runif(1) < p_immune)
    }
    
    indiv_f = draw_gamma(theta_ind)
    frail_total = current_cluster_f * indiv_f
    
    t_ind = 0; event_occured = FALSE
    for (i in seq_along(waves)) {
      lambda0 = -log(1 - wave_attack[i]) / waves[i]
      lambda  = lambda0 * frail_total
      if (group == 1) lambda = lambda * hr_true
      if (immune == 1) lambda = 0
      if (lambda <= tiny_rate) {
        t_wave <- 6.408138e+15       # cannot happen in this wave
      } else {
        t_wave <- rexp(1, rate = lambda)
      }
      # if (is.na(t_wave)){
      #   df = 45
      # }
      #t_wave = rexp(1, rate = pmax(lambda, .Machine$double.eps))
      if (t_wave <= waves[i]) { t_ind = t_ind + t_wave; event_occured = TRUE; break }
      t_ind = t_ind + waves[i]
    }
    
    # time_all[enrolled]    <- t_ind
    # status_all[enrolled]  <- as.integer(event_occured)
    # group_all[enrolled]   <- group
    # cluster_all[enrolled] <- current_cluster_id
    # immune_all[enrolled]  <- as.integer(immune)
    
    # time_all   = c(time_all, t_ind)
    # status_all = c(status_all, as.numeric(event_occured))
    # group_all  = c(group_all, group)
    # cluster_all = c(cluster_all, current_cluster_id)
    # immune_all = c(immune_all, immune) 
    
    # Nite here that total_cases is sero-neg cases
    if (event_occured & immune==0) {
      total_cases = total_cases + 1L
      if (group == 0) cases_placebo = cases_placebo + 1L else cases_vaccine = cases_vaccine + 1L
    }
    if (immune==0){
      if (group == 0) sero_neg_subs_placebo = sero_neg_subs_placebo + 1L else sero_neg_subs_vaccine = sero_neg_subs_vaccine + 1L
    } else {
      if (group == 0) sero_pos_subs_placebo = sero_pos_subs_placebo + 1L else sero_pos_subs_vaccine = sero_pos_subs_vaccine + 1L
    }
    
    if (total_cases %in% looks_at) {
      k = match(total_cases, looks_at)
      ## Use robust SEs by cluster; if it fails (e.g., separation), Firth fallback
      #df = data.frame(time_all, status_all, group_all, cluster_all,immune_all)
      #df = df[1:enrolled, ]
      #df =df[df$immune_all==0,] # remove immune from analysis
      #events_by_arm = tapply(df$status_all, df$group_all, sum)
      #e_v <- events_by_arm["1"]; e_c <- events_by_arm["0"]
      #n_v <- sum(df$group_all == 1); n_c <- sum(df$group_all == 0)
      e_v <- cases_vaccine;e_c <- cases_placebo
      n_v <- sero_neg_subs_vaccine;n_c <- sero_neg_subs_placebo
      if (testof =="poisson"){
        outy = zou_hc0_logrr_se_totals(e_v, n_v, e_c, n_c, ha = TRUE)
        loghr_hat_look = outy$logrr
        se_loghr_look  = outy$se
        ok <- is.finite(loghr_hat_look) && is.finite(se_loghr_look) && (se_loghr_look > 0)
        if (!ok) {
          z_stat_look   <- NA_real_
          U_hr_one_look <- NA_real_
          ve_lower_look <- NA_real_
          hr_hat_look   <- NA_real_
        } else {
          z_stat_look   <- (loghr_hat_look - log(HR0)) / se_loghr_look
          U_hr_one_look <- exp(loghr_hat_look + z_bound[k] * se_loghr_look)
          ve_lower_look <- 1 - U_hr_one_look
          hr_hat_look   <- exp(loghr_hat_look)
        }
        z_seen[k]        = z_stat_look
        ve_lower_seen[k] = ve_lower_look
        hr_seen[k]       = hr_hat_look
      } else if ( testof =="binsplit") {
        # v4: exact discrete group-sequential case-split monitoring.
        # At each look k we have:
        #   N = e_v + e_c total seronegative endpoint cases observed
        #   v = e_v vaccine-arm cases among those N cases
        #
        # The group-sequential stopping rule is now *exact* and discrete:
        #   reject at look k if v <= c_k,
        # where c_k is pre-calibrated (for this nevents) using the exact binomial
        # distribution under H0 and an OBF-style alpha spending plan (strong control).
        #
        # 'cuts' is the list returned by exact_binom_cutoffs_from_nevents().
        
        c_k <- cuts$cutoffs[k]
        # Use the achieved cumulative alpha at this look (conservative due to discreteness)
        # for the supportive one-sided CP-style VE lower bound at this look.
        alpha_k <- cuts$alpha_used[k]
        
        outy <- ve_case_split2(v = e_v, N = e_v + e_c, n_v = n_v, n_c = n_c,
                               VE0 = 1 - HR0, alpha_k = alpha_k, c_k = c_k)
        
        # No Z-statistic in v4 for binsplit (we stop by integer cutoff on v)
        z_stat_look <- NA_real_
        z_seen[k] <- NA_real_
        ve_lower_seen[k] <- outy$VE_lower_onesided
        # Avoid divide-by-zero / Inf when an arm has zero seronegatives
        if (isTRUE(outy$ok) && n_v > 0 && n_c > 0 && e_c > 0) {
          hr_seen[k] <- (e_v / n_v) / (e_c / n_c)
        } else {
          hr_seen[k] <- NA_real_
        }
        #hr_seen[k] <- (e_v / n_v) / (e_c / n_c)
        hr_hat_look <- hr_seen[k]
        ve_lower_look <- ve_lower_seen[k]
      } 
      # Unified stopping decision for both tests
      reject_now <- if (identical(testof, "poisson")) {
        isTRUE(z_stat_look < -z_bound[k])
      } else {
        isTRUE(e_v <= cuts$cutoffs[k])
        isTRUE(outy$reject_H0)
      }
      
      if (k < length(looks_at) && isTRUE(reject_now)) {
        # trim
        return(list(
          reject_one_sided = TRUE,
          stopped_early    = TRUE,
          stop_look        = k,
          z_boundary       = if (identical(testof, "poisson")) z_bound[k] else cuts$cutoffs[k],
          z_stat           = z_stat_look,
          ve_hat           = 1 - hr_hat_look,
          ve_lower         = ve_lower_look,
          hr_hat           = hr_hat_look,
          z_seen           = z_seen,
          ve_lower_seen    = ve_lower_seen,
          hr_seen          = hr_seen,
          total_cases      = total_cases,
          sero_neg_cases_placebo    = cases_placebo,
          sero_neg_cases_vaccine    = cases_vaccine,
          total_enrolled   = enrolled,
          clusters_enrolled= current_cluster_id,
          sero_neg_subs_placebo = sero_neg_subs_placebo,
          sero_neg_subs_vaccine = sero_neg_subs_vaccine,
          sero_pos_subs_placebo = sero_pos_subs_placebo,
          sero_pos_subs_vaccine = sero_pos_subs_vaccine,
          total_vaccinated = sero_neg_subs_vaccine, # for consistency with prvious version Number sro neg vacc for rvf
          Lambda_total     = Lambda_total,
          theta_cluster_used = theta_cluster,
          theta_ind_used     = theta_ind,
          error            = FALSE
        ))
      } else if (k == length(looks_at) ){
        return(list(
          reject_one_sided = reject_now,
          stopped_early    = FALSE,
          stop_look        = k,
          z_boundary       = if (identical(testof, "poisson")) z_bound[k] else cuts$cutoffs[k],
          z_stat           = z_stat_look,
          ve_hat           = 1 - hr_hat_look,
          ve_lower         = ve_lower_look,
          hr_hat           = hr_hat_look,
          z_seen           = z_seen,
          ve_lower_seen    = ve_lower_seen,
          hr_seen          = hr_seen,
          total_cases      = total_cases,
          sero_neg_cases_placebo    = cases_placebo,
          sero_neg_cases_vaccine    = cases_vaccine,
          total_enrolled   = enrolled,
          clusters_enrolled= current_cluster_id,
          sero_neg_subs_placebo = sero_neg_subs_placebo,
          sero_neg_subs_vaccine = sero_neg_subs_vaccine,
          sero_pos_subs_placebo = sero_pos_subs_placebo,
          sero_pos_subs_vaccine = sero_pos_subs_vaccine,
          total_vaccinated = sero_neg_subs_vaccine, # for consistency with prvious version Number sro neg vacc for rvf
          Lambda_total     = Lambda_total,
          theta_cluster_used = theta_cluster,
          theta_ind_used     = theta_ind,
          error            = FALSE
        ))
      }
    }
  }
  
  if (total_cases < nevents) {
    return(list(
      error = TRUE,
      reject_one_sided = NA,
      stopped_early    = NA,
      stop_look        = NA_integer_,
      total_cases      = total_cases,
      total_enrolled   = enrolled,
      clusters_enrolled= current_cluster_id,
      sero_neg_cases_placebo = cases_placebo,
      sero_neg_cases_vaccine = cases_vaccine,
      sero_neg_subs_placebo  = sero_neg_subs_placebo,
      sero_neg_subs_vaccine  = sero_neg_subs_vaccine,
      sero_pos_subs_placebo  = sero_pos_subs_placebo,
      sero_pos_subs_vaccine  = sero_pos_subs_vaccine,
      total_vaccinated = sero_neg_subs_vaccine,  # keep expected field name
      Lambda_total     = Lambda_total,
      theta_cluster_used = theta_cluster,
      theta_ind_used     = theta_ind
    ))
  }
  
}



theta_from_rho = function(Lambda, rho, method) {
  # best to use exact
  stopifnot(rho > 0, rho < 1, Lambda > 0)
  if (method == "first")       return(rho / Lambda) # crapola
  if (method == "second")      return(rho * (1 - Lambda) / Lambda)
  if (method == "engineering") return(rho / (Lambda * (1 - rho))) #crapola
  f = function(theta){
    a = (1 + theta * Lambda)^(-1/theta)
    b = (1 + 2 * theta * Lambda)^(-1/theta)
    (b - a^2) / (a * (1 - a)) - rho
  }
  uniroot(f, lower = 1e-9, upper = 1000)$root
}

#draw_gamma = function(theta) if (theta > 0) rgamma(1, shape = 1/theta, scale = theta) else 1.0
draw_gamma <- function(theta, eps = 1e-12, cap = 1e12) {
  # Robust gamma frailty draw with clamping to avoid 0/Inf/NaN under extreme theta.
  if (!is.finite(theta) || theta <= 0) return(1.0)
  
  u <- rgamma(1, shape = 1/theta, scale = theta)
  
  # Handle underflow/overflow/non-finite
  if (!is.finite(u) || u <= 0) u <- eps
  
  # Clamp to a reasonable range for stability
  # u <- min(max(u, eps), cap)
  if (u<eps) u=eps
  if (u>cap) u=cap
  u
}

zou_hc0_logrr_se_totals <- function(E_v, N_v, E_c, N_c, ha = TRUE) {
  stopifnot(N_v > 0, N_c > 0)
  
  if (ha && (E_v == 0L || E_c == 0L)) {
    e_v <- E_v + 0.5; n_v <- N_v + 1
    e_c <- E_c + 0.5; n_c <- N_c + 1
  } else {
    e_v <- E_v; n_v <- N_v
    e_c <- E_c; n_c <- N_c
  }
  
  p_v <- e_v / n_v
  p_c <- e_c / n_c
  logrr <- log(p_v / p_c)
  
  se =sqrt (1/e_v - 1/n_v + 1/e_c - 1/n_c)
  
  # # HC0 meat pieces: sum (y - mu)^2 within each arm
  # A_c <- e_c * (1 - p_c)^2 + (n_c - e_c) * (p_c)^2
  # A_v <- e_v * (1 - p_v)^2 + (n_v - e_v) * (p_v)^2
  # 
  # # bread X'WX with W = mu (Poisson log-link); sums of mu are e_c, e_v
  # S0 <- e_c + e_v
  # S1 <- e_v
  # XWX <- matrix(c(S0, S1,
  #                 S1, S1), 2, 2, byrow = TRUE)
  # 
  # Meat <- matrix(c(A_c + A_v, A_v,
  #                  A_v,       A_v), 2, 2, byrow = TRUE)
  # 
  # invXWX <- solve(XWX)
  # V <- invXWX %*% Meat %*% invXWX
  
  list(logrr = logrr, se = se ) #sqrt(V[2,2]))
}

# ============================================================
# v4: Exact discrete OBF-style boundaries for the BIN-SPLIT test
# ============================================================

# One-sided OBF alpha spending (cumulative alpha at information fraction t)
alpha_spend_obf_one_sided <- function(alpha_total, t) {
  stopifnot(alpha_total > 0, alpha_total < 1, all(t > 0), all(t <= 1))
  z_alpha <- qnorm(1 - alpha_total)
  1 - pnorm(z_alpha / sqrt(t))
}

# Null theta0 from VE0 and planned allocation counts (ratio matters, not scale)
theta0_from_VE0_nvnc <- function(VE0, n_v, n_c) {
  stopifnot(VE0 >= 0, VE0 < 1, n_v > 0, n_c > 0)
  (n_v * (1 - VE0)) / (n_v * (1 - VE0) + n_c)
}

# Convolution update: add Binomial(dN, theta0) increment to distribution of V
# (See detailed comments in your v4 draft; kept compact here.)
convolve_with_binom_increment <- function(p_prev, dN, theta0) {
  stopifnot(dN >= 0, theta0 > 0, theta0 < 1)
  if (dN == 0) return(p_prev)
  inc <- dbinom(0:dN, size = dN, prob = theta0)
  as.numeric(convolve(p_prev, rev(inc), type = "open"))
}

# Compute exact discrete cutoffs c_k for looks at 50%, 75%, 100% of nevents.
# Returns a list with:
#   cutoffs[k]   integer c_k (reject at look k if V_k <= c_k)
#   alpha_used[k] achieved cumulative alpha by look k (<= OBF targets; strong control)
exact_binom_cutoffs_from_nevents <- function(
    nevents,
    VE0,
    n_v = 1, n_c = 1,          # planned ratio (default 1:1)
    alpha_total = 0.025,
    verbose = FALSE
) {
  stopifnot(nevents >= 3L, alpha_total > 0, alpha_total < 1)
  nevents <- as.integer(nevents)
  
  Nlooks <- looks_from_nevents(nevents)
  K <- length(Nlooks)
  t <- Nlooks / Nlooks[K]
  
  theta0 <- theta0_from_VE0_nvnc(VE0, n_v, n_c)
  A <- alpha_spend_obf_one_sided(alpha_total, t)  # OBF cumulative targets
  
  p_notrej <- 1.0
  alpha_used <- 0.0
  
  cutoffs <- integer(K)
  alpha_used_by_look <- numeric(K)
  
  Nprev <- 0L
  for (k in seq_len(K)) {
    Nk <- Nlooks[k]
    dN <- Nk - Nprev
    
    p_notrej <- convolve_with_binom_increment(p_notrej, dN, theta0)
    cdf <- cumsum(p_notrej)
    
    remaining <- A[k] - alpha_used
    
    if (remaining <= 0) {
      ck <- -1L
      spend <- 0.0
    } else {
      ok <- which(cdf <= remaining)
      if (length(ok) == 0) {
        ck <- -1L
        spend <- 0.0
      } else {
        ck <- max(ok) - 1L
        spend <- cdf[ck + 1L]
      }
    }
    
    cutoffs[k] <- ck
    alpha_used <- alpha_used + spend
    alpha_used_by_look[k] <- alpha_used
    
    if (ck >= 0) p_notrej[1:(ck + 1L)] <- 0
    Nprev <- Nk
  }
  
  if (verbose) {
    cat("\nExact binomial cutoffs:\n")
    cat("  Nlooks:", paste(Nlooks, collapse = ","), "\n")
    cat("  cutoffs:", paste(cutoffs, collapse = ","), "\n")
    cat("  alpha_targets:", paste(sprintf("%.6f", A), collapse = ","), "\n")
    cat("  alpha_used:", paste(sprintf("%.6f", alpha_used_by_look), collapse = ","), "\n")
  }
  
  list(
    Nlooks = Nlooks,
    t = t,
    VE0 = VE0,
    theta0 = theta0,
    alpha_total = alpha_total,
    alpha_targets = A,
    cutoffs = cutoffs,
    alpha_used = alpha_used_by_look
  )
}



# Case-split (conditional) VE inference under unequal randomization
# V | N ~ Binomial(N, theta),   theta(VE) = (n_v*(1-VE)) / (n_v*(1-VE) + n_c)
#
# IMPORTANT: Two different "p-values"/decisions can be reported:
#   1) p_value      = exact one-sided binomial p-value at theta0 (discrete; exact)
#   2) p_value_z    = Normal/score p-value based on zpos (approx; used for OBF looks)
#
# In your main code you pass alpha = alpha_k = 1 - Phi(z_bound[k]),
# so the sequential stopping rule should be based on the score statistic:
#    reject at look k  <=>  zpos >= z_bound[k]  <=>  (1 - Phi(zpos)) <= alpha_k
#

# ============================================================
# v4: Exact case-split inference at a look-specific alpha_k
# ============================================================
#
# Inputs:
#   v       = vaccine cases observed among the N total endpoint cases (sero-neg only)
#   N       = total endpoint cases observed (vaccine + control)
#   n_v,n_c = numbers randomized (sero-neg) to vaccine/control (used for theta0 and VE inversion)
#   VE0     = null threshold (H0: VE <= VE0)
#   alpha_k = one-sided nominal level for this look (use achieved cumulative alpha from
#             exact_binom_cutoffs_from_nevents(), i.e. cuts$alpha_used[k])
#   c_k     = exact discrete rejection cutoff for this look (reject if v <= c_k)
#
# Returns:
#   p_value              = exact one-sided binomial p-value under H0 (theta=theta0)
#   VE_lower_onesided    = one-sided (1-alpha_k) lower bound for VE from CP-type inversion
#   reject_H0            = decision consistent with the discrete cutoff (v <= c_k)
#
ve_case_split2 <- function(v, N, n_v, n_c, VE0 = 0.20, alpha_k = 0.025, c_k) {
  
  # Basic sanity on v,N only (these should always hold)
  if (!(N >= 1 && v >= 0 && v <= N)) {
    return(list(
      theta0 = NA_real_,
      p_value = NA_real_,
      c_k = as.integer(c_k),
      reject_H0 = FALSE,
      thetaU = NA_real_,
      VE_lower_onesided = NA_real_,
      ok = FALSE,
      msg = "Invalid (v,N) inputs"
    ))
  }
  
  # If one arm has zero seronegatives so far, the case-split test is not defined.
  # Do NOT crash the simulation; treat as "cannot reject at this look".
  if (is.na(n_v) || is.na(n_c) || n_v <= 0 || n_c <= 0) {
    return(list(
      theta0 = NA_real_,
      p_value = NA_real_,
      c_k = as.integer(c_k),
      reject_H0 = FALSE,
      thetaU = NA_real_,
      VE_lower_onesided = NA_real_,
      ok = FALSE,
      msg = sprintf("Case-split undefined: n_v=%s n_c=%s", n_v, n_c)
    ))
  }
  
  # --- normal path below ---
  c_k <- as.integer(c_k)
  
  theta0 <- (n_v * (1 - VE0)) / (n_v * (1 - VE0) + n_c)
  p_value <- pbinom(v, size = N, prob = theta0)
  reject_H0 <- (v <= c_k)
  
  if (v == N) {
    thetaU <- 1
    VE_L <- -Inf
  } else {
    thetaU <- qbeta(1 - alpha_k, shape1 = v + 1, shape2 = N - v)
    VE_L <- 1 - (thetaU * n_c) / (n_v * (1 - thetaU))
  }
  
  list(
    theta0 = theta0,
    p_value = p_value,
    c_k = c_k,
    reject_H0 = reject_H0,
    thetaU = thetaU,
    VE_lower_onesided = VE_L,
    ok = TRUE,
    msg = NA_character_
  )
}
looks_from_nevents <- function(nevents) {
  stopifnot(nevents >= 3L)
  a <- max(1L, floor(0.50 * nevents))
  b <- max(2L, floor(0.75 * nevents))
  l <- sort(unique(c(a, b)))
  l <- l[l < nevents]
  c(l, nevents)
}
# ---- Example ----
# res <- ve_case_split(v = 28, N = 100, VE0 = 0.30, alpha = 0.025)
# print(res)
getOBF_cuts = function(){
  library(gsDesign)
  
  alpha  <- 0.025
  timing <- c(0.5, 0.75, 1)
  
  des <- gsDesign(k=3, test.type=1, alpha=alpha, timing=timing, sfu="OF")
  
  b <- as.numeric(des$upper$bound)  # these are b1,b2,b3
  return(b)
  
}
calibrate_a_for_seroprev_scaled <- function(sero_target, b_sero,
                                            logU, mu_logU, sd_logU,
                                            maxit = 60L) {
  # ---------------------------------------------------------------------------
  # PURPOSE
  # ---------------------------------------------------------------------------
  # We want baseline seropositivity ("immune") to be:
  #   (i) correlated with cluster risk (frailty) when b_sero > 0, AND
  #  (ii) to have a prescribed *overall* marginal prevalence sero_target.
  #
  # In the simulation, each cluster has a frailty multiplier U > 0 (higher U = higher
  # infection hazard / "hotter" cluster). We model the probability of being immune
  # at baseline as a logistic function of (standardized) log cluster frailty:
  #
  #   p_immune(U) = logistic( a + b_sero * z ),
  #   where z = (logU - mu_logU) / sd_logU.
  #
  # Here:
  #   - b_sero controls how strongly baseline seropositivity increases with cluster
  #     risk (positive correlation when b_sero > 0).
  #   - a is an intercept that shifts probabilities up/down for *all* clusters.
  #
  # Crucially, changing b_sero changes the average seropositivity across the frailty
  # distribution. For example, if b_sero > 0, hot clusters get higher p_immune, so
  # the intercept a must shift downward to keep the *marginal* prevalence equal to
  # sero_target. This function finds the value of a that enforces:
  #
  #   E[ logistic(a + b_sero * z) ] = sero_target,
  #
  # where the expectation is over the distribution of cluster frailties. We
  # approximate that expectation by Monte Carlo using the supplied logU sample.
  #
  # OUTPUT
  #   Returns the calibrated intercept a_sero such that the simulated average immune
  #   probability matches sero_target (to within bisection tolerance).
  #
  # NUMERICAL / MODELING NOTES
  #   - We standardize logU to z so that b_sero is on a stable scale (per 1 SD change
  #     in log frailty), rather than depending on the raw spread of logU, which varies
  #     with theta_cluster (and hence rho_target and attack rate).
  #   - When theta_cluster is very large (heavy mass near 0), raw log(U) can include
  #     -Inf. Upstream we clamp U away from 0 before taking logs; nevertheless we still
  #     sanitize finite values here.
  #   - The mapping a -> mean(plogis(a + b*z)) is strictly increasing in a, so a unique
  #     solution exists and bisection is safe.
  # ---------------------------------------------------------------------------
  
  
  stopifnot(is.finite(sero_target), sero_target > 0, sero_target < 1)
  stopifnot(is.finite(b_sero))
  
  # sanitize
  logU <- logU[is.finite(logU)]
  if (length(logU) < 100L) stop("calibrate_a_for_seroprev_scaled: logU too short")
  
  if (!is.finite(mu_logU)) mu_logU <- mean(logU)
  if (!is.finite(sd_logU)) sd_logU <- sd(logU)
  
  if (!is.finite(sd_logU) || sd_logU <= 0) {
    # No variation in logU (shouldn't happen unless theta_cluster ~ 0 or Utmp degenerate)
    return(qlogis(sero_target))
  }
  
  z <- (logU - mu_logU) / sd_logU
  z <- z[is.finite(z)]
  if (length(z) < 100L) stop("calibrate_a_for_seroprev_scaled: z has too few finite values")
  
  # Stable computation of mean(plogis(a + b*z))
  mean_plogis <- function(a) {
    eta <- a + b_sero * z
    # plogis is stable for large magnitude eta; still guard:
    p <- plogis(eta)
    mean(p, na.rm = TRUE)
  }
  
  lo <- -50; hi <- 50
  p_lo <- mean_plogis(lo)
  p_hi <- mean_plogis(hi)
  
  if (!is.finite(p_lo) || !is.finite(p_hi)) {
    stop("calibrate_a_for_seroprev_scaled: p_lo/p_hi not finite. ",
         "Check sd_logU=", sd_logU, " range(z)=", paste(range(z), collapse = ","))
  }
  
  # Ensure bracket is valid: plogis(lo) ~ 0, plogis(hi) ~ 1, so target must lie between.
  if (!(p_lo <= sero_target && sero_target <= p_hi)) {
    stop("calibrate_a_for_seroprev_scaled: target not bracketed. ",
         "p_lo=", p_lo, " p_hi=", p_hi, " target=", sero_target)
  }
  
  for (i in seq_len(maxit)) {
    mid <- (lo + hi) / 2
    p_mid <- mean_plogis(mid)
    if (!is.finite(p_mid)) {
      stop("calibrate_a_for_seroprev_scaled: p_mid not finite at iter ", i,
           ". sd_logU=", sd_logU, " mid=", mid)
    }
    if (p_mid > sero_target) hi <- mid else lo <- mid
  }
  
  (lo + hi) / 2
}

clamp01e <- function(x, eps = 1e-12, cap = 1e12) {
  # clamp to [eps, cap] and keep finite
  x <- ifelse(is.finite(x), x, NA_real_)
  x[x<eps]=eps
  x[x>cap]=cap
  # x <- pmax(x, eps)
  # x <- pmin(x, cap)
  x
}


run_test_power9HPC()
