# ============================================================================
# SECTION 0: Reformat Contact Matrix
# ============================================================================
source("reconstruct_contact_matrix_M.R")  # provides M_infant_room, M_household
if (!exists("RUN_ANALYSIS")) RUN_ANALYSIS <- TRUE

# ============================================================================
# SECTION 1: PARAMETERS
# ============================================================================
set.seed(2026178815)
study_start <- "2024-07-01"
study_end   <- "2025-06-30"
sim_dates   <- seq(as.Date(study_start), as.Date(study_end), by = "day")

peak_day <- as.integer(as.Date("2025-01-01") - as.Date(study_start)) + 1L   # ~day 185
sd_days  <- 40
seasonal <- exp(-((seq_len(length(sim_dates)) - peak_day)^2) / (2 * sd_days^2))

surveillance_data <- data.frame(
  date  = sim_dates,
  cases = 0.2 + 100 * seasonal + abs(rnorm(length(sim_dates), 0, 5))
)
plot(surveillance_data$date, surveillance_data$cases, type = "l")

params <- list(
  pop_seed = 20261234L,
  n_sims   = 50,
  # ── DAYCARE COMPOSITION ─────────────────────────────────────────────────
  n_infants_dc      = 90L,
  n_toddlers_dc     = 90L,
  staff_child_ratio = 3L,
  # ── CLUSTER (classroom) STRUCTURE ────────────────────────────────────────
  cluster_size          = 9L,
  cluster_within_weight = 1.00,
  cluster_between_weight_child = 0,
  cluster_between_weight_staff = 1,
  # ── DAYCARE TRANSMISSION SWITCH ──────────────────────────────────────────
  sick_attendance  = 0.7,
  attendance_child = 1.00,
  attendance_staff = 1.00,
  daycare_entry_age_days = 120L,
  daycare_bring_home_sameday = TRUE,
  # ── HOUSEHOLD STRUCTURE ──────────────────────────────────────────────────
  p_infant        = c(0, 1),
  p_toddler       = c(0.7, 0.25, 0.05),
  p_elderly       = c(0.92, 0.07, 0.01),
  p_toddler_only  = c(0.8, 0.2),
  n_parents       = 2L,
  p_two_caregiver = 0.85,
  staff_hh_size = 2L,
  room_mix_cutoff_days  = 540L,
  toddler_age_beta_a = 2,   # Beta shape a  (toddler age skew)
  toddler_age_beta_b = 8,   # Beta shape b  
  # ── VACCINATION (nirsevimab, infants only) ───────────────────────────────
  rsv_season_start_day = 93L,
  rsv_season_months    = c(10L,11L,12L,1L,2L,3L),
  mab_age_max_days     = 240L,
  # ── OUTCOME WINDOW ───────────────────────────────────────────────────────
  outcome_from_day = 93L,
  # ── Age ranges (days) ────────────────────────────────────────────────────
  infant_age_min = 1L,   infant_age_max = 121L,
  toddler_age_min= 366L,  toddler_age_max= 1461L,
  adult_age_min  = 9000L, adult_age_max  = 14000L,
  elderly_age_min= 23000L,elderly_age_max= 29000L,
  staff_age      = 10000L,
  # ── MSIS parameters ──────────────────────────────────────────────────────
  omega_M       = 1 / 112,   # maternal-antibody waning (M -> S)
  omega_Mb      = 1 / 168.6,   
  omega_mab_fix = 1 / 53.6,    
  rho_mab       = 0.2,
  sigma_1 = 0.76, sigma_2 = 0.6, sigma_3 = 0.4,
  rho1 = 0.75, rho2 = 0.51,
  rho  = 1, inf_shape = 3,
  # ── Role-dependent susceptibility (phi) & infectivity (kappa) ────────────
  phi_by_role   = c(infant=3.27, toddler=1.98, adult=1.00, elderly=1.00, staff=1.00),
  kappa_by_role = c(infant=0.78, toddler=1.69, adult=1.00, elderly=1.00, staff=1.00),
  # ── Ct viral dynamics ──────────────────────────────────────────────────────
  params_ct_values = list(
    infant  = list(Cpeak=33.72, r=2.33, d=0.53, t_peak=2.72),
    toddler = list(Cpeak=34.67, r=1.79, d=0.30, t_peak=2.94),
    adult   = list(Cpeak=34.34, r=1.75, d=0.19, t_peak=3.01),
    elderly = list(Cpeak=34.34, r=1.75, d=0.19, t_peak=3.01),
    staff   = list(Cpeak=34.34, r=1.75, d=0.19, t_peak=3.01)
  ),
  vl_midpoint = 32.79, vl_slope = 1.92,
  Ct_LOD = 45, detect_threshold_Ct = 45,
  # ── Transmission parameters ────────────────────────────────────────────────
  beta1 = .02, beta2 = .02,
  C_daycare   = M_infant_room,
  C_household = M_household,
  # ── Community FOI / scaling ─────────────────────────────────────────────────
  lambda_ext_base = 5e-4,
  beta_daycare    = 1.5,
  beta_household  = 1.5,
  daycare_hours   = 1,
  household_hours = 1,
  role_init = list(
    infant  = list(n_init = 0L, p_immune = 0.00),
    toddler = list(n_init = 0L, p_immune = 0.00),
    adult   = list(n_init = 0L, p_immune = 0.00),
    elderly = list(n_init = 0L, p_immune = 0.00),
    staff   = list(n_init = 0L, p_immune = 0.00)
  )
)
# Derived
params$study_start_date <- as.Date(study_start)
params$n_days           <- length(surveillance_data$date)
params$surveillance_data<- surveillance_data$cases

# Remove the spread_even function, .K calculations, and the cat() statement.
params$N_infants  <- params$n_infants_dc
params$N_children <- params$n_infants_dc + params$n_toddlers_dc


params$N_infants     <- params$n_infants_dc
params$N_children    <- params$n_infants_dc + params$n_toddlers_dc



# ============================================================================
# SECTION 2: HELPER FUNCTIONS
# ============================================================================
get_sigma <- function(n_prior, params) {
  if (n_prior == 0) return(1.0)
  if (n_prior == 1) return(params$sigma_1)
  if (n_prior == 2) return(params$sigma_2)
  params$sigma_3
}
get_inf_duration <- function(n_inf, params) {
  mean_dur <- if (n_inf == 1) 10 else if (n_inf == 2) 7 else 5
  dur <- stats::rgamma(1, shape = params$inf_shape, scale = mean_dur / params$inf_shape)
  pmax(1L, ceiling(dur))
}
get_rho <- function(n_inf, params)
  if (n_inf == 1) 1.0 else if (n_inf == 2) params$rho1 else params$rho2

compute_Ct_direct <- function(tau, role, params) {
  ct_p <- params$params_ct_values[[role]]
  if (is.null(ct_p)) ct_p <- params$params_ct_values[["adult"]]
  if (tau < 0) return(params$Ct_LOD)
  ct <- if (tau <= ct_p$t_peak) ct_p$Cpeak + ct_p$r * (ct_p$t_peak - tau)
  else                    ct_p$Cpeak + ct_p$d * (tau - ct_p$t_peak)
  min(ct, params$Ct_LOD)
}
f_Ct <- function(ct, params) 1 / (1 + exp((ct - params$vl_midpoint) / params$vl_slope))
get_phi_role <- function(role, params) {
  val <- params$phi_by_role[role]
  if (is.na(val)) return(unname(params$phi_by_role["adult"]))
  unname(val)
}
get_kappa <- function(role, params) {
  val <- params$kappa_by_role[role]
  if (is.na(val)) return(unname(params$kappa_by_role["adult"]))
  unname(val)
}
lambda_ext_seasonal <- function(day, params) {
  daily_cases <- params$surveillance_data[day]
  params$lambda_ext_base * (daily_cases / max(params$surveillance_data))
}

if (!exists("INFECTOR_ATTRIB")) INFECTOR_ATTRIB <- "proportional"
pick_infector <- function(idx_vec, val_vec) {
  if (length(idx_vec) == 0L) return(NA_integer_)
  if (length(idx_vec) == 1L) return(idx_vec[1L])
  if (identical(INFECTOR_ATTRIB, "proportional"))
    idx_vec[sample.int(length(idx_vec), 1L, prob = val_vec)]
  else idx_vec[which.max(val_vec)]
}

# ============================================================================
# SECTION 3: GENERATE POPULATION  (daycare-first)
# ============================================================================
make_agent <- function(gid, hh, role, age, is_dc, cluster_id = NA_integer_) {
  data.frame(
    global_id = gid, hh_id = hh, role = role, age_days = age,
    is_daycare = is_dc, cluster_id = cluster_id,
    intended_vax = FALSE, state = "S", n_infections = 0L,
    days_in_state = 0L, ever_infected = FALSE, infected_outcome = FALSE,
    breakthrough = FALSE,                       # [MODEL A] in-window ASYMPTOMATIC breakthrough (S_ab -> I_ab)
    infection_source = NA_character_, day_infected = NA_integer_,
    recovery_day = NA_integer_, vax_due_day = NA_integer_,
    stringsAsFactors = FALSE
  )
}

generate_population <- function(params, quiet = TRUE) {
  old_seed <- .Random.seed
  set.seed(params$pop_seed)
  on.exit(.Random.seed <<- old_seed)
  
  Z <- params$staff_child_ratio
  X <- params$n_infants_dc; Y <- params$n_toddlers_dc
  if (X <= 0) stop("n_infants_dc must be > 0.")
  
  mem <- list(); gid <- 0L; hh <- 0L
  infant_gids <- integer(0); toddler_gids <- integer(0)
  add <- function(role, age, is_dc, cl = NA_integer_) {
    gid <<- gid + 1L
    mem[[gid]] <<- make_agent(gid, hh, role, age, is_dc, cl)
    gid
  }
  draw_count <- function(p) sample.int(length(p), 1L, prob = p) - 1L
  add_caregivers_and_elderly <- function() {
    n_par <- if (runif(1) < params$p_two_caregiver) 2L else 1L
    for (p in seq_len(n_par))
      add("adult", sample(params$adult_age_min:params$adult_age_max, 1L), FALSE)
    for (e in seq_len(draw_count(params$p_elderly)))
      add("elderly", sample(params$elderly_age_min:params$elderly_age_max, 1L), FALSE)
  }
  
  
  add_toddler <- function() {
    a <- if (!is.null(params$toddler_age_beta_a)) params$toddler_age_beta_a else 2
    b <- if (!is.null(params$toddler_age_beta_b)) params$toddler_age_beta_b else 4
    u <- rbeta(1, a, b)                                            # in (0,1), skewed young when b>a
    age <- params$toddler_age_min +
      round(u * (params$toddler_age_max - params$toddler_age_min))
    g <- add("toddler", as.integer(age), TRUE)
    toddler_gids <<- c(toddler_gids, g)
  }
  
  ## ---- BATCH 1: infant-anchored households --------------------------------
  for (h in seq_len(X)) {
    hh <- hh + 1L
    add_caregivers_and_elderly()
    g <- add("infant", sample(params$infant_age_min:params$infant_age_max, 1L), TRUE)
    infant_gids <- c(infant_gids, g)
    n_tod <- min(draw_count(params$p_toddler), Y - length(toddler_gids))
    for (t in seq_len(max(0L, n_tod))) add_toddler()
  }
  
  ## ---- BATCH 2: toddler-only households (top up to n_toddlers_dc) ----------
  p_to <- if (!is.null(params$p_toddler_only)) params$p_toddler_only else c(0.8, 0.2)
  guard <- 0L
  while (length(toddler_gids) < Y) {
    hh <- hh + 1L
    add_caregivers_and_elderly()
    n_tod <- min(sample.int(length(p_to), 1L, prob = p_to), Y - length(toddler_gids))
    for (t in seq_len(max(1L, n_tod))) add_toddler()
    guard <- guard + 1L; if (guard > Y + 5L) break
  }
  
  ## ---- ROOMS: age-banded, blended -----------------------------------------
  mix_cut <- if (!is.null(params$room_mix_cutoff_days)) params$room_mix_cutoff_days else 730L
  tod_ages <- if (length(toddler_gids)) sapply(toddler_gids, function(g) mem[[g]]$age_days) else numeric(0)
  young_tod <- toddler_gids[tod_ages <= mix_cut]
  older_tod <- toddler_gids[tod_ages >  mix_cut]
  
  # POOL infants + young toddlers, shuffle, chunk -> proportional blend per room
  young_pool <- sample(c(infant_gids, young_tod))
  older_pool <- if (length(older_tod)) sample(older_tod) else integer(0)
  
  K1 <- if (length(young_pool) > 0) ceiling(length(young_pool) / params$cluster_size) else 0L
  K2 <- if (length(older_pool) > 0) ceiling(length(older_pool) / params$cluster_size) else 0L
  K  <- K1 + K2
  clust_n <- rep(0L, max(K, 1L))
  
  if (K1 > 0) for (k in seq_along(young_pool)) {
    ci <- ((k - 1L) %% K1) + 1L
    mem[[young_pool[k]]]$cluster_id <- ci; clust_n[ci] <- clust_n[ci] + 1L
  }
  if (K2 > 0) for (k in seq_along(older_pool)) {
    ci <- ((k - 1L) %% K2) + 1L + K1
    mem[[older_pool[k]]]$cluster_id <- ci; clust_n[ci] <- clust_n[ci] + 1L
  }
  
  ## ---- staff per room -----------------------------------------------------
  actual_staff_count <- 0L
  for (ci in seq_len(K)) {
    staff_needed <- ceiling(clust_n[ci] / Z)
    actual_staff_count <- actual_staff_count + staff_needed
    for (t in seq_len(staff_needed)) {
      hh <- hh + 1L
      add("staff", params$staff_age, TRUE, ci)
      for (sp in seq_len(params$staff_hh_size - 1L))
        add("adult", sample(params$adult_age_min:params$adult_age_max, 1L), FALSE)
    }
  }
  
  if (!quiet) {
    n_young_tod <- length(young_tod)
    cat(sprintf(
      "Centre (2-batch, blended): %d infants + %d toddlers | young rooms K1=%d (pool %d = %d inf + %d young tod), older rooms K2=%d | %d staff\n",
      length(infant_gids), length(toddler_gids), K1, length(young_pool),
      length(infant_gids), n_young_tod, K2, actual_staff_count))
  }
  list(agents = do.call(rbind, mem), n_staff = actual_staff_count,
       n_clusters = K, n_clusters_g1 = K1)
}

# ============================================================================
# SECTION 4: INITIALIZE AGENTS  (nirsevimab logic)  [MODEL A staging]
# ============================================================================
initialize_agents <- function(params, vax_coverage = 0.0, quiet = TRUE) {
  #agents <- generate_population(params, quiet = quiet)
  pop    <- generate_population(params, quiet = quiet)
  agents <- pop$agents
  params$N_staff       <<- pop$n_staff
  params$n_dc_clusters <<- pop$n_clusters
  params$n_dc_clusters_g1 <<- pop$n_clusters_g1
  
  ss     <- params$study_start_date
  for (i in which(agents$role == "infant")) {
    age0 <- agents$age_days[i]
    agents$n_infections[i] <- 0L
    maternal_state <- function() if (runif(1) < exp(-age0 * params$omega_M)) "M" else "S"
    wants_mab  <- runif(1) < vax_coverage
    birth_date <- ss - age0
    bmonth     <- as.integer(format(birth_date, "%m"))
    in_season  <- bmonth %in% params$rsv_season_months
    if (in_season) {
      dose_date <- birth_date
    } else {
      by  <- as.integer(format(birth_date, "%Y"))
      o1  <- as.Date(sprintf("%d-10-01", by))
      dose_date <- if (birth_date <= o1) o1 else as.Date(sprintf("%d-10-01", by + 1L))
    }
    age_at_dose  <- as.numeric(dose_date - birth_date)
    eligible     <- age_at_dose <= params$mab_age_max_days
    dose_sim_day <- as.integer(dose_date - ss) + 1L
    if (wants_mab && eligible && dose_sim_day <= params$n_days) {
      agents$intended_vax[i] <- TRUE
      if (dose_sim_day < 1L) {
       
        elapsed <- 1L - dose_sim_day
        st <- "M_b"; d_in <- 0L
        # for (k in seq_len(elapsed)) {
        #   if (st == "M_b") {
        #     if (d_in > 0 && runif(1) < (1 - exp(-params$omega_mab_fix))) { st <- "S_ab"; d_in <- 0L }
        #   } else if (st == "S_ab") {
        #     if (runif(1) < (1 - exp(-params$omega_Mb))) { st <- "S"; d_in <- 0L }
        #   }
        #   d_in <- d_in + 1L
        # }
        for (k in seq_len(elapsed)) {
          d_in <- d_in + 1L        # move this to the TOP of the loop
          if (st == "M_b") {
            if (runif(1) < (1 - exp(-params$omega_mab_fix))) { st <- "S_ab"; d_in <- 0L }
          } else if (st == "S_ab") {
            if (runif(1) < (1 - exp(-params$omega_Mb))) { st <- "S"; d_in <- 0L }
          }
        }
        agents$state[i]         <- st
        agents$days_in_state[i] <- d_in
        agents$vax_due_day[i]   <- NA_integer_
      } else {
        agents$state[i]       <- maternal_state()
        agents$vax_due_day[i] <- dose_sim_day
      }
    } else {
      agents$intended_vax[i] <- FALSE
      agents$state[i]        <- maternal_state()
      agents$vax_due_day[i]  <- NA_integer_
    }
  }
  ri <- params$role_init
  for (rl in c("toddler","adult","elderly","staff")) {
    idx <- which(agents$role == rl)
    if (length(idx) == 0) next
    spec <- ri[[rl]]
    if (is.null(spec)) spec <- list(n_init = 0L, p_immune = 0)
    agents$n_infections[idx] <- spec$n_init
    agents$vax_due_day[idx]  <- NA_integer_
    if (!is.null(spec$p_immune) && spec$p_immune > 0) {
      immune <- runif(length(idx)) < spec$p_immune
      agents$n_infections[idx[immune]] <- spec$n_init + 1L
    }
    agents$state[idx] <- "S"
  }
  agents$initial_role <- agents$role
  agents
}

# ============================================================================
# SECTION 5: DIRECTIONAL TRANSMISSION CASCADE  [MODEL A: M_b full, S_ab asymptomatic]
# ============================================================================
transmission_step <- function(agents, day, params) {
  n   <- nrow(agents)
  s0  <- agents$state
  rho_mab <- if (!is.null(params$rho_mab)) params$rho_mab else 1   
  
  ## [MODEL A] susceptibility by state: S and S_ab are fully susceptible (1);
  ## M_b (full protection), MI_ab (immune), M, I, I_ab are not (0).
  ## (theta_mab no longer used — protection is the asymptomatic shunt.)
  sus_of <- function(st) if (st == "S" || st == "S_ab") 1 else 0
  
  beta_dc <- params$beta_daycare   * params$daycare_hours
  beta_hh <- params$beta_household * params$household_hours
  lam_com_base <- lambda_ext_seasonal(day, params)
  w_in <- params$cluster_within_weight
  w_bw_child <- if (!is.null(params$cluster_between_weight_child)) params$cluster_between_weight_child else params$cluster_between_weight
  w_bw_staff <- if (!is.null(params$cluster_between_weight_staff)) params$cluster_between_weight_staff else params$cluster_between_weight
  rm_roles <- rownames(params$C_daycare)
  
  # ---- start-of-day effective infectiousness (prior-day infectious agents) ----

  eff <- numeric(n)
  for (j in which(s0 == "I" | s0 == "I_ab")) {
    ct <- compute_Ct_direct(agents$days_in_state[j], agents$role[j], params)
    eff[j] <- get_kappa(agents$role[j], params) *
      get_rho(agents$n_infections[j], params) *
      (params$beta1 + params$beta2 * f_Ct(ct, params)) *
      (if (s0[j] == "I_ab") rho_mab else 1)         
  }
  
  # ---- mutable working copies ----
  state   <- s0
  source  <- agents$infection_source
  ninf    <- agents$n_infections
  dis     <- agents$days_in_state
  ever    <- agents$ever_infected
  outcome <- agents$infected_outcome
  broke   <- agents$breakthrough
  dayinf  <- agents$day_infected
  recov   <- agents$recovery_day
  newly   <- integer(0)
  newly_by <- integer(0)
  dc_today <- integer(0)
  
  commit <- function(i, src, inf_by = NA_integer_) {
    ## [MODEL A] infection during the mAb window (state S_ab) is ASYMPTOMATIC:
    ## -> state I_ab, flagged breakthrough, NOT counted as infected_outcome,
    ##    but still infectious and confers immunity (recovers to MI_ab).
    ## Infection from an ordinary susceptible (S) is SYMPTOMATIC -> state I.
    asymp <- (state[i] == "S_ab")
    if (asymp) broke[i] <<- TRUE
    state[i]   <<- if (asymp) "I_ab" else "I"
    ninf[i]    <<- ninf[i] + 1L
    dis[i]     <<- 0L
    ever[i]    <<- TRUE
    if (!asymp && day >= params$outcome_from_day) outcome[i] <<- TRUE  
    dayinf[i]  <<- day
    source[i]  <<- src
    recov[i]   <<- day + get_inf_duration(ninf[i], params)
    newly      <<- c(newly, i)
    newly_by   <<- c(newly_by, if (is.na(inf_by)) NA_integer_ else inf_by)
  }
  
  hh_members <- split(seq_len(n), agents$hh_id)
  household_pass <- function(inf_eff, src_label) {
    for (mem in hh_members) {
      if (length(mem) < 2) next
      if (all(inf_eff[mem] == 0)) next
      roles_h <- ifelse(agents$role[mem] == "staff", "adult", agents$role[mem])
      rc <- table(roles_h)
      for (ii in seq_along(mem)) {
        i <- mem[ii]
        sus_i <- sus_of(state[i]); if (sus_i == 0) next
        role_i_h <- roles_h[ii]
        contrib <- numeric(length(mem))
        for (jj in seq_along(mem)) {
          j <- mem[jj]
          if (jj == ii || inf_eff[j] == 0) next
          cc <- params$C_household[role_i_h, roles_h[jj]]
          if (cc > 0) contrib[jj] <- cc * inf_eff[j] / as.numeric(rc[roles_h[jj]])
        }
        lam <- sum(contrib)
        if (lam > 0) {
          phi_i   <- get_phi_role(agents$role[i], params)
          sigma_i <- get_sigma(ninf[i], params)
          if (runif(1) < (1 - exp(-phi_i * sigma_i * sus_i * beta_hh * lam))) {
            pos <- which(contrib > 0)
            commit(i, src_label, pick_infector(mem[pos], contrib[pos]))
          }
        }
      }
    }
  }
  
  ## ---- PHASE 0: COMMUNITY IMPORTATION ----
  for (i in which(state == "S" | state == "S_ab")) {       
    phi_i   <- get_phi_role(agents$role[i], params)
    sigma_i <- get_sigma(ninf[i], params)
    sus_i   <- sus_of(state[i])
    if (sigma_i * sus_i > 0 &&
        runif(1) < (1 - exp(-phi_i * sigma_i * sus_i * lam_com_base)))
      commit(i, "community")
  }
  
  ## ---- PHASE 1: HOME prior-day infectious transmit to household ----
  household_pass(eff, "household")
  
  ## ---- PHASE 2: DAYCARE — infectious children/staff transmit at the center ----
  present <- logical(n)
  for (i in which(agents$is_daycare)) {
    if (agents$role[i] == "infant" &&
        agents$age_days[i] < params$daycare_entry_age_days) next
    well_p <- if (agents$role[i] %in% c("infant","toddler"))
      params$attendance_child else params$attendance_staff
    ## [MODEL A] only SYMPTOMATIC cases stay home at sick_attendance;
    ## asymptomatic I_ab attend normally (well_p) and so transmit at the centre.
    p <- if (s0[i] == "I") params$sick_attendance else well_p
    present[i] <- runif(1) < p
  }
  present_idx     <- which(present)
  present_by_role <- split(present_idx, agents$role[present_idx])
  
  cl <- agents$cluster_id
  for (i in present_idx) {
    sus_i <- sus_of(state[i]); if (sus_i == 0) next
    role_i <- agents$role[i]
    if (!(role_i %in% rm_roles)) next
    cl_i <- cl[i]; lam_dc <- 0
    contrib_idx <- integer(0); contrib_val <- numeric(0)
    for (J in rm_roles) {
      idxJ <- present_by_role[[J]]
      if (is.null(idxJ)) next
      idxJ <- idxJ[idxJ != i]
      if (length(idxJ) == 0) next
      w_bw_eff <- if (role_i == "staff" || J == "staff") w_bw_staff else w_bw_child
      same  <- !is.na(cl_i) & !is.na(cl[idxJ]) & (cl[idxJ] == cl_i)
      S_J   <- sum(same) * w_in + sum(!same) * w_bw_eff
      if (S_J <= 0) next
      effJ  <- eff[idxJ]
      cJ    <- params$C_daycare[role_i, J]
      numer <- w_in * sum(effJ[same]) + w_bw_eff * sum(effJ[!same])
      if (numer > 0)
        lam_dc <- lam_dc + cJ * numer / S_J
      wJ   <- ifelse(same, w_in, w_bw_eff)
      pj   <- cJ * wJ * effJ / S_J
      keep <- effJ > 0 & cJ > 0 & wJ > 0
      if (any(keep)) { contrib_idx <- c(contrib_idx, idxJ[keep]); contrib_val <- c(contrib_val, pj[keep]) }
    }
    if (lam_dc > 0) {
      phi_i   <- get_phi_role(role_i, params)
      sigma_i <- get_sigma(ninf[i], params)
      if (runif(1) < (1 - exp(-phi_i * sigma_i * sus_i * beta_dc * lam_dc))) {
        commit(i, paste0("daycare:C", cl_i), pick_infector(contrib_idx, contrib_val))
        if (role_i %in% c("infant","toddler")) dc_today <- c(dc_today, i)
      }
    }
  }
  
  ## ---- PHASE 3: HOME  — today's daycare infectees carry it home ----
  if (isTRUE(params$daycare_bring_home_sameday) && length(dc_today) > 0) {
    eff_bh <- numeric(n)
    for (j in dc_today) {
      ct <- compute_Ct_direct(0L, agents$role[j], params)
      eff_bh[j] <- get_kappa(agents$role[j], params) *
        get_rho(ninf[j], params) *
        (params$beta1 + params$beta2 * f_Ct(ct, params)) *
        (if (state[j] == "I_ab") rho_mab else 1)     # [MODEL A] asymptomatic breakthroughs shed less
    }
    household_pass(eff_bh, "household:from_daycare")
  }
  
  agents$state            <- state
  agents$infection_source <- source
  agents$n_infections     <- ninf
  agents$days_in_state    <- dis
  agents$ever_infected    <- ever
  agents$infected_outcome <- outcome
  agents$breakthrough     <- broke
  agents$day_infected     <- dayinf
  agents$recovery_day     <- recov
  list(agents = agents, new_infections = newly, new_infectors = newly_by)
}

# ============================================================================
# SECTION 6b: DAILY TRANSITIONS & AGE ADVANCEMENT  [MODEL A staging]
# ============================================================================
daily_transitions <- function(agents, params, day) {
  for (i in seq_len(nrow(agents))) {
    # dosing: an eligible S/M infant receives mAb -> enters full-protection M_b
    if (!is.na(agents$vax_due_day[i]) && day >= agents$vax_due_day[i]) {
      if (agents$state[i] %in% c("S","M")) {
        agents$state[i]         <- "M_b"
        agents$days_in_state[i] <- 0L
        agents$vax_due_day[i]   <- NA_integer_
      }
    }
    s <- agents$state[i]
    if (s == "M") {
      if (runif(1) < (1 - exp(-params$omega_M))) { agents$state[i] <- "S"; agents$days_in_state[i] <- 0L }
      
    } else if (s == "M_b") {                                  # [MODEL A] full protection -> S_ab
      if (agents$days_in_state[i] > 0 && runif(1) < (1 - exp(-params$omega_mab_fix))) {
        agents$state[i] <- "S_ab"; agents$days_in_state[i] <- 0L
      }
      
    } else if (s == "S_ab") {                               
      if (runif(1) < (1 - exp(-params$omega_Mb))) {
        agents$state[i] <- "S"; agents$days_in_state[i] <- 0L
      }
      
    } else if (s == "MI_ab") {                              
      if (runif(1) < (1 - exp(-params$omega_Mb))) {
        agents$state[i] <- "S"; agents$days_in_state[i] <- 0L
      }
      
    } else if (s == "I_ab") {                               
      if (!is.na(agents$recovery_day[i]) && day >= agents$recovery_day[i]) {
        agents$state[i] <- "MI_ab"; agents$days_in_state[i] <- 0L; agents$recovery_day[i] <- NA_integer_
      }
      
    } else if (s == "I") {                                  
      if (!is.na(agents$recovery_day[i]) && day >= agents$recovery_day[i]) {
        agents$state[i] <- "S"; agents$days_in_state[i] <- 0L; agents$recovery_day[i] <- NA_integer_
      }
    }
  }
  agents$days_in_state <- agents$days_in_state + 1L
  age_idx <- which(agents$role %in% c("infant","toddler"))
  if (length(age_idx) > 0) agents$age_days[age_idx] <- agents$age_days[age_idx] + 1L
  agents
}

# ============================================================================
# SECTION 7: COUPLED SIMULATION
# ============================================================================
run_simulation <- function(params, vax_coverage = 0.0, record_panel = FALSE, quiet = TRUE) {
  agents <- initialize_agents(params, vax_coverage, quiet = quiet)
  n_days <- params$n_days
  n_ag   <- nrow(agents)
  
  if (record_panel) {
    inf_mat <- matrix(0L,       n_ag, n_days)
    vl_mat  <- matrix(0,        n_ag, n_days)
    ct_mat  <- matrix(NA_real_, n_ag, n_days)
  }
  
  daily <- data.frame(
    day = 1:n_days, dc_S = 0L, dc_I = 0L, dc_M = 0L, dc_Mb = 0L,
    dc_Sab = 0L, dc_Iab = 0L,
    dc_I_sym = 0L, hh_I_sym = 0L,
    dc_I1 = 0L, dc_I2 = 0L, dc_I3 = 0L, dc_I4plus = 0L, dc_new_inf = 0L,
    hh_S = 0L, hh_I = 0L, hh_new_inf = 0L, total_I = 0L, total_new = 0L,
    mean_Ct_infectious = numeric(n_days),
    inf_dc_I = 0L, tod_dc_I = 0L, staff_I = 0L,
    inf_dc_Isym = 0L, tod_dc_Isym = 0L, hh_Isym = 0L
  )
  
  inf_log_list <- vector("list", n_days * 100); inf_log_idx <- 0L
  log_inf <- function(idx, d, source, infector = NA_integer_) {
    inf_log_idx <<- inf_log_idx + 1L
    role_ct <- if (agents$role[idx] == "staff") "adult" else agents$role[idx]
    inf_log_list[[inf_log_idx]] <<- data.frame(
      global_id = agents$global_id[idx], episode = agents$n_infections[idx],
      day_infected = d, source = source, role_at_infection = agents$role[idx],
      cluster_id = agents$cluster_id[idx], ct_role = role_ct,
      breakthrough = agents$breakthrough[idx],          # [MODEL A] TRUE = asymptomatic in-window
      infector = infector,
      stringsAsFactors = FALSE)
  }
  
  for (d in seq_len(n_days)) {
    step   <- transmission_step(agents, d, params)
    agents <- step$agents
    for (k in seq_along(step$new_infections))
      log_inf(step$new_infections[k], d,
              agents$infection_source[step$new_infections[k]], step$new_infectors[k])
    
    ## [MODEL A] infectious prevalence = symptomatic I + asymptomatic I_ab (both transmit)
    is_I <- agents$state %in% c("I","I_ab")
    is_new <- is_I & agents$days_in_state == 0L
    is_sym <- agents$state == "I"
    daily$dc_S[d] <- sum(agents$state[agents$is_daycare] == "S")
    daily$dc_M[d] <- sum(agents$state[agents$is_daycare] == "M")
    daily$dc_Mb[d] <- sum(agents$state[agents$is_daycare] == "M_b")
    daily$dc_Sab[d] <- sum(agents$state[agents$is_daycare] == "S_ab")
    daily$dc_Iab[d] <- sum(agents$state == "I_ab" & agents$is_daycare)
    daily$dc_I1[d] <- sum(is_I & agents$is_daycare & agents$n_infections == 1)
    daily$dc_I2[d] <- sum(is_I & agents$is_daycare & agents$n_infections == 2)
    daily$dc_I3[d] <- sum(is_I & agents$is_daycare & agents$n_infections == 3)
    daily$dc_I4plus[d] <- sum(is_I & agents$is_daycare & agents$n_infections >= 4)
    daily$hh_S[d] <- sum(agents$state[!agents$is_daycare] == "S")
    daily$total_I[d] <- sum(is_I)
    daily$total_new[d] <- sum(is_new)
    daily$inf_dc_I[d] <- sum(is_I & agents$is_daycare & agents$role == "infant")
    daily$tod_dc_I[d] <- sum(is_I & agents$is_daycare & agents$role == "toddler")
    daily$staff_I[d]  <- sum(is_I & agents$role == "staff")
    
    daily$inf_dc_Isym[d] <- sum(is_sym & agents$is_daycare & agents$role == "infant")
    daily$tod_dc_Isym[d] <- sum(is_sym & agents$is_daycare & agents$role == "toddler")
    daily$hh_Isym[d]      <- sum(is_sym & !agents$is_daycare)
    
    daily$dc_I[d] <- sum(is_I & agents$is_daycare)
    daily$hh_I[d] <- sum(is_I & !agents$is_daycare)
    daily$dc_I_sym[d] <- sum(is_sym & agents$is_daycare)  
    daily$hh_I_sym[d] <- sum(is_sym & !agents$is_daycare)
    daily$dc_new_inf[d] <- sum(is_new & agents$is_daycare)
    daily$hh_new_inf[d] <- sum(is_new & !agents$is_daycare)
    
    if (sum(is_I) > 0) {
      idx_I <- which(is_I)
      cts <- sapply(idx_I, function(idx)
        compute_Ct_direct(agents$days_in_state[idx], agents$role[idx], params))
      daily$mean_Ct_infectious[d] <- mean(cts)
      if (record_panel) {
        inf_mat[idx_I, d] <- 1L
        ct_mat[idx_I, d]  <- cts
        vl_mat[idx_I, d]  <- f_Ct(cts, params)
      }
    } else daily$mean_Ct_infectious[d] <- params$Ct_LOD
    
    agents <- daily_transitions(agents, params, d)
  }
  
  infection_log <- if (inf_log_idx > 0) do.call(rbind, inf_log_list[1:inf_log_idx]) else
    data.frame(global_id=integer(0), episode=integer(0), day_infected=integer(0),
               source=character(0), role_at_infection=character(0),
               cluster_id=integer(0), ct_role=character(0), breakthrough=logical(0),
               infector=integer(0), stringsAsFactors=FALSE)
  
  out <- list(daily = daily, agents_final = agents, infection_log = infection_log)
  if (record_panel) out$panel <- list(inf = inf_mat, vl = vl_mat, ct = ct_mat,
                                      agents = agents)
  out
}

# ============================================================================
# SECTION 8: SCENARIO ANALYSIS  -- unchanged (uses infected_outcome = symptomatic only)
# ============================================================================
run_scenarios <- function(params, vax_levels) {
  n_vc <- length(vax_levels); n_sim <- params$n_sims; base_seed <- 2026178815
  cols <- c("ar_total","ar_vax_infant","ar_unvax_infant","ar_toddler",
            "ar_adult","ar_elderly","ar_staff","ar_daycare_all","ar_household_all",
            "n_infected_at_dc","n_infected_at_hh","n_infected_at_community","hh_secondary",
            "ar_tod_dc","ar_tod_hh","ar_tod_com","ar_tod_dchh")
  results <- data.frame(vax_coverage = rep(vax_levels, each = n_sim),
                        sim = rep(1:n_sim, times = n_vc))
  for (cc in cols) results[[cc]] <- numeric(n_vc * n_sim)
  idx <- 0L
  for (vc in vax_levels) {
    cat(sprintf("  Coverage = %3.0f%%: ", vc * 100))
    for (sim in seq_len(n_sim)) {
      idx <- idx + 1L
      set.seed(base_seed + idx)
      is_very_first_run <- (vc == vax_levels[1] && sim == 1L)
      res <- run_simulation(params, vax_coverage = vc,quiet = !is_very_first_run)
      af  <- res$agents_final
      staff_hh_ids <- unique(af$hh_id[af$role == "staff"])
      core <- !(af$hh_id %in% staff_hh_ids)
      m <- function(mask) if (sum(mask) > 0) mean(af$infected_outcome[mask]) else NA_real_
      results$ar_total[idx]         <- mean(af$infected_outcome)
      results$ar_vax_infant[idx]    <- m(af$intended_vax  & af$role=="infant" & af$is_daycare)
      results$ar_unvax_infant[idx]  <- m(!af$intended_vax & af$role=="infant" & af$is_daycare)
      results$ar_toddler[idx]       <- m(af$role=="toddler")
      results$ar_adult[idx]         <- m(af$role=="adult"   & core)
      results$ar_elderly[idx]       <- m(af$role=="elderly" & core)
      results$ar_staff[idx]         <- m(af$role=="staff")
      results$ar_daycare_all[idx]   <- m(af$is_daycare)
      results$ar_household_all[idx] <- m(!af$is_daycare & core)
      src <- res$infection_log$source
      results$n_infected_at_dc[idx]        <- sum(grepl("^daycare", src))
      results$n_infected_at_hh[idx]        <- sum(src == "household")
      results$n_infected_at_community[idx] <- sum(src == "community")
      inf_inf_hhs <- unique(af$hh_id[af$role=="infant" & af$infected_outcome])
      results$hh_secondary[idx] <- m(!af$is_daycare & core & af$hh_id %in% inf_inf_hhs)
      tod_ids <- af$global_id[af$role == "toddler"]
      n_tod   <- length(tod_ids)
      log_tod <- res$infection_log[
        res$infection_log$role_at_infection == "toddler" &
          res$infection_log$day_infected      >= params$outcome_from_day, ]
      tod_dc  <- unique(log_tod$global_id[grepl("^daycare",   log_tod$source)])
      tod_hh  <- unique(log_tod$global_id[grepl("^household", log_tod$source)])
      tod_com <- unique(log_tod$global_id[log_tod$source == "community"])
      results$ar_tod_dc[idx]   <- if (n_tod > 0) sum(tod_ids %in% tod_dc)             / n_tod else NA_real_
      results$ar_tod_hh[idx]   <- if (n_tod > 0) sum(tod_ids %in% tod_hh)             / n_tod else NA_real_
      results$ar_tod_com[idx]  <- if (n_tod > 0) sum(tod_ids %in% tod_com)            / n_tod else NA_real_
      results$ar_tod_dchh[idx] <- if (n_tod > 0) sum(tod_ids %in% union(tod_dc,tod_hh)) / n_tod else NA_real_
    }
    cat("done\n")
  }
  results
}

# ============================================================================
# SECTION 9: LONGITUDINAL EXPORT  -- unchanged
# ============================================================================
generate_hhbayes_longitudinal <- function(agents_final, infection_log, params,
                                          sample_interval = 4) {
  library(dplyr); library(tidyr)
  sample_days  <- seq(1, params$n_days, by = sample_interval)
  staff_hh_ids <- unique(agents_final$hh_id[agents_final$role == "staff"])
  target_agents <- agents_final %>%
    filter(!(hh_id %in% staff_hh_ids) | is_daycare == TRUE) %>%
    select(hh_id, global_id, role, is_daycare, age_days) %>%
    arrange(hh_id, global_id) %>%
    group_by(hh_id) %>%
    mutate(local_id = row_number(), person_id_str = paste0("HH_", hh_id, "_INV_", local_id)) %>%
    ungroup()
  compute_Ct_vec <- function(tau_vec, role_vec, params) {
    ct_res <- rep(params$Ct_LOD, length(tau_vec))
    for (r in unique(role_vec)) {
      idx <- which(role_vec == r & tau_vec >= 0); if (length(idx) == 0) next
      cp <- params$params_ct_values[[r]]$Cpeak; tp <- params$params_ct_values[[r]]$t_peak
      rv <- params$params_ct_values[[r]]$r;     dv <- params$params_ct_values[[r]]$d
      t_val <- tau_vec[idx]
      ct_res[idx[t_val <= tp]] <- cp + rv * (tp - t_val[t_val <= tp])
      ct_res[idx[t_val >  tp]] <- cp + dv * (t_val[t_val > tp] - tp)
    }
    ct_res[ct_res > params$Ct_LOD] <- params$Ct_LOD; ct_res
  }
  df_base <- expand_grid(global_id = target_agents$global_id, day_index = sample_days) %>%
    mutate(pcr_samples = params$Ct_LOD, active_episode = 0L)
  if (nrow(infection_log) > 0) {
    log_sub <- infection_log %>% filter(global_id %in% target_agents$global_id)
    if (nrow(log_sub) > 0) {
      episode_cts <- log_sub %>% tidyr::crossing(day_index = sample_days) %>%
        mutate(tau = day_index - day_infected) %>% filter(tau >= 0, tau <= 30) %>%
        mutate(ct_val = compute_Ct_vec(tau, ct_role, params))
      episode_min <- episode_cts %>% group_by(global_id, day_index) %>%
        slice_min(ct_val, n = 1, with_ties = FALSE) %>% ungroup() %>%
        select(global_id, day_index, ct_val, episode)
      df_base <- df_base %>% left_join(episode_min, by = c("global_id","day_index")) %>%
        mutate(pcr_samples = ifelse(!is.na(ct_val), ct_val, pcr_samples),
               active_episode = ifelse(!is.na(episode), as.integer(episode), 0L)) %>%
        select(-ct_val, -episode)
    }
  }
  df_computed <- df_base %>%
    mutate(pcr_samples = round(pcr_samples, 2),
           test_result = if_else(pcr_samples < params$detect_threshold_Ct, 1L, 0L),
           episode_id  = if_else(test_result == 1L, active_episode, 0L)) %>%
    left_join(target_agents %>% select(global_id, hh_id, role, is_daycare, age_days, person_id_str),
              by = "global_id") %>% select(-active_episode) %>%
    mutate(age_at_sample = ifelse(role %in% c("infant","toddler"),
                                  age_days - (params$n_days - day_index), NA_integer_))
  df_household <- df_computed %>% filter(!(hh_id %in% staff_hh_ids)) %>%
    mutate(hh_id_str = paste0("HH_", hh_id)) %>%
    select(hh_id = hh_id_str, person_id = person_id_str, role, day_index,
           pcr_samples, test_result, episode_id) %>% arrange(hh_id, person_id, day_index)
  df_daycare <- df_computed %>% filter(is_daycare == TRUE) %>%
    mutate(dc_id = "DC_1", dc_group = role) %>%
    select(dc_id, person_id = person_id_str, role, dc_group, day_index,
           pcr_samples, test_result, episode_id) %>% arrange(dc_id, person_id, day_index)
  list(household = df_household, daycare = df_daycare)
}

# ============================================================================
# SECTION 9b: PER-INDIVIDUAL x PER-DAY PANELS  -- unchanged
# ============================================================================
make_longitudinal <- function(params, vax_coverage = 0.0,
                              sample_interval = 1L, seed = 2026178815) {
  set.seed(seed)
  res  <- run_simulation(params, vax_coverage = vax_coverage, record_panel = TRUE)
  p    <- res$panel; ag <- p$agents
  n_ag <- nrow(ag); n_days <- params$n_days
  days <- seq(1L, n_days, by = sample_interval)
  reps <- length(days)
  local_id  <- as.integer(ave(ag$global_id, ag$hh_id,
                              FUN = function(x) rank(x, ties.method = "first")))
  person_id <- paste0("HH_", ag$hh_id, "_ID_", local_id)
  idx_rep <- rep(seq_len(n_ag), times = reps)
  panel <- data.frame(
    individual_id      = person_id[idx_rep],
    hh_id              = ag$hh_id[idx_rep],
    role               = ag$role[idx_rep],
    is_daycare         = ag$is_daycare[idx_rep],
    day                = rep(days, each = n_ag),
    infection_status   = as.vector(p$inf[, days, drop = FALSE]),
    viral_load         = round(as.vector(p$vl[, days, drop = FALSE]), 4),
    ct_value           = round(as.vector(p$ct[, days, drop = FALSE]), 2),
    vaccination_status = as.integer(ag$intended_vax)[idx_rep],
    stringsAsFactors   = FALSE
  )
  staff_hh <- unique(ag$hh_id[ag$role == "staff"])
  daycare_df <- panel[panel$is_daycare,
                      c("individual_id","role","day","infection_status","viral_load","ct_value","vaccination_status")]
  daycare_df <- daycare_df[order(daycare_df$individual_id, daycare_df$day), ]
  hh_df <- panel[!(panel$hh_id %in% staff_hh),
                 c("hh_id","individual_id","role","day","infection_status","viral_load","ct_value","vaccination_status")]
  hh_df <- hh_df[order(hh_df$hh_id, hh_df$individual_id, hh_df$day), ]
  list(daycare = daycare_df, household = hh_df)
}

# ============================================================================
# SECTION 10-12: EXECUTION, FIGURES, EXPORT  -- unchanged (inherits Model A behaviour)
# ============================================================================
if (RUN_ANALYSIS) {
  vax_levels       <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
  scenario_results <- run_scenarios(params, vax_levels)
  
  set.seed(2026178815)
  snap <- initialize_agents(params, 0)
  staff_hh <- unique(snap$hh_id[snap$role == "staff"]); core <- !(snap$hh_id %in% staff_hh)
  core_hh  <- unique(snap$hh_id[core])
  hh_has_inf <- tapply(snap$role=="infant",  snap$hh_id, any)[as.character(core_hh)]
  hh_has_tod <- tapply(snap$role=="toddler", snap$hh_id, any)[as.character(core_hh)]
  n_infant_dc  <- sum(snap$role=="infant"  & snap$is_daycare)
  n_toddler_dc <- sum(snap$role=="toddler" & snap$is_daycare)
  n_staff_dc   <- sum(snap$role=="staff")
  n_adults     <- sum(snap$role=="adult"   & core)
  n_elderly    <- sum(snap$role=="elderly" & core)
  
  cat("\n================ DEMOGRAPHICS SNAPSHOT ================\n")
  cat(sprintf(" Total agents: %d | core HH: %d | staff HH: %d\n",
              nrow(snap), length(core_hh), length(staff_hh)))
  cat(sprintf(" Infants(DC): %d  Toddlers(DC): %d  Staff: %d  Adults: %d  Elderly: %d\n",
              n_infant_dc, n_toddler_dc, n_staff_dc, n_adults, n_elderly))
  cat(sprintf(" HH types -> infant+toddler: %d | infant-only: %d | toddler-only: %d | no-kid: %d\n",
              sum(hh_has_inf & hh_has_tod), sum(hh_has_inf & !hh_has_tod),
              sum(!hh_has_inf & hh_has_tod), sum(!hh_has_inf & !hh_has_tod)))
  cat("Cluster composition (cluster_id x role) — now MIXED age:\n")
  # print(table(snap$cluster_id[snap$is_daycare], snap$role[snap$is_daycare],
  #             dnn = c("cluster","role")))
  
  # ── Cluster composition table ──────────────────────────────────────────────
  dc_snap <- snap[snap$is_daycare, c("cluster_id", "role")]
  
  # Count each role per cluster
  clust_inf   <- tapply(dc_snap$role == "infant",  dc_snap$cluster_id, sum)
  clust_tod   <- tapply(dc_snap$role == "toddler", dc_snap$cluster_id, sum)
  clust_staff <- tapply(dc_snap$role == "staff",   dc_snap$cluster_id, sum)
  
  all_cids <- sort(unique(dc_snap$cluster_id))
  
  clust_tab <- data.frame(
    cluster    = all_cids,
    group = ifelse(all_cids <= params$n_dc_clusters_g1,
                   "younger (≤18m)", "older (>18m)"),
    infants    = as.integer(clust_inf [as.character(all_cids)]),
    toddlers   = as.integer(clust_tod [as.character(all_cids)]),
    staff      = as.integer(clust_staff[as.character(all_cids)]),
    total_children = as.integer(clust_inf[as.character(all_cids)]) +
      as.integer(clust_tod[as.character(all_cids)])
  )
  clust_tab[is.na(clust_tab)] <- 0L
  
  cat("\n── Cluster composition ──────────────────────────────────────────\n")
  cat(sprintf("%-10s %-16s %-9s %-10s %-7s %-15s\n",
              "Cluster", "Group", "Infants", "Toddlers", "Staff", "Total children"))
  cat(strrep("-", 67), "\n")
  for (r in seq_len(nrow(clust_tab))) {
    cat(sprintf("%-10d %-16s %-9d %-10d %-7d %-15d\n",
                clust_tab$cluster[r], clust_tab$group[r],
                clust_tab$infants[r], clust_tab$toddlers[r],
                clust_tab$staff[r],   clust_tab$total_children[r]))
  }
  cat(strrep("-", 67), "\n")
  cat(sprintf("%-10s %-16s %-9d %-10d %-7d %-15d\n",
              "TOTAL", "",
              sum(clust_tab$infants), sum(clust_tab$toddlers),
              sum(clust_tab$staff),   sum(clust_tab$total_children)))
  cat("─────────────────────────────────────────────────────────────────\n\n")
  
  agg <- function(vc, col) mean(scenario_results[[col]][scenario_results$vax_coverage==vc], na.rm=TRUE)
  ar_cols <- c("ar_total","ar_vax_infant","ar_unvax_infant","ar_toddler",
               "ar_adult","ar_elderly","ar_staff")
  summary_df <- do.call(rbind, lapply(vax_levels, function(vc) {
    row <- data.frame(vax_pct = sprintf("%.0f%%", vc*100), vax_coverage = vc)
    for (cc in ar_cols) row[[cc]] <- agg(vc, cc)
    row$n_dc  <- agg(vc,"n_infected_at_dc")
    row$n_hh  <- agg(vc,"n_infected_at_hh")
    row$n_com <- agg(vc,"n_infected_at_community")
    row
  }))
  
  col_at  <- function(col, vc) scenario_results[[col]][scenario_results$vax_coverage == vc]
  prot_ci <- function(num, den, paired = FALSE, B = 2000) {
    if (paired) { ok <- !is.na(num) & !is.na(den); num <- num[ok]; den <- den[ok] }
    else        { num <- num[!is.na(num)]; den <- den[!is.na(den)] }
    if (length(num) == 0 || length(den) == 0 || mean(den) == 0)
      return(c(est = NA, lo = NA, hi = NA))
    bs <- numeric(B)
    for (b in seq_len(B)) {
      if (paired) { i <- sample.int(length(num), replace = TRUE); bn <- num[i]; bd <- den[i] }
      else        { bn <- sample(num, replace = TRUE);            bd <- sample(den, replace = TRUE) }
      bs[b] <- 1 - mean(bn) / mean(bd)
    }
    c(est = 1 - mean(num)/mean(den),
      lo  = unname(quantile(bs, 0.025, na.rm = TRUE)),
      hi  = unname(quantile(bs, 0.975, na.rm = TRUE)))
  }
  
  set.seed(99)
  prot_list <- lapply(vax_levels, function(vc) {
    d  <- prot_ci(col_at("ar_vax_infant", vc),   col_at("ar_unvax_infant", vc), paired = TRUE)
    t  <- prot_ci(col_at("ar_vax_infant", vc),   col_at("ar_unvax_infant", 0))
    it      <- prot_ci(col_at("ar_toddler",    vc), col_at("ar_toddler",    0))
    it_dc   <- prot_ci(col_at("ar_tod_dc",     vc), col_at("ar_tod_dc",     0))
    it_hh   <- prot_ci(col_at("ar_tod_hh",     vc), col_at("ar_tod_hh",     0))
    it_dchh <- prot_ci(col_at("ar_tod_dchh",   vc), col_at("ar_tod_dchh",   0))
    ia      <- prot_ci(col_at("ar_adult",       vc), col_at("ar_adult",       0))
    ie      <- prot_ci(col_at("ar_elderly",     vc), col_at("ar_elderly",     0))
    ov      <- prot_ci(col_at("ar_total",       vc), col_at("ar_total",       0))
    data.frame(vax_pct = sprintf("%.0f%%", vc*100), vax_coverage = vc,
               direct_est=d[1],       direct_lo=d[2],       direct_hi=d[3],
               total_est=t[1],        total_lo=t[2],        total_hi=t[3],
               ind_tod_est=it[1],     ind_tod_lo=it[2],     ind_tod_hi=it[3],
               ind_tod_dc_est=it_dc[1],     ind_tod_dc_lo=it_dc[2],     ind_tod_dc_hi=it_dc[3],
               ind_tod_hh_est=it_hh[1],     ind_tod_hh_lo=it_hh[2],     ind_tod_hh_hi=it_hh[3],
               ind_tod_dchh_est=it_dchh[1], ind_tod_dchh_lo=it_dchh[2], ind_tod_dchh_hi=it_dchh[3],
               ind_adult_est=ia[1],   ind_adult_lo=ia[2],   ind_adult_hi=ia[3],
               ind_eldly_est=ie[1],   ind_eldly_lo=ie[2],   ind_eldly_hi=ie[3],
               overall_est=ov[1],     overall_lo=ov[2],     overall_hi=ov[3],
               row.names = NULL)
  })
  prot_df <- do.call(rbind, prot_list)
  
  cat("\n================ RESULTS: Attack Rates (mean) ================\n")
  cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-8s %-8s\n",
              "Cov","Tot","VaxInf","UVInf","Toddler","Adult","Eldly"))
  for (i in seq_len(nrow(summary_df))) {
    f <- function(x) ifelse(is.na(x)," N/A   ", sprintf("%-8.3f", x))
    cat(sprintf("%-6s %s%s%s%s%s%s\n", summary_df$vax_pct[i],
                f(summary_df$ar_total[i]), f(summary_df$ar_vax_infant[i]), f(summary_df$ar_unvax_infant[i]),
                f(summary_df$ar_toddler[i]), f(summary_df$ar_adult[i]), f(summary_df$ar_elderly[i])))
  }
  
  fmtp <- function(e,l,h) if (is.na(e)) "     N/A        " else
    sprintf("%5.1f [%5.1f,%5.1f]", e*100, l*100, h*100)
  cat("\n--- Protection %, est [95% CI] ---\n")
  for (i in seq_len(nrow(prot_df))) {
    cat(sprintf("\nCoverage %s:\n", prot_df$vax_pct[i]))
    cat(sprintf("  DIRECT   (infant: vax vs unvax)     : %s\n", fmtp(prot_df$direct_est[i], prot_df$direct_lo[i], prot_df$direct_hi[i])))
    cat(sprintf("  TOTAL    (infant: vax vs baseline)  : %s\n", fmtp(prot_df$total_est[i],  prot_df$total_lo[i],  prot_df$total_hi[i])))
    cat(sprintf("  INDIRECT (toddler, all-source)      : %s\n",
                fmtp(prot_df$ind_tod_est[i],      prot_df$ind_tod_lo[i],      prot_df$ind_tod_hi[i])))
    cat(sprintf("  INDIRECT (toddler, daycare-acquired): %s\n",
                fmtp(prot_df$ind_tod_dc_est[i],   prot_df$ind_tod_dc_lo[i],   prot_df$ind_tod_dc_hi[i])))
    cat(sprintf("  INDIRECT (toddler, hh-acquired)     : %s\n",
                fmtp(prot_df$ind_tod_hh_est[i],   prot_df$ind_tod_hh_lo[i],   prot_df$ind_tod_hh_hi[i])))
    cat(sprintf("  INDIRECT (toddler, dc+hh combined)  : %s\n",
                fmtp(prot_df$ind_tod_dchh_est[i], prot_df$ind_tod_dchh_lo[i], prot_df$ind_tod_dchh_hi[i])))
    cat(sprintf("  INDIRECT (adult)                    : %s\n",
                fmtp(prot_df$ind_adult_est[i],    prot_df$ind_adult_lo[i],    prot_df$ind_adult_hi[i])))
    cat(sprintf("  INDIRECT (older adult)              : %s\n", fmtp(prot_df$ind_eldly_est[i],prot_df$ind_eldly_lo[i],prot_df$ind_eldly_hi[i])))
    cat(sprintf("  OVERALL  (whole population)         : %s\n", fmtp(prot_df$overall_est[i],  prot_df$overall_lo[i],  prot_df$overall_hi[i])))
  }
  
  cat("\n--- Infection source (mean per sim) ---\n")
  cat(sprintf("%-6s %-12s %-12s %-12s\n","Cov","Daycare","Household","Community"))
  for (i in seq_len(nrow(summary_df)))
    cat(sprintf("%-6s %-12.1f %-12.1f %-12.1f\n", summary_df$vax_pct[i],
                summary_df$n_dc[i], summary_df$n_hh[i], summary_df$n_com[i]))
  
  # ============================================================================
  # SECTION 11: FIGURES
  # ============================================================================
  out_dir <- getwd()
  N_plot      <- params$n_sims
  fig_cov_hi  <- 1
  n_days      <- params$n_days
  base_seed_p <- 2026178815
  days        <- 1:n_days
  month_starts<- c(1,32,63,93,124,154,184,215,243,274,304,335)
  month_labels<- c("Jul","Aug","Sep","Oct","Nov","Dec",
                   "Jan","Feb","Mar","Apr","May","Jun")
  
  sim_curves <- function(vc) {
    dc <- hh <- mb <- matrix(0, N_plot, n_days)
    for (r in seq_len(N_plot)) {
      set.seed(base_seed_p + r)
      tmp <- run_simulation(params, vax_coverage = vc)
      dc[r,] <- tmp$daily$dc_I_sym   # symptomatic only (asymptomatic I_ab excluded)
      hh[r,] <- tmp$daily$hh_I_sym   # symptomatic only
      mb[r,] <- tmp$daily$dc_Mb + tmp$daily$dc_Sab + tmp$daily$dc_Iab
    }
    list(dc = dc, hh = hh, mb = mb)
  }
  
  cat(sprintf("\n[Figure 1] averaging %d sims at 0%% and %.0f%% coverage ...\n",
              N_plot, fig_cov_hi*100))
  c0 <- sim_curves(0); c1 <- sim_curves(fig_cov_hi)
  q_lo <- function(x) apply(x, 2, quantile, 0.025)
  q_hi <- function(x) apply(x, 2, quantile, 0.975)
  ymax_I  <- max(q_hi(c0$dc), q_hi(c0$hh),
                 q_hi(c1$dc), q_hi(c1$hh)) + 1
  ymax_Mb <- max(q_hi(c0$mb), q_hi(c1$mb)) + 1
  if (ymax_Mb < 10) ymax_Mb <- 10
  
  draw_panel <- function(cc, ttl) {
    plot(days, colMeans(cc$dc), type="l", col="#f0027f", lwd=2.5, xaxt="n",
         xlab="Month", ylab="Symptomatic cases",
         main=ttl, ylim=c(0, ymax_I))
    polygon(c(days,rev(days)), c(q_lo(cc$dc),rev(q_hi(cc$dc))),
            col=adjustcolor("#f0027f", 0.15), border=NA)
    lines(days, colMeans(cc$hh), col="#386cb0", lwd=2.5)
    polygon(c(days,rev(days)), c(q_lo(cc$hh),rev(q_hi(cc$hh))),
            col=adjustcolor("#386cb0", 0.15), border=NA)
    axis(1, at=month_starts, labels=month_labels, cex.axis=0.8)
    par(new=TRUE)
    plot(days, colMeans(cc$mb), type="l", col="#005a32", lwd=2.5, lty=2,
         axes=FALSE, xlab="", ylab="", ylim=c(0, ymax_Mb))
    polygon(c(days,rev(days)), c(q_lo(cc$mb),rev(q_hi(cc$mb))),
            col=adjustcolor("#005a32", 0.15), border=NA)
    axis(4, col="#005a32", col.axis="#005a32", las=1)
    mtext("Protected infants (M_b)", side=4, line=2.5,
          col="#005a32", cex=0.8, font=2)
    legend("topright",
           legend=c("Daycare symptomatic cases","Household symptomatic cases",
                    "Protected (M_b)"),
           col=c("#f0027f","#386cb0","#005a32"),
           lwd=2.5, lty=c(1,1,2), cex=0.65, bty="n")
  }
  
  f1 <- file.path(out_dir, "results1/RSV_fig1_epidemic_curves.pdf")
  pdf(f1, width = 11, height = 4.5)
  par(mfrow = c(1,2), mar = c(4.5,4.5,3,4))
  draw_panel(c0, "A) 0% coverage")
  draw_panel(c1, sprintf("B) %.0f%% coverage", fig_cov_hi*100))
  dev.off()
  cat(sprintf("[Figure 1] saved -> %s\n", f1))
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    pull <- function(df, est, lo, hi, label) data.frame(
      vax = df$vax_coverage * 100, metric = label,
      est = df[[est]] * 100, lo = df[[lo]] * 100, hi = df[[hi]] * 100)
    inf_df <- rbind(
      pull(prot_df, "direct_est", "direct_lo", "direct_hi", "Direct (infant)"),
      pull(prot_df, "total_est",  "total_lo",  "total_hi",  "Total (infant)"))
    inf_df <- inf_df[!is.na(inf_df$est), ]
    tod_df <- rbind(
      pull(prot_df, "ind_tod_est",    "ind_tod_lo",    "ind_tod_hi",    "All-source"),
      pull(prot_df, "ind_tod_dc_est", "ind_tod_dc_lo", "ind_tod_dc_hi", "Daycare-acquired"),
      pull(prot_df, "ind_tod_hh_est", "ind_tod_hh_lo", "ind_tod_hh_hi", "Household-acquired"),
      pull(prot_df, "ind_tod_dchh_est","ind_tod_dchh_lo","ind_tod_dchh_hi","DC + HH combined"))
    tod_df <- tod_df[!is.na(tod_df$est), ]
    ind_df <- rbind(
      pull(prot_df, "ind_adult_est", "ind_adult_lo", "ind_adult_hi", "Adult"),
      pull(prot_df, "ind_eldly_est", "ind_eldly_lo", "ind_eldly_hi", "Older adult"))
    ind_df <- ind_df[!is.na(ind_df$est), ]
    dodge    <- position_dodge(width = 5)
    base_thm <- theme_classic(base_size = 11) +
      theme(legend.position = "top", legend.title = element_blank(),
            panel.grid.major.y = element_line(colour = "grey92"),
            legend.text = element_text(size = 9))
    tod_cols <- c(
      "All-source"       = "grey60",
      "Daycare-acquired"                     = "#e6550d",
      "Household-acquired"                   = "#3182bd",
      "DC + HH combined"   = "#31a354")
    mk <- function(d, ttl, col_map = NULL) {
      p <- ggplot(d, aes(vax, est, colour = metric, group = metric)) +
        geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
        geom_line(position = dodge, linewidth = 0.9) +
        geom_errorbar(aes(ymin = lo, ymax = hi), width = 3.5, position = dodge) +
        geom_point(size = 2.8, position = dodge) +
        labs(title = ttl, x = "Infant mAb coverage (%)",
             y = "Protection (%)  [mean, 95% CI]") +
        base_thm
      if (!is.null(col_map)) p <- p + scale_colour_manual(values = col_map)
      p
    }
    p_inf <- mk(inf_df, "A) Infant protection (direct / total)")
    p_tod <- mk(tod_df, "B) Toddler indirect — by transmission route", tod_cols)
    p_ind <- mk(ind_df, "C) Indirect protection — adults & older adults")
    f2 <- file.path(out_dir, "results1/RSV_fig2_protection.pdf")
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      ggsave(f2, (p_inf / p_tod / p_ind), width = 9, height = 13)
    } else {
      ggsave(sub("\\.pdf$", "_infant.pdf",  f2), p_inf, width = 7, height = 4.5)
      ggsave(sub("\\.pdf$", "_toddler.pdf", f2), p_tod, width = 7, height = 5.0)
      ggsave(sub("\\.pdf$", "_indirect.pdf",f2), p_ind, width = 7, height = 4.5)
    }
    cat(sprintf("[Figure 2] saved -> %s\n", f2))
  } else {
    cat("[Figure 2] ggplot2 not installed; skipped (install.packages('ggplot2')).\n")
  }
  
  # ============================================================================
  # SECTION 12: EXCEL EXPORT
  # ============================================================================
  pct <- function(x) round(100 * x, 2)
  ar_tab <- data.frame(
    `Vaccination Coverage` = summary_df$vax_pct,
    `N Vax Infants`        = round(n_infant_dc * summary_df$vax_coverage),
    `N Unvax Infants`      = round(n_infant_dc * (1 - summary_df$vax_coverage)),
    `AR Vax Infants (%)`   = pct(summary_df$ar_vax_infant),
    `AR Unvax Infants (%)` = pct(summary_df$ar_unvax_infant),
    `Exp. Infected Vax`    = round(summary_df$ar_vax_infant   * n_infant_dc * summary_df$vax_coverage, 1),
    `Exp. Infected Unvax`  = round(summary_df$ar_unvax_infant * n_infant_dc * (1 - summary_df$vax_coverage), 1),
    `N Toddlers`           = n_toddler_dc,
    `AR Toddlers (%)`      = pct(summary_df$ar_toddler),
    `N Adults`             = n_adults,
    `AR Adults (%)`        = pct(summary_df$ar_adult),
    `N Elderly`            = n_elderly,
    `AR Elderly (%)`       = pct(summary_df$ar_elderly),
    `N Staff`              = n_staff_dc,
    `AR Staff (%)`         = pct(summary_df$ar_staff),
    `AR Total (%)`         = pct(summary_df$ar_total),
    check.names = FALSE)
  ci_str <- function(e,l,h) ifelse(is.na(e), "N/A",
                                   sprintf("%.1f [%.1f, %.1f]", e*100, l*100, h*100))
  prot_tab <- data.frame(
    `Vaccination Coverage`              = prot_df$vax_pct,
    `Direct - infant (%, 95% CI)`       = ci_str(prot_df$direct_est, prot_df$direct_lo, prot_df$direct_hi),
    `Total - infant (%, 95% CI)`        = ci_str(prot_df$total_est,  prot_df$total_lo,  prot_df$total_hi),
    `Indirect - toddler (%, 95% CI)`    = ci_str(prot_df$ind_tod_est,prot_df$ind_tod_lo,prot_df$ind_tod_hi),
    `Indirect - adult (%, 95% CI)`      = ci_str(prot_df$ind_adult_est,prot_df$ind_adult_lo,prot_df$ind_adult_hi),
    `Indirect - older adult (%, 95% CI)`= ci_str(prot_df$ind_eldly_est,prot_df$ind_eldly_lo,prot_df$ind_eldly_hi),
    `Overall (%, 95% CI)`               = ci_str(prot_df$overall_est,prot_df$overall_lo,prot_df$overall_hi),
    check.names = FALSE)
  xlsx_path <- file.path(out_dir, "results1/RSV_Simulation_Summary.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    library(openxlsx)
    wb <- createWorkbook()
    addWorksheet(wb, "Attack Rates"); writeDataTable(wb, 1, ar_tab,   tableStyle = "TableStyleLight9")
    addWorksheet(wb, "Protection");   writeDataTable(wb, 2, prot_tab, tableStyle = "TableStyleLight9")
    freezePane(wb, 1, firstActiveRow = 2); setColWidths(wb, 1, 1:ncol(ar_tab),   "auto")
    freezePane(wb, 2, firstActiveRow = 2); setColWidths(wb, 2, 1:ncol(prot_tab), "auto")
    saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    cat(sprintf("[Excel] saved -> %s\n", xlsx_path))
  } else {
    write.csv(ar_tab,   file.path(out_dir, "results1/RSV_attack_rates.csv"), row.names = FALSE)
    write.csv(prot_tab, file.path(out_dir, "results1/RSV_protection.csv"),   row.names = FALSE)
    cat("[Excel] openxlsx not installed; wrote CSVs instead.\n")
  }
  
  # ============================================================================
  # SECTION 13: PER-INDIVIDUAL LONGITUDINAL PANELS
  # ============================================================================
  longitudinal_coverage <- 0.60
  panel_interval        <- 1L
  long <- make_longitudinal(params, vax_coverage = longitudinal_coverage,
                            sample_interval = panel_interval)
  write.csv(long$daycare,   file.path(out_dir, "results1/RSV_daycare_panel.csv"),   row.names = FALSE)
  write.csv(long$household, file.path(out_dir, "results1/RSV_household_panel.csv"), row.names = FALSE)
}  # end if (RUN_ANALYSIS)

