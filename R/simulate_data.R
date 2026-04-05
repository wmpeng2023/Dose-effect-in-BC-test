# simulate_data.R
# Simulates breast cancer patient dataset with dose levels, Cu/Zn ratio,
# and dysgeusia toxicity events (CTCAE grading).
#
# Hypothesis encoded:
#   higher dose -> higher Cu/Zn ratio -> earlier onset + higher grade dysgeusia
# Reference Cu/Zn range for breast cancer patients: ~1.5-3.5
# (Kucharska et al., Maj et al., published oncology Cu/Zn literature)

simulate_bc_data <- function(n_patients = 100, seed = 42) {
  set.seed(seed)

  dose_levels <- c(1, 2, 4, 8, 15, 30)
  n_doses <- length(dose_levels)

  # Distribute patients evenly across dose levels
  n_per_dose <- ceiling(n_patients / n_doses)
  dose_assignments <- rep(dose_levels, n_per_dose)[1:n_patients]

  patient_data <- data.frame(
    patient_id = 1:n_patients,
    dose_mg    = sample(dose_assignments)  # randomize assignment order
  )

  dose_idx <- match(patient_data$dose_mg, dose_levels)  # 1-6

  # Cu/Zn ratio: increases with dose (treatment effect + noise)
  # At dose 1mg: mean ~1.5; at dose 30mg: mean ~2.6
  cu_zn_mean <- 1.50 + 0.22 * (dose_idx - 1)
  patient_data$cu_zn_ratio <- round(
    rnorm(n_patients, mean = cu_zn_mean, sd = 0.28),
    digits = 2
  )
  patient_data$cu_zn_ratio <- pmax(patient_data$cu_zn_ratio, 1.0)  # floor at 1.0

  # --- Dysgeusia event generation ---
  # Probability of reaching each grade (cumulative, must be non-increasing)
  # Grade 1 incidence: 20% (1mg) to 75% (30mg)
  # Grade 4 incidence: 2% (1mg) to 15% (30mg)
  p_g1 <- pmin(pmax(0.20 + 0.11 * (dose_idx - 1) + rnorm(n_patients, 0, 0.04), 0.05), 0.90)
  p_g2 <- pmin(pmax(0.12 + 0.09 * (dose_idx - 1) + rnorm(n_patients, 0, 0.04), 0.02), p_g1)
  p_g3 <- pmin(pmax(0.05 + 0.05 * (dose_idx - 1) + rnorm(n_patients, 0, 0.03), 0.01), p_g2)
  p_g4 <- pmin(pmax(0.01 + 0.025 * (dose_idx - 1) + rnorm(n_patients, 0, 0.02), 0.005), p_g3)

  p_grade <- cbind(p_g1, p_g2, p_g3, p_g4)

  event_rows <- vector("list", n_patients * 4)
  event_count <- 0

  for (i in seq_len(n_patients)) {
    pid    <- patient_data$patient_id[i]
    d_idx  <- dose_idx[i]
    cu_zn  <- patient_data$cu_zn_ratio[i]

    # Onset rate: higher dose and higher Cu/Zn -> faster onset
    onset_rate <- 0.04 + 0.015 * (d_idx - 1) + 0.008 * max(cu_zn - 1.5, 0)

    prev_time <- 0

    for (g in 1:4) {
      if (runif(1) > p_grade[i, g]) break  # Stops at first "no" — grade hierarchy enforced

      if (g == 1) {
        time_g <- max(1L, round(rexp(1, rate = onset_rate)))
      } else {
        inter_event <- max(2L, round(runif(1, 3, 14)))
        time_g <- prev_time + inter_event
      }

      event_count <- event_count + 1
      event_rows[[event_count]] <- data.frame(
        patient_id         = pid,
        dysgeusia_grade    = g,
        time_to_first_onset = time_g
      )
      prev_time <- time_g
    }
  }

  # Build event data frame
  if (event_count > 0) {
    event_data <- do.call(rbind, event_rows[seq_len(event_count)])
  } else {
    event_data <- data.frame(
      patient_id          = integer(0),
      dysgeusia_grade     = integer(0),
      time_to_first_onset = numeric(0)
    )
  }

  # Build full merged dataset
  # Event patients: one row per grade (grade >= 1)
  full_events <- merge(patient_data, event_data, by = "patient_id")

  # Event-free patients: one row with grade = 0, time = NA
  event_pids <- unique(event_data$patient_id)
  no_event_pids <- setdiff(patient_data$patient_id, event_pids)

  if (length(no_event_pids) > 0) {
    no_event_rows <- patient_data[patient_data$patient_id %in% no_event_pids, ]
    no_event_rows$dysgeusia_grade    <- 0L
    no_event_rows$time_to_first_onset <- NA_real_
    full_data <- rbind(full_events, no_event_rows)
  } else {
    full_data <- full_events
  }

  full_data <- full_data[order(full_data$patient_id, full_data$dysgeusia_grade), ]
  rownames(full_data) <- NULL

  # Label dose as factor for plotting
  full_data$dose_label <- factor(
    paste0(full_data$dose_mg, " mg"),
    levels = paste0(dose_levels, " mg")
  )

  list(
    patient_data = patient_data,
    event_data   = event_data,
    full_data    = full_data
  )
}
