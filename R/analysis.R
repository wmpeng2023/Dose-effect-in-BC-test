# analysis.R
# Statistical analysis functions:
#   1. Summary statistics
#   2. Dose-response (logistic regression, Spearman correlation)
#   3. Cu/Zn ratio vs dysgeusia grade (linear regression)
#   4. Time to dysgeusia (Cox proportional hazards, KM curves)

library(dplyr)
library(ggplot2)
library(plotly)
library(survival)
library(broom)

CENSORING_DAY <- 180  # Patients without events are censored at this day

# Helper: patient-level summary from full_data
patient_level <- function(full_data, dose_filter) {
  full_data %>%
    filter(dose_mg %in% as.numeric(dose_filter)) %>%
    group_by(patient_id, dose_mg, cu_zn_ratio) %>%
    summarise(
      max_grade    = max(dysgeusia_grade, na.rm = TRUE),
      any_dysgeusia = as.integer(max(dysgeusia_grade, na.rm = TRUE) > 0),
      time_to_any  = {
        t <- time_to_first_onset[dysgeusia_grade > 0]
        if (length(t) > 0 && any(!is.na(t))) min(t, na.rm = TRUE) else NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      surv_time  = ifelse(is.na(time_to_any), CENSORING_DAY, time_to_any),
      surv_event = any_dysgeusia
    )
}

# ---- 1. Summary Statistics ----
summarize_data <- function(data_list, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(data_list$full_data, dose_filter)

  cat("=== DATASET SUMMARY ===\n\n")
  cat("Total patients analyzed:", nrow(pl), "\n")
  cat("Dose levels included:   ", paste(sort(unique(pl$dose_mg)), collapse = ", "), "mg\n")
  cat("Patients with any dysgeusia:", sum(pl$any_dysgeusia), "\n")
  cat("Overall incidence:      ", round(mean(pl$any_dysgeusia) * 100, 1), "%\n\n")

  cat("--- Incidence & Grade by Dose Level ---\n")
  dose_tbl <- pl %>%
    group_by(`Dose (mg)` = dose_mg) %>%
    summarise(
      N                = n(),
      `N Dysgeusia`    = sum(any_dysgeusia),
      `Incidence (%)`  = round(mean(any_dysgeusia) * 100, 1),
      `Mean Cu/Zn`     = round(mean(cu_zn_ratio), 2),
      `SD Cu/Zn`       = round(sd(cu_zn_ratio), 2),
      `Mean Max Grade` = round(mean(max_grade), 2),
      .groups = "drop"
    )
  print(as.data.frame(dose_tbl), row.names = FALSE)

  cat("\n--- Cu/Zn Ratio Distribution (all patients) ---\n")
  cuzn_stats <- data.frame(
    Statistic = c("Min", "Q1", "Median", "Mean", "Q3", "Max"),
    Value     = round(quantile(pl$cu_zn_ratio, c(0, 0.25, 0.5, NA, 0.75, 1), na.rm = TRUE), 2)
  )
  cuzn_stats$Value[4] <- round(mean(pl$cu_zn_ratio), 2)
  print(cuzn_stats, row.names = FALSE)
}

# ---- 2. Dose-Response Analysis ----
analyze_dose_response <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)

  cat("=== DOSE-RESPONSE ANALYSIS ===\n\n")

  cat("--- Logistic Regression: log(Dose) -> Dysgeusia Incidence ---\n")
  lr <- glm(any_dysgeusia ~ log(dose_mg), data = pl, family = binomial)
  print(tidy(lr, exponentiate = TRUE, conf.int = TRUE), digits = 3)
  cat("\n")

  cat("--- Spearman Correlation: Dose vs Maximum Grade ---\n")
  ct <- cor.test(pl$dose_mg, pl$max_grade, method = "spearman", exact = FALSE)
  cat(sprintf("  rho = %.3f,  p-value = %.4f\n\n", ct$estimate, ct$p.value))

  cat("--- ANOVA: Max Grade across Dose Groups ---\n")
  aov_fit <- aov(max_grade ~ factor(dose_mg), data = pl)
  print(summary(aov_fit))
}

# ---- 3. Cu/Zn vs Dysgeusia Grade ----
analyze_cuzn_grade <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)

  cat("=== CU/ZN RATIO vs DYSGEUSIA GRADE ===\n\n")

  cat("--- Spearman Correlation: Cu/Zn vs Max Grade ---\n")
  ct <- cor.test(pl$cu_zn_ratio, pl$max_grade, method = "spearman", exact = FALSE)
  cat(sprintf("  rho = %.3f,  p-value = %.4f\n\n", ct$estimate, ct$p.value))

  cat("--- Simple Linear Regression: Cu/Zn -> Max Grade ---\n")
  lm1 <- lm(max_grade ~ cu_zn_ratio, data = pl)
  print(tidy(lm1, conf.int = TRUE), digits = 3)
  cat(sprintf("  R² = %.3f\n\n", summary(lm1)$r.squared))

  cat("--- Multiple Regression: log(Dose) + Cu/Zn -> Max Grade ---\n")
  lm2 <- lm(max_grade ~ log(dose_mg) + cu_zn_ratio, data = pl)
  print(tidy(lm2, conf.int = TRUE), digits = 3)
  cat(sprintf("  R² = %.3f\n\n", summary(lm2)$r.squared))

  cat("Note: Cu/Zn ratio is collinear with dose; interpret multiple regression with caution.\n")
}

# ---- 4. Survival / Time to Dysgeusia ----
run_cox_analysis <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)

  cat("=== TIME TO FIRST DYSGEUSIA — COX REGRESSION ===\n\n")
  cat(sprintf("Events: %d / %d patients  (censored at %d days)\n\n",
              sum(pl$surv_event), nrow(pl), CENSORING_DAY))

  cox_fit <- coxph(
    Surv(surv_time, surv_event) ~ log(dose_mg) + cu_zn_ratio,
    data = pl
  )

  cat("--- Hazard Ratios ---\n")
  print(tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE), digits = 3)

  cat(sprintf("\nConcordance (C-index): %.3f\n", cox_fit$concordance["concordance"]))
  cat("Interpretation: HR > 1 means higher exposure -> faster dysgeusia onset.\n")
}

# ---- Plot: Dose-Response Bar + Line ----
plot_dose_response <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)

  summary_df <- pl %>%
    group_by(dose_mg) %>%
    summarise(
      incidence  = mean(any_dysgeusia) * 100,
      mean_grade = mean(max_grade),
      .groups    = "drop"
    )

  p <- ggplot(summary_df, aes(x = factor(dose_mg))) +
    geom_col(aes(y = incidence), fill = "#1976D2", alpha = 0.85, width = 0.6) +
    geom_line(aes(y = mean_grade * 20, group = 1), color = "#F44336", linewidth = 1.2) +
    geom_point(aes(y = mean_grade * 20), color = "#F44336", size = 3.5) +
    scale_y_continuous(
      name     = "Dysgeusia Incidence (%)",
      limits   = c(0, 100),
      sec.axis = sec_axis(~ . / 20, name = "Mean Maximum Grade", breaks = 0:4)
    ) +
    labs(
      x       = "Dose Level (mg)",
      title   = "Dose-Response: Incidence (bars) & Mean Maximum Grade (red line)",
      caption = "Blue bars = incidence %; Red line = mean max CTCAE grade"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggplotly(p)
}

# ---- Plot: Cu/Zn vs Max Grade scatter ----
plot_cuzn_grade <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)

  p <- ggplot(pl,
              aes(x    = cu_zn_ratio,
                  y    = max_grade,
                  color = factor(dose_mg),
                  text  = paste0("Patient: ", patient_id,
                                 "<br>Dose: ", dose_mg, " mg",
                                 "<br>Cu/Zn: ", round(cu_zn_ratio, 2),
                                 "<br>Max Grade: ", max_grade))) +
    geom_jitter(size = 2.8, alpha = 0.70, height = 0.08, width = 0) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_color_brewer(palette = "Set1", name = "Dose (mg)") +
    scale_y_continuous(breaks = 0:4) +
    labs(
      x       = "Baseline Cu/Zn Ratio",
      y       = "Maximum Dysgeusia Grade",
      title   = "Cu/Zn Ratio vs Maximum Dysgeusia Grade",
      caption = "Dashed line = overall linear trend"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggplotly(p, tooltip = "text")
}

# ---- Plot: KM curves by dose ----
plot_km <- function(full_data, dose_filter = c(1, 2, 4, 8, 15, 30)) {
  pl <- patient_level(full_data, dose_filter)
  dose_levels_used <- sort(as.numeric(dose_filter))

  km_list <- lapply(dose_levels_used, function(d) {
    d_sub <- pl %>% filter(dose_mg == d)
    if (nrow(d_sub) == 0) return(NULL)
    fit <- survfit(Surv(surv_time, surv_event) ~ 1, data = d_sub)
    data.frame(
      time       = c(0, fit$time),
      cum_inc    = c(0, 1 - fit$surv),
      dose_group = paste0(d, " mg")
    )
  })

  km_data <- do.call(rbind, Filter(Negate(is.null), km_list))

  if (nrow(km_data) == 0) {
    return(plotly::plot_ly() %>% plotly::layout(title = "No data available"))
  }

  km_data$dose_group <- factor(km_data$dose_group,
                                levels = paste0(dose_levels_used, " mg"))

  p <- ggplot(km_data, aes(x = time, y = cum_inc * 100, color = dose_group)) +
    geom_step(linewidth = 1.0) +
    scale_color_brewer(palette = "Set1", name = "Dose Group") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(
      x       = "Days from Treatment Start",
      y       = "Cumulative Dysgeusia Incidence (%)",
      title   = "Time to First Dysgeusia by Dose Group",
      caption = paste0("Patients without events censored at day ", CENSORING_DAY)
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggplotly(p)
}
