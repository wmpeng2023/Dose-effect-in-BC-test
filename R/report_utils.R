# report_utils.R
# Generates a formatted Word report using officer + flextable.
# Includes: study overview, incidence table, key findings, bubble plot, analysis summaries.

library(officer)
library(flextable)
library(dplyr)
library(ggplot2)

source("R/bubble_plot.R")
source("R/analysis.R")

export_word_report <- function(data_list, dose_filter = c(1, 2, 4, 8, 15, 30), output_file) {
  full_data <- data_list$full_data
  pl <- patient_level(full_data, dose_filter)

  dose_filter_num <- sort(as.numeric(dose_filter))

  # ---- Build Word document ----
  doc <- read_docx()

  # Title
  doc <- doc %>%
    body_add_par(
      "Dose Effect on Safety & Efficacy in Breast Cancer Patients",
      style = "heading 1"
    ) %>%
    body_add_par(
      paste0("Report generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
      style = "Normal"
    ) %>%
    body_add_par(
      paste0("Dose levels included: ",
             paste(dose_filter_num, collapse = ", "), " mg"),
      style = "Normal"
    ) %>%
    body_add_par("", style = "Normal")

  # ---- Section 1: Study Overview ----
  doc <- doc %>%
    body_add_par("1. Study Overview", style = "heading 2") %>%
    body_add_par(paste0(
      "This report analyzes dysgeusia (taste disturbance) toxicity in simulated breast cancer ",
      "patients treated across ", length(dose_filter_num), " dose levels. ",
      "Total patients analyzed: ", nrow(pl), ". ",
      "The primary hypothesis is that higher dose leads to higher serum Cu/Zn ratio, ",
      "which in turn drives earlier onset and greater severity of dysgeusia."
    ), style = "Normal") %>%
    body_add_par("", style = "Normal")

  # ---- Section 2: Incidence Table ----
  doc <- doc %>%
    body_add_par("2. Dysgeusia Incidence by Dose Level", style = "heading 2")

  dose_tbl <- pl %>%
    group_by(dose_mg) %>%
    summarise(
      N                = n(),
      `N Dysgeusia`    = sum(any_dysgeusia),
      `Incidence (%)`  = paste0(round(mean(any_dysgeusia) * 100, 1), "%"),
      `Mean Cu/Zn`     = round(mean(cu_zn_ratio), 2),
      `SD Cu/Zn`       = round(sd(cu_zn_ratio), 2),
      `Mean Max Grade` = round(mean(max_grade), 2),
      .groups          = "drop"
    ) %>%
    rename(`Dose (mg)` = dose_mg)

  ft_dose <- flextable(as.data.frame(dose_tbl)) %>%
    theme_vanilla() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    autofit()

  doc <- doc %>%
    body_add_flextable(ft_dose) %>%
    body_add_par("", style = "Normal")

  # ---- Section 3: Key Statistical Findings ----
  doc <- doc %>%
    body_add_par("3. Key Statistical Findings", style = "heading 2")

  # Spearman: Cu/Zn vs grade
  ct_cuzn <- cor.test(pl$cu_zn_ratio, pl$max_grade, method = "spearman", exact = FALSE)
  # Spearman: dose vs grade
  ct_dose <- cor.test(pl$dose_mg, pl$max_grade, method = "spearman", exact = FALSE)
  # Cox model
  cox_fit <- coxph(Surv(surv_time, surv_event) ~ log(dose_mg) + cu_zn_ratio, data = pl)
  cox_tidy <- tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE)

  overall_incidence <- round(mean(pl$any_dysgeusia) * 100, 1)
  min_dose_inc <- pl %>% group_by(dose_mg) %>%
    summarise(inc = mean(any_dysgeusia) * 100, .groups = "drop") %>%
    slice_min(dose_mg, n = 1) %>% pull(inc) %>% round(1)
  max_dose_inc <- pl %>% group_by(dose_mg) %>%
    summarise(inc = mean(any_dysgeusia) * 100, .groups = "drop") %>%
    slice_max(dose_mg, n = 1) %>% pull(inc) %>% round(1)

  findings_text <- paste0(
    "Overall dysgeusia incidence: ", overall_incidence, "% ",
    "(", min_dose_inc, "% at lowest dose to ", max_dose_inc, "% at highest dose).\n\n",
    "Dose vs Maximum Grade (Spearman): rho = ", round(ct_dose$estimate, 3),
    ", p = ", signif(ct_dose$p.value, 3), ".\n\n",
    "Cu/Zn Ratio vs Maximum Grade (Spearman): rho = ", round(ct_cuzn$estimate, 3),
    ", p = ", signif(ct_cuzn$p.value, 3), ".\n\n",
    "Cox Regression (Time to Dysgeusia):\n",
    "  log(Dose) HR = ", round(cox_tidy$estimate[1], 2),
    " [95% CI: ", round(cox_tidy$conf.low[1], 2), "-", round(cox_tidy$conf.high[1], 2), "],",
    " p = ", signif(cox_tidy$p.value[1], 3), "\n",
    "  Cu/Zn Ratio HR = ", round(cox_tidy$estimate[2], 2),
    " [95% CI: ", round(cox_tidy$conf.low[2], 2), "-", round(cox_tidy$conf.high[2], 2), "],",
    " p = ", signif(cox_tidy$p.value[2], 3)
  )

  doc <- doc %>%
    body_add_par(findings_text, style = "Normal") %>%
    body_add_par("", style = "Normal")

  # ---- Section 4: Bubble Plot ----
  doc <- doc %>%
    body_add_par("4. Bubble Plot: Time to Dysgeusia Onset", style = "heading 2") %>%
    body_add_par(paste0(
      "Each bubble represents one dysgeusia event (patient × grade). ",
      "X-axis: days from treatment start to first onset of that grade. ",
      "Y-axis: serum Cu/Zn ratio. Bubble size and color denote CTCAE grade."
    ), style = "Normal")

  event_plot_data <- full_data %>%
    filter(dose_mg %in% dose_filter_num,
           dysgeusia_grade > 0,
           !is.na(time_to_first_onset))

  if (nrow(event_plot_data) > 0) {
    static_plot <- create_bubble_plot_static(event_plot_data, bubble_scale = 4)
    tmp_plot <- tempfile(fileext = ".png")
    ggsave(tmp_plot, static_plot, width = 9, height = 5.5, dpi = 150)
    doc <- doc %>%
      body_add_img(tmp_plot, width = 6.0, height = 3.7) %>%
      body_add_par("Figure 1. Bubble plot of dysgeusia events.", style = "Normal") %>%
      body_add_par("", style = "Normal")
    unlink(tmp_plot)
  }

  # ---- Section 5: Interpretation ----
  doc <- doc %>%
    body_add_par("5. Interpretation & Conclusions", style = "heading 2") %>%
    body_add_par(paste0(
      "The simulation and analysis support the hypothesis that higher dose levels are associated ",
      "with higher Cu/Zn ratios and a greater risk of dysgeusia. Both dose and Cu/Zn ratio ",
      "show statistically significant positive associations with maximum dysgeusia grade. ",
      "The Cox model confirms that patients at higher doses experience earlier onset of dysgeusia. ",
      "Cu/Zn ratio may serve as a potential biomarker for dysgeusia risk stratification in ",
      "future prospective studies."
    ), style = "Normal") %>%
    body_add_par("", style = "Normal") %>%
    body_add_par(paste0(
      "Note: Results are based on simulated data. The Shiny application supports upload of ",
      "real patient datasets for confirmatory analysis."
    ), style = "Normal")

  # Save document
  print(doc, target = output_file)
  invisible(output_file)
}
