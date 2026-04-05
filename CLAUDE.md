# CLAUDE.md — Dose Effect on Safety & Efficacy in Breast Cancer Patients

## Project Overview
An R Shiny application analyzing the relationship between dose level, Copper/Zinc (Cu/Zn) ratio, and dysgeusia toxicity in breast cancer patients.

## Research Question
Does higher dose level lead to higher Cu/Zn ratio, and does that drive earlier onset and greater severity of dysgeusia (taste disturbance)?

## Hypothesis
Higher dose → higher Cu/Zn ratio → earlier onset + higher CTCAE grade of dysgeusia

## Data Structure

### Patient Table (`patient_data`)
One row per patient:
| Column | Type | Description |
|--------|------|-------------|
| patient_id | integer | Unique patient identifier |
| dose_mg | numeric | Dose level: 1, 2, 4, 8, 15, or 30 mg |
| cu_zn_ratio | numeric | Baseline serum Cu/Zn ratio |

### Event Table (`event_data`) — Long Format
One row per patient × grade:
| Column | Type | Description |
|--------|------|-------------|
| patient_id | integer | Links to patient table |
| dysgeusia_grade | integer | CTCAE grade (1–4) at first occurrence |
| time_to_first_onset | integer | Days from treatment start to first report of this grade |

### Full/Merged Table (`full_data`)
- Event patients: multiple rows (one per grade observed, grade ≥ 1)
- Event-free patients: one row with dysgeusia_grade = 0, time_to_first_onset = NA

## Key Rules
- A patient can have grade 1 without grade 2, but NOT grade 2 without grade 1
- `time_to_first_onset` is the FIRST occurrence of each grade (not a repeat)
- Grade 0 = no dysgeusia (event-free / censored)
- Censoring time = 180 days for survival analysis

## File Structure
```
Project 1/
├── CLAUDE.md               ← This file
├── .gitignore
├── app.R                   ← Main Shiny entry point
├── R/
│   ├── simulate_data.R     ← Data simulation logic
│   ├── bubble_plot.R       ← Bubble plot functions (ggplot2 + plotly)
│   ├── analysis.R          ← Statistical analysis functions
│   └── report_utils.R      ← Word export / reporting
├── data/
│   └── simulated_data.csv  ← Generated on app start (auto-saved)
├── reports/
│   └── summary_report.Rmd  ← Rmarkdown report template
└── www/
    └── custom.css          ← App styling
```

## R Package Dependencies
- **Shiny ecosystem**: `shiny`, `bslib`, `DT`, `shinyWidgets`
- **Visualization**: `ggplot2`, `plotly`
- **Data**: `dplyr`, `tidyr`
- **Statistics**: `survival`, `broom`
- **Reporting**: `officer`, `flextable`
- **Optional**: `survminer`

## Shiny App Architecture
- **Data tab**: View/upload data, download template, regenerate simulation
- **Bubble Plot tab**: Interactive plot with dose/grade filters
- **Analysis tab**: Summary stats, dose-response, Cu/Zn vs grade, survival analysis
- **Report**: Export to Word from Analysis tab

## Clinical Reference: Cu/Zn Ratio in Cancer
- Healthy adults: ~0.7–1.3
- Breast cancer patients: typically 1.5–3.0+
- Higher Cu/Zn associated with oxidative stress and treatment toxicity
- Source: Published literature on serum copper/zinc in oncology

## GitHub
- Repository: `Dose effect in BC (test)` (private)
- Owner: Wmpeng2023
