# Dose Effect on Safety & Efficacy in Breast Cancer Patients

An R Shiny application analyzing the relationship between dose level, serum Copper/Zinc (Cu/Zn) ratio, and dysgeusia (taste disturbance) toxicity in breast cancer patients.

## Research Question

> Does higher dose lead to higher Cu/Zn ratio, and does that drive earlier onset and greater severity of dysgeusia?

## App Features

- **Simulated dataset** — 100 patients across 6 dose levels (1, 2, 4, 8, 15, 30 mg)
- **Interactive bubble plot** — time vs Cu/Zn ratio, bubble size/color by CTCAE dysgeusia grade
- **Dose & grade filters** — explore subsets interactively
- **Statistical analysis** — dose-response, Cu/Zn vs grade, Cox regression, KM curves
- **Upload your own data** — CSV template provided in the app
- **Export to Word** — download a formatted summary report

## Access Options

### Option A — Live Web App (no installation needed)
> **Click to open:** `https://YOUR_SHINYAPPS_USERNAME.shinyapps.io/Dose-effect-in-BC-test/`
>
> *(Link will be active after deployment — see [Deploying to shinyapps.io](#deploying-to-shinyappsio) below)*

Anyone with the link can open and use the app in their browser — no R required.

### Option B — Run Locally (requires R)
Follow the Quick Start steps below.

---

## Requirements (Option B only)

- R >= 4.1.0 ([download](https://cran.r-project.org/))
- RStudio (recommended — [download](https://posit.co/download/rstudio-desktop/))

## Quick Start

### Step 1 — Clone the repository

```bash
git clone https://github.com/wmpeng2023/Dose-effect-in-BC-test.git
cd Dose-effect-in-BC-test
```

### Step 2 — Restore all R packages (one time only)

Open the project in RStudio, then run in the R console:

```r
install.packages("renv")   # if not already installed
renv::restore()            # installs exact package versions from renv.lock
```

This will install all dependencies automatically. Takes a few minutes on first run.

### Step 3 — Launch the app

```r
shiny::runApp()
```

The app will open in your browser at `http://127.0.0.1:XXXX`.

## Data Structure

The app uses a **long-format** dataset where each row is one patient × dysgeusia grade:

| Column | Type | Description |
|--------|------|-------------|
| `patient_id` | integer | Unique patient ID |
| `dose_mg` | numeric | Dose level (1, 2, 4, 8, 15, or 30 mg) |
| `cu_zn_ratio` | numeric | Baseline serum Cu/Zn ratio |
| `dysgeusia_grade` | integer | CTCAE grade (0 = none, 1–4 = severity) |
| `time_to_first_onset` | numeric | Days from treatment start to first occurrence of this grade (NA if grade 0) |

- Event patients have **multiple rows** (one per grade experienced)
- Event-free patients have **one row** with grade = 0 and time = NA
- Download the CSV template from the **Data** tab to upload your own dataset

## Deploying to shinyapps.io

To publish the app so anyone can open it via a browser link:

1. Sign up for a free account at [shinyapps.io](https://www.shinyapps.io)
2. Go to **Account → Tokens → Show** and copy your `name`, `token`, and `secret`
3. Open `deploy.R`, fill in your credentials, and run it in RStudio:

```r
source("deploy.R")
```

Your app will be live at:
```
https://YOUR_USERNAME.shinyapps.io/Dose-effect-in-BC-test/
```

Share that URL — anyone can click it, no R needed.

> **Note:** The free tier allows 25 active hours/month and up to 5 apps.
> Credentials in `deploy.R` are **never committed to GitHub** (excluded via `.gitignore`).

---

## File Structure

```
├── app.R                  # Main Shiny entry point
├── deploy.R               # shinyapps.io deployment script
├── DESCRIPTION            # Package dependency declarations
├── renv.lock              # Exact package versions (used by renv::restore())
├── R/
│   ├── simulate_data.R    # Patient data simulation
│   ├── bubble_plot.R      # Bubble plot (plotly + ggplot2)
│   ├── analysis.R         # Statistical analyses
│   └── report_utils.R     # Word report export
├── www/
│   └── custom.css         # App styling
└── data/                  # Auto-generated simulated dataset
```

## Clinical Background

Serum Cu/Zn ratio is elevated in breast cancer patients (~1.5–3.0+) compared to healthy adults (~0.7–1.3). Higher ratios are associated with increased oxidative stress and may predict treatment-related toxicity such as dysgeusia.

## License

MIT
