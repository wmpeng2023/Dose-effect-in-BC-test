# app.R — Dose Effect on Safety & Efficacy in Breast Cancer Patients
# Shiny application for analyzing dysgeusia toxicity in relation to
# dose level and Copper/Zinc ratio.
#
# Usage: shiny::runApp("app.R")
# Required packages: shiny, bslib, plotly, dplyr, DT, ggplot2,
#                    survival, broom, officer, flextable

library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(DT)
library(ggplot2)

source("R/simulate_data.R")
source("R/bubble_plot.R")
source("R/analysis.R")
source("R/report_utils.R")

DOSE_LEVELS <- c(1, 2, 4, 8, 15, 30)

# Generate default data once at startup
default_data <- simulate_bc_data(n_patients = 100, seed = 42)

# Save a copy of the simulated data for reference
if (!file.exists("data/simulated_data.csv")) {
  write.csv(default_data$full_data, "data/simulated_data.csv", row.names = FALSE)
}

# ============================================================
# UI
# ============================================================
ui <- page_navbar(
  title = "Dose Effect: Safety & Efficacy in Breast Cancer",
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#1976D2"
  ),
  includeCSS("www/custom.css"),

  # ---- Tab 1: Data Management ----
  nav_panel(
    title = "Data",
    icon  = icon("table"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h4("Data Source"),
        radioButtons(
          "data_source", NULL,
          choices  = c("Simulated Data" = "simulated", "Upload CSV" = "upload"),
          selected = "simulated"
        ),

        conditionalPanel(
          condition = "input.data_source == 'simulated'",
          hr(),
          numericInput("n_patients", "Number of Patients:", value = 100, min = 20, max = 500, step = 10),
          numericInput("sim_seed",   "Random Seed:", value = 42, min = 1),
          actionButton("regenerate", "Regenerate Data", class = "btn-primary btn-sm", width = "100%")
        ),

        conditionalPanel(
          condition = "input.data_source == 'upload'",
          hr(),
          fileInput("upload_file", "Upload CSV:", accept = ".csv"),
          downloadButton("download_template", "Download Template", class = "btn-sm btn-secondary")
        ),

        hr(),
        h4("Export"),
        downloadButton("download_data", "Download Data (CSV)", class = "btn-sm btn-outline-primary")
      ),

      # Main panel: data table
      card(
        card_header("Patient Dataset"),
        DTOutput("data_table")
      )
    )
  ),

  # ---- Tab 2: Bubble Plot ----
  nav_panel(
    title = "Bubble Plot",
    icon  = icon("chart-bubble"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260,
        h4("Filter by Dose"),
        checkboxGroupInput(
          "filter_dose", NULL,
          choices  = paste0(DOSE_LEVELS, " mg"),
          selected = paste0(DOSE_LEVELS, " mg")
        ),
        hr(),
        h4("Filter by Grade"),
        checkboxGroupInput(
          "filter_grade", NULL,
          choices  = paste0("Grade ", 1:4),
          selected = paste0("Grade ", 1:4)
        ),
        hr(),
        h4("Time Range (days)"),
        sliderInput("time_range", NULL, min = 0, max = 365, value = c(0, 180)),
        hr(),
        h4("Display"),
        sliderInput("bubble_scale", "Bubble Size:", min = 1, max = 10, value = 4, step = 0.5),
        hr(),
        downloadButton("download_plot", "Save Plot (PNG)", class = "btn-sm btn-outline-primary")
      ),
      card(
        card_header("Dysgeusia Onset: Cu/Zn Ratio vs Time"),
        plotlyOutput("bubble_plot", height = "580px")
      )
    )
  ),

  # ---- Tab 3: Analysis ----
  nav_panel(
    title = "Analysis",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 260,
        h4("Include Dose Levels"),
        checkboxGroupInput(
          "analysis_dose", NULL,
          choices  = paste0(DOSE_LEVELS, " mg"),
          selected = paste0(DOSE_LEVELS, " mg")
        ),
        hr(),
        h4("Export"),
        downloadButton("download_report", "Export Report (Word)",
                       class = "btn-primary btn-sm", style = "width:100%")
      ),
      navset_card_tab(
        nav_panel("Summary Statistics",
          verbatimTextOutput("summary_stats")
        ),
        nav_panel("Dose-Response",
          plotlyOutput("dose_response_plot", height = "380px"),
          hr(),
          verbatimTextOutput("dose_response_text")
        ),
        nav_panel("Cu/Zn vs Grade",
          plotlyOutput("cuzn_grade_plot", height = "380px"),
          hr(),
          verbatimTextOutput("cuzn_analysis_text")
        ),
        nav_panel("Time to Dysgeusia",
          plotlyOutput("km_plot", height = "380px"),
          hr(),
          verbatimTextOutput("cox_text")
        )
      )
    )
  ),

  # ---- Tab 4: About ----
  nav_panel(
    title = "About",
    icon  = icon("circle-info"),
    card(
      card_header("About This Application"),
      card_body(
        h4("Dose Effect on Safety & Efficacy in Breast Cancer Patients"),
        p("This Shiny application analyzes the relationship between dose level, serum
           Copper/Zinc (Cu/Zn) ratio, and dysgeusia (taste disturbance) toxicity in
           breast cancer patients."),
        h4("Hypothesis"),
        tags$blockquote(
          "Higher dose → Higher Cu/Zn ratio → Earlier onset and greater severity of dysgeusia"
        ),
        h4("Variables"),
        tags$ul(
          tags$li(tags$b("Dose Level:"), " 1, 2, 4, 8, 15, 30 mg"),
          tags$li(tags$b("Cu/Zn Ratio:"), " Single baseline serum measurement (cancer patients: ~1.5–3.5)"),
          tags$li(tags$b("Dysgeusia Grade:"), " CTCAE v5 grading (0 = none, 1 = mild, 2 = moderate, 3 = severe, 4 = life-threatening)"),
          tags$li(tags$b("Time to Dysgeusia:"), " Days from treatment start to first occurrence of each grade")
        ),
        h4("Data"),
        p("Default data is simulated with realistic clinical parameters. Use the ",
          tags$b("Data"), " tab to upload your own dataset. Download the CSV template
          for the required column format."),
        h4("Reference: Cu/Zn in Breast Cancer"),
        p("Serum Cu/Zn ratio is elevated in cancer patients (typically 1.5–3.0+)
           compared to healthy adults (~0.7–1.3). Higher ratios reflect increased
           oxidative stress and may predict treatment-related toxicity."),
        hr(),
        p(em("Built with R Shiny | GitHub: Wmpeng2023/Dose effect in BC (test)"))
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: current dataset ----
  current_data <- reactiveVal(default_data)

  # Regenerate simulated data
  observeEvent(input$regenerate, {
    withProgress(message = "Simulating data...", {
      new_data <- simulate_bc_data(
        n_patients = input$n_patients,
        seed       = input$sim_seed
      )
      current_data(new_data)
    })
    showNotification("Data regenerated successfully.", type = "message", duration = 3)
  })

  # Upload custom data
  observeEvent(input$upload_file, {
    req(input$upload_file)
    tryCatch({
      uploaded <- read.csv(input$upload_file$datapath, stringsAsFactors = FALSE)

      required_cols <- c("patient_id", "dose_mg", "cu_zn_ratio", "dysgeusia_grade")
      missing_cols  <- setdiff(required_cols, names(uploaded))

      if (length(missing_cols) > 0) {
        showNotification(
          paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
          type = "error", duration = 6
        )
        return()
      }

      # Build data_list structure from uploaded CSV
      patient_data <- uploaded %>%
        select(patient_id, dose_mg, cu_zn_ratio) %>%
        distinct()

      event_data <- uploaded %>%
        filter(dysgeusia_grade > 0) %>%
        select(patient_id, dysgeusia_grade, any_of("time_to_first_onset"))

      # Add dose_label if missing
      if (!"dose_label" %in% names(uploaded)) {
        all_doses <- sort(unique(uploaded$dose_mg))
        uploaded$dose_label <- factor(
          paste0(uploaded$dose_mg, " mg"),
          levels = paste0(all_doses, " mg")
        )
      }

      current_data(list(
        patient_data = patient_data,
        event_data   = event_data,
        full_data    = uploaded
      ))
      showNotification("Data uploaded successfully.", type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 6)
    })
  })

  # ---- Helper: parse filter inputs ----
  selected_doses <- reactive({
    as.numeric(gsub(" mg", "", input$filter_dose))
  })

  analysis_doses <- reactive({
    as.numeric(gsub(" mg", "", input$analysis_dose))
  })

  selected_grades <- reactive({
    as.integer(gsub("Grade ", "", input$filter_grade))
  })

  # ---- Filtered data for bubble plot ----
  filtered_data <- reactive({
    fd <- current_data()$full_data

    fd %>%
      filter(
        dose_mg         %in% selected_doses(),
        dysgeusia_grade %in% selected_grades(),
        is.na(time_to_first_onset) |
          (time_to_first_onset >= input$time_range[1] &
             time_to_first_onset <= input$time_range[2])
      )
  })

  # ---- Data table ----
  output$data_table <- renderDT({
    datatable(
      current_data()$full_data,
      options = list(
        pageLength = 15,
        scrollX    = TRUE,
        dom        = "Bfrtip"
      ),
      rownames = FALSE,
      filter   = "top"
    )
  })

  # ---- Bubble plot ----
  output$bubble_plot <- renderPlotly({
    create_bubble_plot(filtered_data(), bubble_scale = input$bubble_scale)
  })

  # ---- Analysis outputs ----
  output$summary_stats <- renderPrint({
    summarize_data(current_data(), analysis_doses())
  })

  output$dose_response_plot <- renderPlotly({
    plot_dose_response(current_data()$full_data, analysis_doses())
  })

  output$dose_response_text <- renderPrint({
    analyze_dose_response(current_data()$full_data, analysis_doses())
  })

  output$cuzn_grade_plot <- renderPlotly({
    plot_cuzn_grade(current_data()$full_data, analysis_doses())
  })

  output$cuzn_analysis_text <- renderPrint({
    analyze_cuzn_grade(current_data()$full_data, analysis_doses())
  })

  output$km_plot <- renderPlotly({
    plot_km(current_data()$full_data, analysis_doses())
  })

  output$cox_text <- renderPrint({
    run_cox_analysis(current_data()$full_data, analysis_doses())
  })

  # ---- Downloads ----
  output$download_data <- downloadHandler(
    filename = function() paste0("bc_dose_dysgeusia_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(current_data()$full_data, file, row.names = FALSE)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("bubble_plot_", Sys.Date(), ".png"),
    content  = function(file) {
      p <- create_bubble_plot_static(filtered_data(), bubble_scale = input$bubble_scale)
      ggsave(file, p, width = 12, height = 7, dpi = 300)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() paste0("dysgeusia_report_", Sys.Date(), ".docx"),
    content  = function(file) {
      withProgress(message = "Generating Word report...", {
        export_word_report(current_data(), analysis_doses(), file)
      })
    }
  )

  output$download_template <- downloadHandler(
    filename = "data_upload_template.csv",
    content  = function(file) {
      template <- data.frame(
        patient_id          = c(1L, 1L, 2L, 3L),
        dose_mg             = c(1,  1,  2,  4),
        cu_zn_ratio         = c(1.8, 1.8, 2.1, 1.6),
        dysgeusia_grade     = c(1L, 2L, 1L, 0L),
        time_to_first_onset = c(5, 12, 8, NA)
      )
      write.csv(template, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
