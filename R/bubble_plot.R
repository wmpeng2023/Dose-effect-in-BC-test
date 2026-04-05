# bubble_plot.R
# Creates interactive and static bubble plots:
#   X: time from treatment start (days)
#   Y: Cu/Zn ratio
#   Bubble size: CTCAE dysgeusia grade
#   Bubble color: CTCAE dysgeusia grade

library(ggplot2)
library(plotly)
library(dplyr)

# Grade color palette (colorblind-friendly)
GRADE_COLORS <- c(
  "Grade 0" = "#AAAAAA",
  "Grade 1" = "#4CAF50",
  "Grade 2" = "#FF9800",
  "Grade 3" = "#F44336",
  "Grade 4" = "#9C27B0"
)

# ---- Interactive plotly bubble plot (for Shiny) ----
create_bubble_plot <- function(data, bubble_scale = 5) {
  if (is.null(data) || nrow(data) == 0) {
    return(
      plotly::plot_ly() %>%
        plotly::layout(title = "No data to display with current filters")
    )
  }

  plot_data <- data %>%
    filter(!is.na(time_to_first_onset), dysgeusia_grade > 0) %>%
    mutate(
      grade_label = factor(
        paste0("Grade ", dysgeusia_grade),
        levels = paste0("Grade ", 1:4)
      ),
      bubble_size = dysgeusia_grade * bubble_scale,
      tooltip_text = paste0(
        "<b>Patient:</b> ", patient_id, "<br>",
        "<b>Dose:</b> ", dose_mg, " mg<br>",
        "<b>Cu/Zn Ratio:</b> ", round(cu_zn_ratio, 2), "<br>",
        "<b>Dysgeusia Grade:</b> ", dysgeusia_grade, "<br>",
        "<b>Day of Onset:</b> ", time_to_first_onset
      )
    )

  if (nrow(plot_data) == 0) {
    return(
      plotly::plot_ly() %>%
        plotly::layout(title = "No dysgeusia events in selected filters")
    )
  }

  p <- ggplot(plot_data,
              aes(x     = time_to_first_onset,
                  y     = cu_zn_ratio,
                  size  = dysgeusia_grade,
                  color = grade_label,
                  text  = tooltip_text)) +
    geom_point(alpha = 0.78) +
    scale_size_continuous(
      range  = c(bubble_scale * 0.7, bubble_scale * 2.8),
      name   = "Grade",
      breaks = 1:4
    ) +
    scale_color_manual(
      values = GRADE_COLORS[paste0("Grade ", 1:4)],
      name   = "Dysgeusia Grade",
      drop   = FALSE
    ) +
    labs(
      title    = "Dysgeusia Onset: Cu/Zn Ratio vs Time from Treatment Start",
      x        = "Time from Treatment Start (Days)",
      y        = "Serum Copper/Zinc Ratio",
      caption  = "Bubble size and color reflect CTCAE dysgeusia grade (1 = mild, 4 = severe)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position  = "right",
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggplotly(p, tooltip = "text") %>%
    layout(legend = list(orientation = "v"))
}

# ---- Static ggplot version (for download / Word export) ----
create_bubble_plot_static <- function(data, bubble_scale = 5) {
  plot_data <- data %>%
    filter(!is.na(time_to_first_onset), dysgeusia_grade > 0) %>%
    mutate(
      grade_label = factor(
        paste0("Grade ", dysgeusia_grade),
        levels = paste0("Grade ", 1:4)
      )
    )

  if (nrow(plot_data) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No dysgeusia events to display") +
        theme_void()
    )
  }

  ggplot(plot_data,
         aes(x     = time_to_first_onset,
             y     = cu_zn_ratio,
             size  = dysgeusia_grade,
             color = grade_label)) +
    geom_point(alpha = 0.78) +
    scale_size_continuous(
      range  = c(bubble_scale * 0.7, bubble_scale * 2.8),
      name   = "Grade",
      breaks = 1:4
    ) +
    scale_color_manual(
      values = GRADE_COLORS[paste0("Grade ", 1:4)],
      name   = "Dysgeusia Grade",
      drop   = FALSE
    ) +
    labs(
      title   = "Dysgeusia Onset: Cu/Zn Ratio vs Time from Treatment Start",
      x       = "Time from Treatment Start (Days)",
      y       = "Serum Copper/Zinc Ratio",
      caption = "Bubble size and color reflect CTCAE dysgeusia grade"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "right",
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}
