# =============================================================================
# 10_utils.R — Shared Utility Functions
# =============================================================================

library(sf)
library(tidyverse)
library(glue)

# --- Theme for consistent plots ----------------------------------------------

theme_spatial <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey40"),
      legend.position = "right",
      panel.grid    = element_line(color = "grey95")
    )
}

# --- Safe file reader (check existence) --------------------------------------

safe_read_sf <- function(path, ...) {
  if (!file.exists(path)) {
    warning(glue("File not found: {path}"))
    return(NULL)
  }
  st_read(path, quiet = TRUE, ...)
}

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) {
    warning(glue("File not found: {path}"))
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE, ...)
}

# --- Accuracy metrics --------------------------------------------------------

calc_metrics <- function(actual, predicted) {
  #' Compute standard regression metrics.
  valid <- !is.na(actual) & !is.na(predicted)
  a <- actual[valid]
  p <- predicted[valid]
  
  tibble(
    n    = length(a),
    R2   = cor(a, p)^2,
    RMSE = sqrt(mean((a - p)^2)),
    MAE  = mean(abs(a - p)),
    MAPE = mean(abs((a - p) / a)) * 100,
    ME   = mean(p - a)  # bias
  )
}

# --- Grid helpers ------------------------------------------------------------

grid_centroid_coords <- function(grid_sf) {
  #' Extract centroid coordinates as x, y columns.
  coords <- st_coordinates(st_centroid(grid_sf))
  grid_sf$x <- coords[, 1]
  grid_sf$y <- coords[, 2]
  grid_sf
}

# --- Pretty printing ---------------------------------------------------------

print_section <- function(title) {
  width <- 70
  line <- paste(rep("=", width), collapse = "")
  message(glue("\n{line}\n  {title}\n{line}\n"))
}

# =============================================================================