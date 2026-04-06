# =============================================================================
# 00_setup.R — Package Installation & Project Configuration
# Spatial ML Project: Air Pollution in Poland
# Authors: Ondřej Marvan & Adam Jaworski (WNE UW, 2026)
# =============================================================================

# --- 1. Install required packages -------------------------------------------
# We avoid loading packages here to prevent side effects (e.g. terra
# overriding base::c). Each script loads only what it needs.

# CORE packages (from CRAN) — needed for data collection & basic analysis
core_packages <- base::c(
  "sf", "terra", "stars",
  "httr2", "jsonlite", "osmdata",
  "tidyverse", "lubridate", "glue",
  "viridis", "corrplot"
)

# ANALYSIS packages (from CRAN) — needed for modelling & kriging
analysis_packages <- base::c(
  "spdep", "spatialreg", "ranger", "caret", "nnet",
  "gstat", "automap",
  "dbscan", "cluster",
  "rgee", "reticulate"
)

install_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  message(paste0("Installing ", pkg, "..."))
  tryCatch(
    utils::install.packages(pkg, dependencies = TRUE),
    warning = function(w) message(paste0("  WARNING: ", w$message)),
    error   = function(e) message(paste0("  FAILED: ", e$message))
  )
  # Check if it actually installed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("  >>> ", pkg, " not available. Install manually or skip."))
  }
}

message("--- Installing core packages ---")
for (pkg in core_packages) install_pkg(pkg)

message("\n--- Installing analysis packages ---")
for (pkg in analysis_packages) install_pkg(pkg)

# SPECIAL INSTALLS (not on CRAN for all R versions):
#
# SpatialML (geographically weighted random forest):
#   install.packages("remotes")
#   remotes::install_github("SpatialML/SpatialML")
#
# tmap (optional, for interactive maps — has many dependencies):
#   install.packages("tmap")
#
# If terra fails due to network, retry:
#   install.packages("terra", repos = "https://cran.rstudio.com/")

# Quick status check
installed <- sapply(base::c(core_packages, analysis_packages), function(p) {
  requireNamespace(p, quietly = TRUE)
})
n_ok <- sum(installed)
n_fail <- sum(!installed)
message(paste0("\nPackage status: ", n_ok, " installed, ", n_fail, " missing"))
if (n_fail > 0) {
  message("Missing: ", paste(names(installed[!installed]), collapse = ", "))
  message("Some may be optional. Re-run or install manually as needed.")
}

# --- 2. Project configuration -----------------------------------------------

CONFIG <- base::list(
  # CRS
  crs_pl    = 2180L,
  crs_wgs84 = 4326L,

  # Temporal scope
  date_start = "2024-01-01",
  date_end   = "2024-12-31",

  # Poland bounding box (WGS84) — using base::c explicitly
  bbox = base::c(xmin = 14.07, ymin = 49.00, xmax = 24.15, ymax = 54.84),

  # Grid resolution (meters)
  grid_res = 5000L,

  # Directories (relative to project root)
  dir_raw       = "data/raw",
  dir_processed = "data/processed",
  dir_output    = "data/output",
  dir_figures   = "figures"
)

# Create directories
for (d in base::c(
  CONFIG$dir_raw,
  file.path(CONFIG$dir_raw, "tropomi"),
  file.path(CONFIG$dir_raw, "gios"),
  file.path(CONFIG$dir_raw, "gus"),
  file.path(CONFIG$dir_raw, "osm"),
  CONFIG$dir_processed,
  CONFIG$dir_output,
  CONFIG$dir_figures
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- 3. Shared helper functions ----------------------------------------------

#' Safely read an sf file (returns NULL if missing)
safe_read_sf <- function(path, ...) {
  if (!file.exists(path)) { warning(paste0("Not found: ", path)); return(NULL) }
  sf::st_read(path, quiet = TRUE, ...)
}

#' Safely read CSV (returns NULL if missing)
safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) { warning(paste0("Not found: ", path)); return(NULL) }
  readr::read_csv(path, show_col_types = FALSE, ...)
}

#' Regression accuracy metrics
calc_metrics <- function(actual, predicted) {
  ok <- !is.na(actual) & !is.na(predicted)
  a <- actual[ok]; p <- predicted[ok]
  if (length(a) < 3) return(NULL)
  data.frame(
    n = length(a), R2 = cor(a, p)^2,
    RMSE = sqrt(mean((a - p)^2)),
    MAE  = mean(abs(a - p))
  )
}

message("Setup complete.")
message("Year: ", format(as.Date(CONFIG$date_start), "%Y"))
message("Grid: ", CONFIG$grid_res / 1000, " km")

# =============================================================================
# Google Earth Engine — one-time setup (run interactively):
#   rgee::ee_install()
#   rgee::ee_Initialize()
#
# If rgee doesn't work, use S5P-PAL portal: https://maps.s5p-pal.com/no2/
# =============================================================================
