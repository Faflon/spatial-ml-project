# =============================================================================
# 00_setup.R — Package Installation & Project Configuration
# Spatial ML Project: Air Pollution in Poland
# Authors: Ondřej Marvan & Adam Jaworski (WNE UW, 2026)
# =============================================================================

# --- 1. Install required packages -------------------------------------------

required_packages <- c(
  # Spatial core
  "sf", "terra", "stars",
  # Google Earth Engine
  "rgee", "reticulate",
  # Data APIs
  "httr2", "jsonlite",
  # OSM data
  "osmdata",
  # Data wrangling
  "tidyverse", "lubridate", "glue",
  # Spatial analysis (will need later)
  "spdep", "spatialreg", "ranger", "caret",
  # Visualization
  "ggplot2", "tmap", "viridis"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing ", pkg, "..."))
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))

# bdl package (GUS BDL API wrapper) — install separately if needed:
# install.packages("bdl")
# If not on CRAN for your R version, the API can be called directly via httr2
# (see 03_collect_gus.R for both approaches)

# --- 2. Project configuration -----------------------------------------------

# All paths relative to project root
CONFIG <- list(
  # CRS
  crs_pl     = 2180,        # EPSG:2180 — Poland CS92 (metric, for spatial ops)
  crs_wgs84  = 4326,        # EPSG:4326 — WGS84 (for API queries & GEE)
  
  # Temporal scope — pick a full calendar year for clean analysis
  date_start = "2024-01-01",
  date_end   = "2024-12-31",
  
  # Spatial scope — Poland bounding box (WGS84)
  # NOTE: base::c() used explicitly because terra overrides c() after loading
  bbox_poland = base::c(xmin = 14.07, ymin = 49.00, xmax = 24.15, ymax = 54.84),
  
  # Grid resolution (meters) — roughly matching TROPOMI pixel size
  
  grid_res = 5000,
  
  # Output directories
  dir_raw    = "data/raw",
  dir_processed = "data/processed",
  dir_output = "data/output",
  dir_figures = "figures"
)

# Create directory structure
dirs <- base::c(
  CONFIG$dir_raw,
  file.path(CONFIG$dir_raw, "tropomi"),
  file.path(CONFIG$dir_raw, "gios"),
  file.path(CONFIG$dir_raw, "gus"),
  file.path(CONFIG$dir_raw, "osm"),
  CONFIG$dir_processed,
  CONFIG$dir_output,
  CONFIG$dir_figures
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Setup complete. Directory structure created.")
message("Temporal scope: ", CONFIG$date_start, " to ", CONFIG$date_end)

# --- 3. Google Earth Engine initialization -----------------------------------
# NOTE: First-time setup requires:
#   1. A Google account with GEE access (signup: https://earthengine.google.com/)
#   2. Python environment with 'earthengine-api' installed
#
# Run once:
#   rgee::ee_install()          # installs Python env
#   rgee::ee_Initialize()       # authenticates with Google
#
# If you have issues with rgee, the fallback is downloading GeoTIFFs
# from the S5P-PAL portal: https://maps.s5p-pal.com/no2/
# =============================================================================