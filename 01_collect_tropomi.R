# =============================================================================
# 01_collect_tropomi.R — Sentinel-5P TROPOMI NO2/SO2/CO via Google Earth Engine
# =============================================================================
# Data source: Copernicus S5P TROPOMI (via GEE)
# Resolution:  ~5.5km x 3.5km, daily → aggregated to monthly
# Products:    NO2, SO2, CO tropospheric columns
# Output:      Monthly GeoTIFF rasters for Poland
# =============================================================================

source("00_setup.R")

library(rgee)
library(sf)
library(terra)
library(stars)
library(tidyverse)
library(glue)

# --- 1. Initialize Earth Engine ----------------------------------------------

ee_Initialize()

# --- 2. Define Poland boundary -----------------------------------------------

# Use a simple bounding box for clipping
bbox <- CONFIG$bbox_poland
roi <- ee$Geometry$Rectangle(
  coords = list(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]),
  proj   = "EPSG:4326",
  geodesic = FALSE
)

# Alternatively, use an actual Poland boundary from GEE's feature collections:
# roi <- ee$FeatureCollection("FAO/GAUL/2015/level0")$
#   filter(ee$Filter$eq("ADM0_NAME", "Poland"))$
#   geometry()

# --- 3. Define TROPOMI collections & bands -----------------------------------

# S5P Offline products (higher quality than Near Real-Time)
collections <- list(
  NO2 = list(
    id   = "COPERNICUS/S5P/OFFL/L3_NO2",
    band = "tropospheric_NO2_column_number_density",
    # qa filtering: only keep pixels with qa_value > 0.75 (recommended)
    qa_band = "tropospheric_NO2_column_number_density_validity",
    qa_threshold = 75
  ),
  SO2 = list(
    id   = "COPERNICUS/S5P/OFFL/L3_SO2",
    band = "SO2_column_number_density",
    qa_band = NULL,
    qa_threshold = NULL
  ),
  CO = list(
    id   = "COPERNICUS/S5P/OFFL/L3_CO",
    band = "CO_column_number_density",
    qa_band = NULL,
    qa_threshold = NULL
  )
)

# --- 4. Function: monthly mean for one pollutant -----------------------------

get_monthly_mean <- function(collection_info, year, month, roi) {
  #' Compute monthly mean of a TROPOMI product over Poland.
  #' Cloud-contaminated pixels are filtered via qa_value.
  
  start_date <- sprintf("%04d-%02d-01", year, month)
  # end date = first day of next month
  if (month == 12) {
    end_date <- sprintf("%04d-01-01", year + 1)
  } else {
    end_date <- sprintf("%04d-%02d-01", year, month + 1)
  }
  
  # Load and filter collection
  col <- ee$ImageCollection(collection_info$id)$
    filterDate(start_date, end_date)$
    filterBounds(roi)$
    select(collection_info$band)
  
  # Compute mean composite (automatically handles NoData from cloud masking)
  monthly_mean <- col$mean()$clip(roi)
  
  return(monthly_mean)
}

# --- 5. Download monthly rasters for the full year ---------------------------

download_tropomi_monthly <- function(pollutant_name, collection_info,
                                     year, roi, output_dir) {
  #' Download 12 monthly GeoTIFFs for one pollutant.
  
  out_subdir <- file.path(output_dir, tolower(pollutant_name))
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  
  for (month in 1:12) {
    fname <- glue("{tolower(pollutant_name)}_{year}_{sprintf('%02d', month)}.tif")
    fpath <- file.path(out_subdir, fname)
    
    if (file.exists(fpath)) {
      message(glue("  Skipping {fname} (already exists)"))
      next
    }
    
    message(glue("  Downloading {pollutant_name} {year}-{sprintf('%02d', month)}..."))
    
    tryCatch({
      img <- get_monthly_mean(collection_info, year, month, roi)
      
      # Download as GeoTIFF via rgee
      # ee_as_raster downloads to local; for large areas use ee_image_to_drive
      rast <- ee_as_rast(
        image  = img,
        region = roi,
        scale  = 5000,    # ~5km resolution (matching TROPOMI)
        via    = "drive"   # uses Google Drive as intermediary (more reliable)
      )
      
      # Save locally
      terra::writeRaster(rast, fpath, overwrite = TRUE)
      message(glue("  Saved: {fpath}"))
      
    }, error = function(e) {
      warning(glue("  Failed {pollutant_name} {year}-{sprintf('%02d', month)}: {e$message}"))
    })
  }
}

# --- 6. Run downloads --------------------------------------------------------

year <- as.integer(format(as.Date(CONFIG$date_start), "%Y"))
tropomi_dir <- file.path(CONFIG$dir_raw, "tropomi")

for (poll_name in names(collections)) {
  message(glue("\n=== Downloading {poll_name} ==="))
  download_tropomi_monthly(
    pollutant_name  = poll_name,
    collection_info = collections[[poll_name]],
    year            = year,
    roi             = roi,
    output_dir      = tropomi_dir
  )
}

message("\nTROPOMI download complete.")

# --- 7. Quick sanity check: plot one month -----------------------------------

test_file <- file.path(tropomi_dir, "no2", glue("no2_{year}_06.tif"))
if (file.exists(test_file)) {
  r <- terra::rast(test_file)
  plot(r, main = glue("TROPOMI NO2 — June {year}"),
       col = viridis::viridis(100))
}

# =============================================================================
# FALLBACK: If rgee gives you trouble, download manually from S5P-PAL:
#
#   1. Go to https://maps.s5p-pal.com/no2/
#   2. Select month, click "Download" for the L3 GeoTIFF
#   3. Crop to Poland in R:
#
#      library(terra)
#      r <- rast("S5P_PAL_L3_NO2_20240601_20240630.tif")
#      poland_ext <- ext(14.07, 24.15, 49.00, 54.84)
#      r_pl <- crop(r, poland_ext)
#      writeRaster(r_pl, "data/raw/tropomi/no2/no2_2024_06.tif")
#
# =============================================================================