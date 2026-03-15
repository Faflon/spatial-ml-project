# =============================================================================
# 02_collect_gios.R — GIOŚ Ground Station Air Quality Data
# =============================================================================
# Data source: GIOŚ "Jakość Powietrza" API (v1)
# Docs:        https://powietrze.gios.gov.pl/pjp/content/api
# OpenAPI:     https://dev.api.gios.gov.pl/pjp-api/swagger-ui/#/
# What:        Hourly NO2, SO2, PM10, PM2.5, CO, O3 from ~250 stations
# Output:      Station metadata (sf) + measurement time series (data.frame)
# =============================================================================

source("00_setup.R")

library(httr2)
library(jsonlite)
library(sf)
library(tidyverse)
library(glue)
library(lubridate)

# --- API Configuration -------------------------------------------------------

BASE_URL <- "https://api.gios.gov.pl/pjp-api/v1/rest"

# Target pollutants (matching what TROPOMI measures + PM for health context)
TARGET_PARAMS <- c("NO2", "SO2", "CO", "O3", "PM10", "PM2.5")

# Rate limits: metadata endpoints = 2 req/min, data endpoints = 1500 req/min
# We add small delays to be safe.

# --- Helper: GET with retry --------------------------------------------------

gios_get <- function(endpoint, delay = 0.1) {
  url <- paste0(BASE_URL, endpoint)
  
  resp <- tryCatch({
    request(url) |>
      req_headers("Accept" = "application/json") |>
      req_retry(max_tries = 3, backoff = ~2) |>
      req_perform()
  }, error = function(e) {
    warning(glue("API error for {endpoint}: {e$message}"))
    return(NULL)
  })
  
  if (is.null(resp)) return(NULL)
  
  Sys.sleep(delay)
  resp_body_json(resp, simplifyVector = TRUE)
}

# --- 1. Get all stations -----------------------------------------------------

message("Fetching station list...")
stations_raw <- gios_get("/station/findAll", delay = 1)

# The response is a list/data.frame with columns:
# id, stationName, gegrLat, gegrLon, city (nested), addressStreet, ...

stations <- stations_raw |>
  as_tibble() |>
  mutate(
    station_id = id,
    name       = stationName,
    lat        = as.numeric(gegrLat),
    lon        = as.numeric(gegrLon),
    city_name  = city$name,
    commune    = city$commune$communeName,
    district   = city$commune$districtName,
    province   = city$commune$provinceName
  ) |>
  select(station_id, name, lat, lon, city_name, commune, district, province) |>
  filter(!is.na(lat), !is.na(lon))

message(glue("Found {nrow(stations)} stations"))

# Convert to sf object
stations_sf <- st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)

# Save station metadata
write_csv(stations, file.path(CONFIG$dir_raw, "gios", "stations_metadata.csv"))
st_write(stations_sf, file.path(CONFIG$dir_raw, "gios", "stations.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
message("Station metadata saved.")

# --- 2. Get sensors (measurement devices) per station ------------------------

message("\nFetching sensor info for all stations...")
Sys.sleep(30)  # respect 2 req/min for metadata endpoints

all_sensors <- list()

for (i in seq_len(nrow(stations))) {
  sid <- stations$station_id[i]
  
  if (i > 1) Sys.sleep(30)  # 2 req/min limit on metadata endpoints
  
  sensors <- gios_get(glue("/station/sensors/{sid}"), delay = 0)
  
  if (is.null(sensors) || length(sensors) == 0) next
  
  sensors_df <- as_tibble(sensors) |>
    mutate(
      sensor_id   = id,
      station_id  = stationId,
      param_name  = param$paramName,
      param_code  = param$paramCode,
      param_formula = param$paramFormula
    ) |>
    select(sensor_id, station_id, param_name, param_code, param_formula)
  
  all_sensors[[i]] <- sensors_df
  
  if (i %% 10 == 0) message(glue("  Processed {i}/{nrow(stations)} stations"))
}

sensors_df <- bind_rows(all_sensors)
message(glue("Found {nrow(sensors_df)} sensors across all stations"))

# Filter to target pollutants
sensors_target <- sensors_df |>
  filter(param_formula %in% TARGET_PARAMS | param_code %in% TARGET_PARAMS)

message(glue("Target sensors (NO2, SO2, CO, O3, PM10, PM2.5): {nrow(sensors_target)}"))

write_csv(sensors_df, file.path(CONFIG$dir_raw, "gios", "sensors_all.csv"))
write_csv(sensors_target, file.path(CONFIG$dir_raw, "gios", "sensors_target.csv"))

# --- 3. Download measurement data -------------------------------------------
# NOTE: The GIOŚ real-time API only provides last ~72h of data.
# For historical data, use the "Bank danych pomiarowych" bulk download:
#   https://powietrze.gios.gov.pl/pjp/archives
#
# The bulk CSV files contain hourly data for full years.
# Below we show BOTH approaches.

# --- 3a. APPROACH A: Bulk download (RECOMMENDED for historical data) ---------

message("\n=== APPROACH A: Bulk Historical Data ===")
message("For full-year historical data, download ZIP files manually from:")
message("  https://powietrze.gios.gov.pl/pjp/archives")
message("  → 'Przygotowane dane do pobrania' → select year")
message("  → Download 'Wyniki pomiarów [year]'")
message("  → Place ZIP in: data/raw/gios/")
message("")
message("The ZIP contains CSV files with hourly measurements from all stations.")
message("Columns: station code, date, hour, value (µg/m³)")

# After manual download, parse with:
parse_gios_bulk <- function(zip_path, pollutant_code = "NO2") {
  #' Parse GIOŚ bulk measurement ZIP file for a specific pollutant.
  #' Returns a tidy data.frame with: station_code, datetime, value
  
  # List CSV files in the ZIP
  csv_files <- unzip(zip_path, list = TRUE)$Name
  target_csv <- csv_files[grepl(pollutant_code, csv_files, ignore.case = TRUE)]
  
  if (length(target_csv) == 0) {
    warning(glue("No CSV found for {pollutant_code} in {zip_path}"))
    return(NULL)
  }
  
  # Extract and read
  tmp_dir <- tempdir()
  unzip(zip_path, files = target_csv, exdir = tmp_dir)
  
  df <- read_csv(file.path(tmp_dir, target_csv[1]),
                 show_col_types = FALSE,
                 locale = locale(encoding = "Windows-1250"))
  
  message(glue("  Parsed {pollutant_code}: {nrow(df)} rows"))
  return(df)
}

# --- 3b. APPROACH B: API for recent/current data ----------------------------

message("\n=== APPROACH B: API (recent data, last 72h) ===")

download_sensor_data <- function(sensor_id) {
  #' Download recent measurement data for one sensor via API.
  #' Returns data.frame with datetime and value columns.
  
  data <- gios_get(glue("/data/getData/{sensor_id}"), delay = 0.05)
  
  if (is.null(data) || is.null(data$values)) return(NULL)
  
  tibble(
    sensor_id = sensor_id,
    datetime  = as.POSIXct(data$values$date, format = "%Y-%m-%d %H:%M:%S"),
    value     = as.numeric(data$values$value)
  ) |>
    filter(!is.na(value))
}

# Download for all target sensors (this uses the 1500 req/min data endpoint)
message("Downloading recent data for target sensors...")

recent_data <- list()
for (i in seq_len(nrow(sensors_target))) {
  sid <- sensors_target$sensor_id[i]
  d <- download_sensor_data(sid)
  if (!is.null(d) && nrow(d) > 0) {
    d$param_code <- sensors_target$param_code[i]
    d$station_id <- sensors_target$station_id[i]
    recent_data[[i]] <- d
  }
  if (i %% 100 == 0) message(glue("  Downloaded {i}/{nrow(sensors_target)}"))
}

recent_df <- bind_rows(recent_data)
message(glue("Recent data: {nrow(recent_df)} observations"))

write_csv(recent_df, file.path(CONFIG$dir_raw, "gios", "recent_measurements.csv"))

# --- 4. Aggregate to monthly station means -----------------------------------

aggregate_to_monthly <- function(measurements_df) {
  #' Aggregate hourly measurements to monthly means per station & pollutant.
  
  measurements_df |>
    mutate(year_month = floor_date(datetime, "month")) |>
    group_by(station_id, param_code, year_month) |>
    summarise(
      mean_value   = mean(value, na.rm = TRUE),
      median_value = median(value, na.rm = TRUE),
      max_value    = max(value, na.rm = TRUE),
      n_obs        = sum(!is.na(value)),
      coverage_pct = n_obs / (days_in_month(first(year_month)) * 24) * 100,
      .groups = "drop"
    ) |>
    # Keep only months with >= 75% data coverage
    filter(coverage_pct >= 75)
}

# --- 5. Summary --------------------------------------------------------------

message("\n=== GIOŚ Data Collection Summary ===")
message(glue("Stations:       {nrow(stations)}"))
message(glue("Target sensors: {nrow(sensors_target)}"))
message(glue("Files saved in: {file.path(CONFIG$dir_raw, 'gios')}"))
message("")
message("NEXT STEPS:")
message("  1. Download bulk historical ZIP from powietrze.gios.gov.pl/pjp/archives")
message("  2. Run parse_gios_bulk() on the ZIP")
message("  3. Run aggregate_to_monthly() to get station-month means")
message("  4. Join with stations_sf for spatial analysis")

# =============================================================================