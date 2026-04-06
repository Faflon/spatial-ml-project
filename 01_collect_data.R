# =============================================================================
# 01_collect_data.R — Data Collection (All Sources)
# =============================================================================
# Run sections independently or all at once.
# Each section saves outputs to data/raw/ and prints a status summary.
#
# Sections:
#   A. Sentinel-5P TROPOMI (via Google Earth Engine)
#   B. GIOŚ ground stations (via REST API + bulk download)
#   C. GUS BDL demographics (via REST API)
#   D. OpenStreetMap features (via osmdata / Overpass API)
# =============================================================================

source("00_setup.R")

library(sf)
library(tidyverse)
library(httr2)
library(jsonlite)
library(glue)

# --- Progress tracker --------------------------------------------------------

progress_tracker <- function(total, label = "Progress") {
  #' Creates a progress tracker that shows elapsed time, ETA, and % done.
  #' Usage:
  #'   tick <- progress_tracker(100, "Downloading")
  #'   for (i in 1:100) { do_work(); tick(i) }
  
  start_time <- proc.time()["elapsed"]
  
  function(i, extra_info = "") {
    elapsed <- proc.time()["elapsed"] - start_time
    pct <- i / total * 100
    per_item <- elapsed / i
    eta <- per_item * (total - i)
    
    # Format times
    fmt_time <- function(seconds) {
      seconds <- round(seconds)
      if (seconds < 60) return(paste0(seconds, "s"))
      if (seconds < 3600) return(sprintf("%dm %02ds", seconds %/% 60, seconds %% 60))
      sprintf("%dh %02dm", seconds %/% 3600, (seconds %% 3600) %/% 60)
    }
    
    msg <- sprintf("\r  [%3.0f%%] %d/%d | Elapsed: %s | ETA: %s",
                   pct, i, total, fmt_time(elapsed), fmt_time(eta))
    if (nchar(extra_info) > 0) msg <- paste0(msg, " | ", extra_info)
    
    cat(msg)
    if (i == total) cat("\n")
    flush.console()
  }
}


# #############################################################################
# A. SENTINEL-5P TROPOMI                                                     #
# #############################################################################

collect_tropomi <- function() {
  
  library(rgee)
  library(terra)
  
  message("\n========== A. TROPOMI (Sentinel-5P) ==========\n")
  
  # --- A1. Initialize GEE ---
  ee_Initialize()
  
  # --- A2. Region of interest (Poland bbox) ---
  roi <- ee$Geometry$Rectangle(
    coords = as.list(as.numeric(CONFIG$bbox)),
    proj   = "EPSG:4326",
    geodesic = FALSE
  )
  
  # --- A3. Define products ---
  products <- base::list(
    no2 = base::list(
      collection = "COPERNICUS/S5P/OFFL/L3_NO2",
      band       = "tropospheric_NO2_column_number_density"
    ),
    so2 = base::list(
      collection = "COPERNICUS/S5P/OFFL/L3_SO2",
      band       = "SO2_column_number_density"
    ),
    co = base::list(
      collection = "COPERNICUS/S5P/OFFL/L3_CO",
      band       = "CO_column_number_density"
    )
  )
  
  year <- as.integer(format(as.Date(CONFIG$date_start), "%Y"))
  
  # --- A4. Download monthly composites ---
  total_downloads <- length(products) * 12
  current <- 0
  tick <- progress_tracker(total_downloads, "TROPOMI")
  
  for (poll_name in names(products)) {
    poll <- products[[poll_name]]
    out_dir <- file.path(CONFIG$dir_raw, "tropomi", poll_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (month in 1:12) {
      current <- current + 1
      fname <- sprintf("%s_%d_%02d.tif", poll_name, year, month)
      fpath <- file.path(out_dir, fname)
      
      if (file.exists(fpath)) {
        tick(current, glue("Skip: {fname}"))
        next
      }
      
      start_d <- sprintf("%04d-%02d-01", year, month)
      end_d   <- if (month == 12) {
        sprintf("%04d-01-01", year + 1)
      } else {
        sprintf("%04d-%02d-01", year, month + 1)
      }
      
      message(glue("\n  Downloading {poll_name} {year}-{sprintf('%02d', month)}..."))
      
      tryCatch({
        img <- ee$ImageCollection(poll$collection)$
          filterDate(start_d, end_d)$
          filterBounds(roi)$
          select(poll$band)$
          mean()$
          clip(roi)
        
        r <- ee_as_rast(
          image  = img,
          region = roi,
          scale  = 5000,
          via    = "drive"
        )
        
        terra::writeRaster(r, fpath, overwrite = TRUE)
        tick(current, glue("OK: {fname}"))
        
      }, error = function(e) {
        tick(current, glue("FAIL: {fname}"))
        warning(glue("\n  FAILED {fname}: {e$message}"))
      })
    }
  }
  
  # --- A5. Status ---
  tifs <- list.files(file.path(CONFIG$dir_raw, "tropomi"),
                     pattern = "\\.tif$", recursive = TRUE)
  message(glue("\nTROPOMI done: {length(tifs)} GeoTIFFs saved"))
  message("Location: data/raw/tropomi/{no2,so2,co}/")
}

# FALLBACK if rgee doesn't work:
# 1. Go to https://maps.s5p-pal.com/no2/
# 2. Download monthly L3 GeoTIFFs
# 3. In R:
#    library(terra)
#    r <- rast("downloaded_file.tif")
#    r_pl <- crop(r, ext(14.07, 24.15, 49.00, 54.84))
#    writeRaster(r_pl, "data/raw/tropomi/no2/no2_2024_01.tif")


# #############################################################################
# B. GIOŚ GROUND STATIONS                                                    #
# #############################################################################

collect_gios <- function() {
  
  message("\n========== B. GIOŚ Ground Stations ==========\n")
  
  gios_start <- proc.time()["elapsed"]
  BASE_URL <- "https://api.gios.gov.pl/pjp-api/v1/rest"
  
  # Helper: GET with retry
  gios_get <- function(endpoint, delay = 0.5) {
    url <- paste0(BASE_URL, endpoint)
    resp <- tryCatch({
      request(url) |>
        req_headers("Accept" = "application/json") |>
        req_retry(max_tries = 3, backoff = ~ 2) |>
        req_perform()
    }, error = function(e) {
      warning(paste0("GIOŚ API error: ", e$message))
      return(NULL)
    })
    Sys.sleep(delay)
    if (is.null(resp)) return(NULL)
    resp_body_json(resp, simplifyVector = TRUE)
  }
  
  # --- B1. Station list ---
  message("Fetching stations...")
  raw <- gios_get("/station/findAll", delay = 1)
  
  if (is.null(raw) || length(raw) == 0) {
    stop("Failed to fetch station list. Check network / API status.")
  }
  
  stations <- tibble(
    station_id = raw$id,
    name       = raw$stationName,
    lat        = as.numeric(raw$gegrLat),
    lon        = as.numeric(raw$gegrLon),
    city       = raw$city$name,
    commune    = raw$city$commune$communeName,
    province   = raw$city$commune$provinceName
  ) |>
    filter(!is.na(lat), !is.na(lon))
  
  message(glue("  Found {nrow(stations)} stations"))
  
  # Save as CSV + GeoPackage
  write_csv(stations, file.path(CONFIG$dir_raw, "gios", "stations.csv"))
  
  stations_sf <- st_as_sf(stations, coords = base::c("lon", "lat"), crs = 4326L)
  st_write(stations_sf,
           file.path(CONFIG$dir_raw, "gios", "stations.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  
  # --- B2. Sensors per station ---
  # NOTE: metadata endpoint has 2 req/min limit — this is SLOW (~2h for all)
  # For a quick start, skip this and use bulk CSV download instead.
  
  message("\nFetching sensors (2 req/min limit — this takes a while)...")
  message("TIP: For faster start, download bulk data from:")
  message("  https://powietrze.gios.gov.pl/pjp/archives")
  
  all_sensors <- base::list()
  tick_sensors <- progress_tracker(nrow(stations), "Sensors")
  
  for (i in seq_len(nrow(stations))) {
    sid <- stations$station_id[i]
    
    sensors <- gios_get(glue("/station/sensors/{sid}"), delay = 31)
    
    if (!is.null(sensors) && length(sensors) > 0) {
      df <- tibble(
        sensor_id    = sensors$id,
        station_id   = sensors$stationId,
        param_name   = sensors$param$paramName,
        param_code   = sensors$param$paramFormula
      )
      all_sensors[[i]] <- df
    }
    
    tick_sensors(i, stations$name[i])
  }
  
  sensors_df <- bind_rows(all_sensors)
  write_csv(sensors_df, file.path(CONFIG$dir_raw, "gios", "sensors.csv"))
  
  target_params <- base::c("NO2", "SO2", "CO", "O3", "PM10", "PM2.5")
  sensors_target <- sensors_df |> filter(param_code %in% target_params)
  
  message(glue("\nSensors total: {nrow(sensors_df)}"))
  message(glue("Target sensors (NO2/SO2/CO/O3/PM): {nrow(sensors_target)}"))
  
  # --- B3. Download measurement data (recent, last 72h) ---
  message("\nDownloading recent measurements...")
  
  recent <- base::list()
  tick_meas <- progress_tracker(nrow(sensors_target), "Measurements")
  
  for (i in seq_len(nrow(sensors_target))) {
    sid <- sensors_target$sensor_id[i]
    d <- gios_get(glue("/data/getData/{sid}"), delay = 0.05)
    
    if (!is.null(d) && !is.null(d$values)) {
      vals <- tibble(
        sensor_id  = sid,
        station_id = sensors_target$station_id[i],
        param_code = sensors_target$param_code[i],
        datetime   = as.POSIXct(d$values$date, format = "%Y-%m-%d %H:%M:%S"),
        value      = as.numeric(d$values$value)
      ) |>
        filter(!is.na(value))
      
      if (nrow(vals) > 0) recent[[i]] <- vals
    }
    
    tick_meas(i, sensors_target$param_code[i])
  }
  
  recent_df <- bind_rows(recent)
  write_csv(recent_df, file.path(CONFIG$dir_raw, "gios", "recent_measurements.csv"))
  
  message(glue("\nRecent data: {nrow(recent_df)} observations"))
  
  # --- B4. Status ---
  gios_elapsed <- round(proc.time()["elapsed"] - gios_start)
  message(glue("\nGIOŚ done. Total time: {gios_elapsed %/% 60}m {gios_elapsed %% 60}s"))
  message("For HISTORICAL data (full year), download bulk CSV from:")
  message("  https://powietrze.gios.gov.pl/pjp/archives")
  message("  → 'Wyniki pomiarów [year]' → save ZIP to data/raw/gios/")
}

# Helper: parse bulk ZIP file (call after manual download)
parse_gios_bulk <- function(zip_path) {
  #' Parse a GIOŚ bulk measurement ZIP into a tidy data.frame.
  #' Call: parse_gios_bulk("data/raw/gios/2024_Wyniki_pomiarow.zip")
  
  tmp <- tempdir()
  files <- unzip(zip_path, exdir = tmp)
  message(glue("Extracted {length(files)} files from ZIP"))
  
  all_data <- base::list()
  for (f in files) {
    tryCatch({
      df <- read_csv(f, show_col_types = FALSE,
                     locale = locale(encoding = "Windows-1250"))
      all_data[[basename(f)]] <- df
      message(glue("  Parsed: {basename(f)} ({nrow(df)} rows)"))
    }, error = function(e) {
      warning(glue("  Failed: {basename(f)}: {e$message}"))
    })
  }
  
  result <- bind_rows(all_data, .id = "source_file")
  message(glue("\nTotal: {nrow(result)} rows"))
  return(result)
}


# #############################################################################
# C. GUS BDL DEMOGRAPHICS                                                    #
# #############################################################################

collect_gus <- function() {
  
  message("\n========== C. GUS BDL Demographics ==========\n")
  
  BDL_BASE <- "https://bdl.stat.gov.pl/api/v1"
  
  bdl_get <- function(endpoint, params = base::list(), delay = 1.5) {
    url <- paste0(BDL_BASE, endpoint)
    req <- request(url) |>
      req_url_query(!!!params, format = "json") |>
      req_headers("Accept" = "application/json")
    
    resp <- tryCatch({
      req |> req_retry(max_tries = 3, backoff = ~ 3) |> req_perform()
    }, error = function(e) {
      warning(paste0("BDL error: ", e$message))
      return(NULL)
    })
    
    Sys.sleep(delay)
    if (is.null(resp)) return(NULL)
    resp_body_json(resp, simplifyVector = TRUE)
  }
  
  # --- C1. Search for relevant variables ---
  message("Searching for relevant BDL variables...")
  
  search_terms <- base::c(
    "ludność",
    "gęstość zaludnienia",
    "pojazdy samochodowe",
    "powierzchnia"
  )
  
  search_results <- base::list()
  for (term in search_terms) {
    message(glue("  Searching: '{term}'"))
    result <- bdl_get("/variables/search", params = base::list(name = term))
    if (!is.null(result) && !is.null(result$results) && length(result$results) > 0) {
      df <- as_tibble(result$results)
      df$search_term <- term
      search_results[[term]] <- df
    }
    Sys.sleep(2)  # be gentle with API
  }
  
  all_vars <- bind_rows(search_results)
  
  if (nrow(all_vars) > 0) {
    # Keep useful columns (varies by API response, so select safely)
    cols_keep <- intersect(
      names(all_vars),
      base::c("id", "n1", "n2", "n3", "subjectId", "measureUnitName",
              "level", "search_term")
    )
    all_vars <- all_vars[, cols_keep]
    write_csv(all_vars, file.path(CONFIG$dir_raw, "gus", "bdl_variable_search.csv"))
    message(glue("\nFound {nrow(all_vars)} variables. Saved to data/raw/gus/bdl_variable_search.csv"))
    message("Review this file and pick variable IDs for download.")
  } else {
    message("No variables found — check API connectivity.")
  }
  
  # --- C2. Download data for a specific variable ---
  # After reviewing the search results, fill in actual variable IDs here
  # and run this section.
  
  download_bdl_variable <- function(var_id, var_name, year = 2024) {
    message(glue("  Downloading '{var_name}' (id={var_id}) for {year}..."))
    
    all_data <- base::list()
    page <- 0
    page_size <- 100
    
    repeat {
      params <- base::list(year = year, `unit-level` = 6,
                           `page-size` = page_size, page = page)
      result <- bdl_get(glue("/data/by-variable/{var_id}"), params)
      
      if (is.null(result) || is.null(result$results) || length(result$results) == 0) break
      
      all_data <- base::c(all_data, base::list(as_tibble(result$results)))
      page <- page + 1
      
      total <- result$totalRecords %||% 0
      fetched <- page * page_size
      if (fetched >= total) break
    }
    
    if (length(all_data) == 0) {
      message(glue("    No data for {var_name}"))
      return(NULL)
    }
    
    df <- bind_rows(all_data)
    df$variable_name <- var_name
    message(glue("    Got {nrow(df)} rows"))
    return(df)
  }
  
  # EXAMPLE (uncomment and fill with real IDs after reviewing search results):
  # pop_data <- download_bdl_variable("72305", "population_total", 2024)
  # density_data <- download_bdl_variable("60559", "population_density", 2024)
  # gus_all <- bind_rows(pop_data, density_data)
  # write_csv(gus_all, file.path(CONFIG$dir_raw, "gus", "bdl_data.csv"))
  
  # --- C3. Gmina boundaries ---
  message("\nFor gmina boundaries, download from GADM:")
  message("  library(geodata)")
  message('  pol <- gadm(country = "POL", level = 2, path = "data/raw/gus")')
  message("  pol_sf <- st_as_sf(pol)")
  message('  st_write(pol_sf, "data/raw/gus/gminy.gpkg")')
  
  message("\nGUS done.")
}


# #############################################################################
# D. OPENSTREETMAP FEATURES                                                  #
# #############################################################################

collect_osm <- function() {
  
  library(osmdata)
  
  message("\n========== D. OpenStreetMap Features ==========\n")
  
  osm_dir <- file.path(CONFIG$dir_raw, "osm")
  osm_start <- proc.time()["elapsed"]
  
  # Helper: query OSM with error handling and timing
  osm_get <- function(key, value, geom_type = "osm_lines", timeout = 180) {
    t0 <- proc.time()["elapsed"]
    message(glue("  Querying: {key}={value}..."), appendLF = FALSE)
    
    result <- tryCatch({
      opq(bbox = as.numeric(CONFIG$bbox),
          timeout = timeout) |>
        add_osm_feature(key = key, value = value) |>
        osmdata_sf()
    }, error = function(e) {
      warning(glue("  OSM query failed: {e$message}"))
      return(NULL)
    })
    
    elapsed <- round(proc.time()["elapsed"] - t0)
    
    if (is.null(result)) { message(glue(" FAILED ({elapsed}s)")); return(NULL) }
    geom <- result[[geom_type]]
    if (is.null(geom) || nrow(geom) == 0) {
      message(glue(" empty ({elapsed}s)"))
      return(NULL)
    }
    message(glue(" {nrow(geom)} features ({elapsed}s)"))
    return(geom)
  }
  
  # --- D1. Major roads ---
  message("[1/4] Roads:")
  road_types <- base::c("motorway", "trunk", "primary", "secondary")
  road_list <- base::list()
  
  for (rt in road_types) {
    r <- osm_get("highway", rt, "osm_lines")
    if (!is.null(r)) {
      road_list[[rt]] <- r |>
        transmute(osm_id = osm_id, road_class = rt) |>
        st_transform(CONFIG$crs_pl)
    }
  }
  
  if (length(road_list) > 0) {
    roads <- bind_rows(road_list)
    st_write(roads, file.path(osm_dir, "major_roads.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message(glue("  Saved {nrow(roads)} road segments"))
  }
  
  # --- D2. Railways ---
  message("\n[2/4] Railways:")
  rail <- osm_get("railway", "rail", "osm_lines")
  if (!is.null(rail)) {
    rail <- rail |>
      transmute(osm_id = osm_id) |>
      st_transform(CONFIG$crs_pl)
    st_write(rail, file.path(osm_dir, "railways.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message(glue("  Saved {nrow(rail)} segments"))
  }
  
  # --- D3. Industrial areas ---
  message("\n[3/4] Industrial areas:")
  ind <- osm_get("landuse", "industrial", "osm_polygons")
  if (!is.null(ind)) {
    ind <- ind |>
      transmute(osm_id = osm_id, name = name) |>
      st_transform(CONFIG$crs_pl)
    st_write(ind, file.path(osm_dir, "industrial.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message(glue("  Saved {nrow(ind)} polygons"))
  }
  
  # --- D4. Power plants ---
  message("\n[4/4] Power plants:")
  pp <- osm_get("power", "plant", "osm_polygons")
  if (!is.null(pp)) {
    pp <- pp |>
      transmute(osm_id = osm_id, name = name) |>
      st_transform(CONFIG$crs_pl) |>
      st_centroid()
    st_write(pp, file.path(osm_dir, "power_plants.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message(glue("  Saved {nrow(pp)} plants"))
  }
  
  # --- D5. Status ---
  osm_elapsed <- round(proc.time()["elapsed"] - osm_start)
  gpkg_files <- list.files(osm_dir, pattern = "\\.gpkg$")
  message(glue("\nOSM done: {length(gpkg_files)} layers saved to data/raw/osm/"))
  message(glue("Total time: {osm_elapsed %/% 60}m {osm_elapsed %% 60}s"))
}


# #############################################################################
# RUN SECTIONS                                                                #
# #############################################################################
# Uncomment the sections you want to run:

# collect_tropomi()   # ~30 min (needs rgee + GEE auth)
# collect_gios()      # ~2h (slow metadata API) or quick with bulk download
# collect_gus()       # ~5 min (search), then manual variable selection
collect_osm()       # ~10 min (Overpass API queries)

message("\n=== 01_collect_data.R ready ===")
message("Uncomment collect_*() calls above to run each section.")
message("Recommended order: collect_osm() → collect_gios() → collect_gus() → collect_tropomi()")

# =============================================================================