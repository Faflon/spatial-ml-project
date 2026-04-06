# =============================================================================
# 02_build_grid.R — Create Analysis Grid & Integrate All Data Sources
# =============================================================================
# 1. Creates 5km grid over Poland (EPSG:2180)
# 2. Extracts TROPOMI raster values to grid cells
# 3. Assigns GIOŚ stations to grid cells
# 4. Computes OSM spatial features per cell (road density, industrial %)
# 5. Joins GUS demographics via gmina boundaries
# 6. Outputs a single analysis-ready GeoPackage
# =============================================================================

source("00_setup.R")

library(sf)
library(terra)
library(tidyverse)
library(glue)

# =============================================================================
# 1. BUILD GRID
# =============================================================================

message("\n========== 1. Building grid ==========\n")

# Create Poland bbox in projected CRS
bbox_wgs <- st_bbox(
  base::c(xmin = unname(CONFIG$bbox["xmin"]),
          ymin = unname(CONFIG$bbox["ymin"]),
          xmax = unname(CONFIG$bbox["xmax"]),
          ymax = unname(CONFIG$bbox["ymax"])),
  crs = st_crs(CONFIG$crs_wgs84)
)

bbox_pl <- bbox_wgs |>
  st_as_sfc() |>
  st_transform(CONFIG$crs_pl) |>
  st_bbox()

# Make grid
grid <- st_make_grid(
  st_as_sfc(bbox_pl),
  cellsize = base::c(CONFIG$grid_res, CONFIG$grid_res),
  what = "polygons"
) |>
  st_as_sf() |>
  mutate(cell_id = row_number())

# Rename geometry column if needed
if (!"geometry" %in% names(grid)) {
  grid <- grid |> rename(geometry = x)
}

message(glue("Grid: {nrow(grid)} cells at {CONFIG$grid_res/1000} km"))

# Optional: clip to actual Poland boundary
gminy_path <- file.path(CONFIG$dir_raw, "gus", "gminy.gpkg")
if (file.exists(gminy_path)) {
  message("Clipping grid to Poland boundary...")
  poland <- st_read(gminy_path, quiet = TRUE) |>
    st_transform(CONFIG$crs_pl) |>
    st_union()

  grid <- grid[st_intersects(grid, poland, sparse = FALSE)[, 1], ]
  message(glue("After clip: {nrow(grid)} cells"))
}

# Add centroid coordinates (useful as features in ML)
centroids <- st_coordinates(st_centroid(grid))
grid$x_coord <- centroids[, 1]
grid$y_coord <- centroids[, 2]


# =============================================================================
# 2. EXTRACT TROPOMI VALUES
# =============================================================================

message("\n========== 2. TROPOMI extraction ==========\n")

tropomi_dir <- file.path(CONFIG$dir_raw, "tropomi")
tropomi_all <- base::list()

for (poll in base::c("no2", "so2", "co")) {
  tif_files <- list.files(file.path(tropomi_dir, poll),
                          pattern = "\\.tif$", full.names = TRUE)

  if (length(tif_files) == 0) {
    message(glue("  No TIFFs for {poll} — skipping"))
    next
  }

  message(glue("  Extracting {poll}: {length(tif_files)} files"))

  for (f in tif_files) {
    # Parse date from filename: e.g. no2_2024_06.tif
    parts <- regmatches(basename(f),
                        regexec("(\\w+)_(\\d{4})_(\\d{2})\\.tif", basename(f)))[[1]]

    if (length(parts) < 4) next

    year_month <- as.Date(paste0(parts[3], "-", parts[4], "-01"))

    tryCatch({
      r <- rast(f)

      # Reproject if needed
      r_crs <- crs(r, describe = TRUE)$code
      if (!is.null(r_crs) && r_crs != as.character(CONFIG$crs_pl)) {
        r <- project(r, paste0("EPSG:", CONFIG$crs_pl))
      }

      # Extract mean value per grid cell
      vals <- extract(r, vect(grid), fun = mean, na.rm = TRUE)

      tropomi_all[[basename(f)]] <- tibble(
        cell_id    = grid$cell_id,
        year_month = year_month,
        pollutant  = toupper(poll),
        value      = vals[[2]]
      )
    }, error = function(e) {
      warning(glue("  Failed {basename(f)}: {e$message}"))
    })
  }
}

if (length(tropomi_all) > 0) {
  tropomi_df <- bind_rows(tropomi_all)
  write_csv(tropomi_df, file.path(CONFIG$dir_processed, "tropomi_grid.csv"))
  message(glue("  Saved: {nrow(tropomi_df)} cell-month-pollutant observations"))

  # Add annual NO2 mean to grid
  no2_annual <- tropomi_df |>
    filter(pollutant == "NO2") |>
    group_by(cell_id) |>
    summarise(tropomi_no2 = mean(value, na.rm = TRUE), .groups = "drop")

  grid <- left_join(grid, no2_annual, by = "cell_id")

  # Same for SO2 and CO
  for (p in base::c("SO2", "CO")) {
    annual <- tropomi_df |>
      filter(pollutant == p) |>
      group_by(cell_id) |>
      summarise(val = mean(value, na.rm = TRUE), .groups = "drop")

    col_name <- paste0("tropomi_", tolower(p))
    names(annual)[2] <- col_name
    grid <- left_join(grid, annual, by = "cell_id")
  }
} else {
  message("  No TROPOMI data available yet")
}


# =============================================================================
# 3. GIOŚ STATIONS
# =============================================================================

message("\n========== 3. GIOŚ station assignment ==========\n")

stations_path <- file.path(CONFIG$dir_raw, "gios", "stations.gpkg")

if (file.exists(stations_path)) {
  stations_sf <- st_read(stations_path, quiet = TRUE) |>
    st_transform(CONFIG$crs_pl)

  # Find which grid cell each station belongs to
  station_grid <- st_join(stations_sf, grid |> select(cell_id), join = st_within)

  n_matched <- sum(!is.na(station_grid$cell_id))
  message(glue("  Matched {n_matched}/{nrow(stations_sf)} stations to grid cells"))

  st_write(station_grid,
           file.path(CONFIG$dir_processed, "stations_gridded.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
} else {
  message("  No station data yet (run collect_gios() first)")
}


# =============================================================================
# 4. OSM FEATURES
# =============================================================================

message("\n========== 4. OSM feature engineering ==========\n")

osm_dir <- file.path(CONFIG$dir_raw, "osm")

# --- 4a. Road density (km of major roads per cell) ---
roads_path <- file.path(osm_dir, "major_roads.gpkg")

if (file.exists(roads_path)) {
  message("  Computing road density...")
  roads <- st_read(roads_path, quiet = TRUE)

  # Intersect roads with grid, compute length per cell
  tryCatch({
    road_in_grid <- st_intersection(roads, grid |> select(cell_id))

    road_stats <- road_in_grid |>
      mutate(length_m = as.numeric(st_length(geometry))) |>
      st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(
        total_road_km = sum(length_m) / 1000,
        motorway_km   = sum(length_m[road_class == "motorway"]) / 1000,
        .groups = "drop"
      )

    grid <- left_join(grid, road_stats, by = "cell_id")
    grid$total_road_km[is.na(grid$total_road_km)] <- 0
    grid$motorway_km[is.na(grid$motorway_km)] <- 0

    message(glue("  Road density: done ({nrow(road_stats)} cells with roads)"))
  }, error = function(e) {
    warning(glue("  Road density failed: {e$message}"))
  })
} else {
  message("  No road data (run collect_osm() first)")
}

# --- 4b. Industrial area coverage (% of cell) ---
ind_path <- file.path(osm_dir, "industrial.gpkg")

if (file.exists(ind_path)) {
  message("  Computing industrial coverage...")
  tryCatch({
    industrial <- st_read(ind_path, quiet = TRUE)
    ind_in_grid <- st_intersection(industrial, grid |> select(cell_id))

    ind_stats <- ind_in_grid |>
      mutate(area_m2 = as.numeric(st_area(geometry))) |>
      st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(industrial_m2 = sum(area_m2), .groups = "drop")

    grid <- left_join(grid, ind_stats, by = "cell_id")
    cell_area <- as.numeric(st_area(grid[1, ]))
    grid$industrial_pct <- ifelse(
      is.na(grid$industrial_m2), 0,
      grid$industrial_m2 / cell_area * 100
    )
    grid$industrial_m2 <- NULL

    message("  Industrial coverage: done")
  }, error = function(e) {
    warning(glue("  Industrial coverage failed: {e$message}"))
  })
} else {
  message("  No industrial data")
}

# --- 4c. Power plant count (within 10km of cell centroid) ---
pp_path <- file.path(osm_dir, "power_plants.gpkg")

if (file.exists(pp_path)) {
  message("  Counting nearby power plants...")
  tryCatch({
    pp <- st_read(pp_path, quiet = TRUE)
    buffers <- st_buffer(st_centroid(grid), 10000)
    grid$power_plants_10km <- lengths(st_intersects(buffers, pp))
    message("  Power plant count: done")
  }, error = function(e) {
    warning(glue("  Power plant count failed: {e$message}"))
  })
} else {
  message("  No power plant data")
}

# --- 4d. Railway density ---
rail_path <- file.path(osm_dir, "railways.gpkg")

if (file.exists(rail_path)) {
  message("  Computing railway density...")
  tryCatch({
    rail <- st_read(rail_path, quiet = TRUE)
    rail_in_grid <- st_intersection(rail, grid |> select(cell_id))

    rail_stats <- rail_in_grid |>
      mutate(length_m = as.numeric(st_length(geometry))) |>
      st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(railway_km = sum(length_m) / 1000, .groups = "drop")

    grid <- left_join(grid, rail_stats, by = "cell_id")
    grid$railway_km[is.na(grid$railway_km)] <- 0
    message("  Railway density: done")
  }, error = function(e) {
    warning(glue("  Railway density failed: {e$message}"))
  })
} else {
  message("  No railway data")
}


# =============================================================================
# 5. GUS DEMOGRAPHICS (if available)
# =============================================================================

message("\n========== 5. GUS demographics ==========\n")

gus_path <- file.path(CONFIG$dir_raw, "gus", "bdl_data.csv")

if (file.exists(gus_path) && file.exists(gminy_path)) {
  message("  Joining GUS data via gmina boundaries...")
  # Spatial join: each grid centroid gets the gmina it falls in
  # Then join BDL data by gmina ID (TERYT)
  # (Implementation depends on actual BDL column structure)
  message("  TODO: implement after BDL variable IDs are selected")
} else {
  message("  GUS data not ready yet")
}


# =============================================================================
# 6. SAVE FINAL GRID
# =============================================================================

message("\n========== 6. Saving ==========\n")

st_write(grid,
         file.path(CONFIG$dir_processed, "analysis_grid.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

# Print summary
feature_cols <- setdiff(names(grid),
                        base::c("geometry", "cell_id", "x", "x_coord", "y_coord"))
message(glue("Grid saved: {nrow(grid)} cells"))
message(glue("Features ({length(feature_cols)}): {paste(feature_cols, collapse = ', ')}"))
message(glue("Output: {CONFIG$dir_processed}/analysis_grid.gpkg"))

# Quick summary stats
for (col in feature_cols) {
  if (is.numeric(grid[[col]])) {
    vals <- grid[[col]][!is.na(grid[[col]])]
    if (length(vals) > 0) {
      message(sprintf("  %-20s  mean=%.4g  sd=%.4g  n=%d",
                      col, mean(vals), sd(vals), length(vals)))
    }
  }
}

message("\n02_build_grid.R complete.")

# =============================================================================
