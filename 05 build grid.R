# =============================================================================
# 05_build_grid.R — Create Analysis Grid & Integrate All Data Sources
# =============================================================================
# This script:
#   1. Creates a regular 5km grid over Poland (EPSG:2180)
#   2. Extracts TROPOMI raster values to grid cells
#   3. Joins GIOŚ station data to nearest grid cells
#   4. Joins GUS demographic data via spatial overlay
#   5. Computes OSM-derived spatial features per cell
#   6. Outputs a single analysis-ready sf object
# =============================================================================

source("00_setup.R")

library(sf)
library(terra)
library(stars)
library(tidyverse)
library(glue)

# --- 1. Create regular grid over Poland -------------------------------------

message("=== Building 5km analysis grid ===")

# Poland boundary — use GADM or similar
# Option A: simple bbox-based grid
poland_bbox_2180 <- st_bbox(
  st_transform(
    st_as_sfc(st_bbox(c(
      xmin = CONFIG$bbox_poland["xmin"],
      ymin = CONFIG$bbox_poland["ymin"],
      xmax = CONFIG$bbox_poland["xmax"],
      ymax = CONFIG$bbox_poland["ymax"]
    ), crs = 4326)),
    CONFIG$crs_pl
  )
)

# Create grid
grid <- st_make_grid(
  st_as_sfc(poland_bbox_2180),
  cellsize = c(CONFIG$grid_res, CONFIG$grid_res),  # 5km x 5km
  what = "polygons"
) |>
  st_as_sf() |>
  rename(geometry = x) |>
  mutate(cell_id = row_number())

message(glue("Grid created: {nrow(grid)} cells ({CONFIG$grid_res/1000}km resolution)"))

# Option B (better): clip grid to actual Poland boundary
# Uncomment after downloading Poland boundary:
#
# poland_boundary <- st_read("data/raw/gus/gadm41_POL_0.gpkg") |>
#   st_transform(CONFIG$crs_pl)
# grid <- grid[st_intersects(grid, poland_boundary, sparse = FALSE)[,1], ]
# message(glue("After clipping to Poland: {nrow(grid)} cells"))

# --- 2. Extract TROPOMI values to grid cells ---------------------------------

message("\n=== Extracting TROPOMI raster values ===")

extract_tropomi_to_grid <- function(raster_dir, pollutant, grid_sf) {
  #' Read monthly GeoTIFFs and extract mean values per grid cell.
  #' Returns data.frame: cell_id, year_month, pollutant, value
  
  tif_files <- list.files(
    file.path(raster_dir, pollutant),
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  if (length(tif_files) == 0) {
    message(glue("  No TIFF files found for {pollutant}"))
    return(NULL)
  }
  
  results <- list()
  
  for (f in tif_files) {
    # Parse year-month from filename (e.g., no2_2024_06.tif)
    fname <- basename(f)
    parts <- str_match(fname, "(\\w+)_(\\d{4})_(\\d{2})\\.tif")
    ym <- paste0(parts[3], "-", parts[4], "-01")
    
    message(glue("  Extracting {pollutant} {parts[3]}-{parts[4]}..."))
    
    r <- terra::rast(f)
    
    # Reproject raster to match grid CRS if needed
    if (terra::crs(r, describe = TRUE)$code != as.character(CONFIG$crs_pl)) {
      r <- terra::project(r, paste0("EPSG:", CONFIG$crs_pl))
    }
    
    # Extract mean value per grid cell
    vals <- terra::extract(r, terra::vect(grid_sf), fun = mean, na.rm = TRUE)
    
    results[[fname]] <- tibble(
      cell_id    = grid_sf$cell_id,
      year_month = as.Date(ym),
      pollutant  = toupper(pollutant),
      value      = vals[[2]]  # first column is ID, second is value
    )
  }
  
  bind_rows(results)
}

tropomi_dir <- file.path(CONFIG$dir_raw, "tropomi")

tropomi_grid <- bind_rows(
  extract_tropomi_to_grid(tropomi_dir, "no2", grid),
  extract_tropomi_to_grid(tropomi_dir, "so2", grid),
  extract_tropomi_to_grid(tropomi_dir, "co", grid)
)

if (!is.null(tropomi_grid) && nrow(tropomi_grid) > 0) {
  message(glue("TROPOMI extraction: {nrow(tropomi_grid)} cell-month-pollutant obs"))
} else {
  message("No TROPOMI data extracted yet (run 01_collect_tropomi.R first)")
}

# --- 3. Assign GIOŚ stations to grid cells ----------------------------------

message("\n=== Assigning GIOŚ stations to grid ===")

stations_path <- file.path(CONFIG$dir_raw, "gios", "stations.gpkg")

if (file.exists(stations_path)) {
  stations_sf <- st_read(stations_path, quiet = TRUE) |>
    st_transform(CONFIG$crs_pl)
  
  # Find which grid cell each station falls in
  station_grid <- st_join(stations_sf, grid, join = st_within)
  message(glue("Matched {sum(!is.na(station_grid$cell_id))} stations to grid cells"))
  
  # Save mapping
  st_write(station_grid,
           file.path(CONFIG$dir_processed, "stations_with_grid.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
} else {
  message("No station data yet (run 02_collect_gios.R first)")
}

# --- 4. Join GUS demographic data to grid ------------------------------------

message("\n=== Joining GUS data to grid ===")
message("This requires gmina boundaries + BDL data from script 03.")
message("After 03_collect_gus.R is complete:")
message('  gminy_sf <- st_read("data/raw/gus/gadm41_POL_2.gpkg") |>')
message('    st_transform(CONFIG$crs_pl)')
message('  # Spatial join: assign gmina attributes to each grid cell')
message('  # (using the gmina that covers the cell centroid)')
message('  grid_centroids <- st_centroid(grid)')
message('  grid_with_gus <- st_join(grid_centroids, gminy_sf, join = st_within)')

# --- 5. Compute OSM-derived features ----------------------------------------

message("\n=== Computing OSM-derived spatial features ===")

osm_dir <- file.path(CONFIG$dir_raw, "osm")
roads_path <- file.path(osm_dir, "major_roads.gpkg")

if (file.exists(roads_path)) {
  message("Computing road density per grid cell (this may take a while)...")
  major_roads <- st_read(roads_path, quiet = TRUE)
  
  # Compute total road length intersecting each grid cell
  road_lengths <- st_intersection(major_roads, grid) |>
    mutate(length_m = as.numeric(st_length(geometry))) |>
    st_drop_geometry() |>
    group_by(cell_id) |>
    summarise(
      total_road_km = sum(length_m) / 1000,
      motorway_km   = sum(length_m[road_class == "motorway"]) / 1000,
      .groups = "drop"
    )
  
  grid <- grid |> left_join(road_lengths, by = "cell_id")
  grid$total_road_km[is.na(grid$total_road_km)] <- 0
  grid$motorway_km[is.na(grid$motorway_km)] <- 0
  
  message(glue("Road density computed for {nrow(grid)} cells"))
} else {
  message("No road data yet (run 04_collect_osm.R first)")
}

# Industrial area coverage
industrial_path <- file.path(osm_dir, "industrial_areas.gpkg")
if (file.exists(industrial_path)) {
  message("Computing industrial area coverage per grid cell...")
  industrial <- st_read(industrial_path, quiet = TRUE)
  
  ind_coverage <- st_intersection(industrial, grid) |>
    mutate(area_m2 = as.numeric(st_area(geometry))) |>
    st_drop_geometry() |>
    group_by(cell_id) |>
    summarise(industrial_area_m2 = sum(area_m2), .groups = "drop")
  
  grid <- grid |>
    left_join(ind_coverage, by = "cell_id") |>
    mutate(
      cell_area_m2 = as.numeric(st_area(geometry)),
      industrial_pct = coalesce(industrial_area_m2 / cell_area_m2 * 100, 0)
    )
  
  message("Industrial coverage computed")
}

# --- 6. Save the integrated grid ---------------------------------------------

message("\n=== Saving analysis grid ===")

st_write(grid, file.path(CONFIG$dir_processed, "analysis_grid.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

if (!is.null(tropomi_grid) && nrow(tropomi_grid) > 0) {
  write_csv(tropomi_grid,
            file.path(CONFIG$dir_processed, "tropomi_grid_values.csv"))
}

message(glue("\nGrid saved: {nrow(grid)} cells"))
message(glue("Columns: {paste(names(grid), collapse = ', ')}"))

message("\n=== Integration Summary ===")
message(glue("Grid:     {nrow(grid)} cells at {CONFIG$grid_res/1000}km"))
message(glue("CRS:      EPSG:{CONFIG$crs_pl}"))
message(glue("Output:   {CONFIG$dir_processed}/analysis_grid.gpkg"))
message(glue("          {CONFIG$dir_processed}/tropomi_grid_values.csv"))

# =============================================================================