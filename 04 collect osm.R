# =============================================================================
# 04_collect_osm.R — OpenStreetMap Spatial Features
# =============================================================================
# Data source: OpenStreetMap via Overpass API (osmdata package)
# What:        Roads, industrial sites, railways, fuel stations, etc.
# Output:      sf objects saved as .gpkg for spatial feature engineering
# =============================================================================

source("00_setup.R")

library(osmdata)
library(sf)
library(tidyverse)
library(glue)

# --- Configuration -----------------------------------------------------------

osm_dir <- file.path(CONFIG$dir_raw, "osm")

# Poland bounding box for queries
bb <- CONFIG$bbox_poland
poland_bb <- opq(bbox = c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]))

# NOTE: Large Overpass queries for all of Poland can time out.
# Strategy: query by województwo or use a larger timeout.
# The osmdata package default timeout is 25s — we increase it.

# --- Helper: query with timeout and caching ----------------------------------

osm_query <- function(feature_key, feature_value = NULL, geometry_type = "osm_lines",
                      timeout = 180) {
  #' Query OSM data for Poland with extended timeout.
  #' Returns sf object or NULL on failure.
  
  q <- poland_bb |>
    opq(timeout = timeout)
  
  if (is.null(feature_value)) {
    q <- q |> add_osm_feature(key = feature_key)
  } else {
    q <- q |> add_osm_feature(key = feature_key, value = feature_value)
  }
  
  message(glue("  Querying OSM: {feature_key}={feature_value %||% '*'}..."))
  
  result <- tryCatch({
    osmdata_sf(q)
  }, error = function(e) {
    warning(glue("  OSM query failed: {e$message}"))
    return(NULL)
  })
  
  if (is.null(result)) return(NULL)
  
  # Extract the requested geometry type
  geom <- result[[geometry_type]]
  
  if (is.null(geom) || nrow(geom) == 0) {
    message(glue("  No {geometry_type} found"))
    return(NULL)
  }
  
  message(glue("  Got {nrow(geom)} features"))
  return(geom)
}

# --- 1. Major roads ----------------------------------------------------------
# Road density is a strong predictor of NO2 (traffic emissions)

message("\n=== 1. Roads ===")

# Motorways and trunk roads (highest traffic / emissions)
motorways <- osm_query("highway", "motorway", "osm_lines")
trunk_roads <- osm_query("highway", "trunk", "osm_lines")
primary_roads <- osm_query("highway", "primary", "osm_lines")
secondary_roads <- osm_query("highway", "secondary", "osm_lines")

# Combine all major roads
major_roads <- bind_rows(
  motorways |> mutate(road_class = "motorway") |> select(osm_id, road_class, geometry),
  trunk_roads |> mutate(road_class = "trunk") |> select(osm_id, road_class, geometry),
  primary_roads |> mutate(road_class = "primary") |> select(osm_id, road_class, geometry),
  secondary_roads |> mutate(road_class = "secondary") |> select(osm_id, road_class, geometry)
) |>
  st_transform(CONFIG$crs_pl)

st_write(major_roads, file.path(osm_dir, "major_roads.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
message(glue("Saved {nrow(major_roads)} road segments"))

# --- 2. Railways -------------------------------------------------------------

message("\n=== 2. Railways ===")

railways <- osm_query("railway", "rail", "osm_lines")
if (!is.null(railways)) {
  railways <- railways |>
    select(osm_id, geometry) |>
    st_transform(CONFIG$crs_pl)
  st_write(railways, file.path(osm_dir, "railways.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  message(glue("Saved {nrow(railways)} railway segments"))
}

# --- 3. Industrial sites -----------------------------------------------------
# Industrial land use = potential emission sources

message("\n=== 3. Industrial areas ===")

industrial <- osm_query("landuse", "industrial", "osm_polygons")
if (!is.null(industrial)) {
  industrial <- industrial |>
    select(osm_id, name, geometry) |>
    st_transform(CONFIG$crs_pl)
  st_write(industrial, file.path(osm_dir, "industrial_areas.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  message(glue("Saved {nrow(industrial)} industrial polygons"))
}

# --- 4. Power plants & factories (point sources) ----------------------------

message("\n=== 4. Power plants ===")

power_plants <- osm_query("power", "plant", "osm_polygons")
if (!is.null(power_plants)) {
  # Convert polygons to centroids for point-based analysis
  pp_points <- power_plants |>
    select(osm_id, name, geometry) |>
    st_transform(CONFIG$crs_pl) |>
    st_centroid()
  st_write(pp_points, file.path(osm_dir, "power_plants.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  message(glue("Saved {nrow(pp_points)} power plants"))
}

# --- 5. Residential / urban areas -------------------------------------------

message("\n=== 5. Residential areas ===")

residential <- osm_query("landuse", "residential", "osm_polygons", timeout = 300)
if (!is.null(residential)) {
  residential <- residential |>
    select(osm_id, geometry) |>
    st_transform(CONFIG$crs_pl)
  st_write(residential, file.path(osm_dir, "residential_areas.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)
  message(glue("Saved {nrow(residential)} residential polygons"))
}

# --- 6. Feature engineering helper -------------------------------------------
# After collecting raw OSM, we need to compute density metrics per grid cell.
# This will be done in 05_build_grid.R, but here's the logic:

compute_road_density <- function(roads_sf, grid_sf, buffer_m = 5000) {
  #' Compute total road length (km) within buffer of each grid cell centroid.
  #' Returns grid_sf with added road_density_km column.
  
  centroids <- st_centroid(grid_sf)
  buffers <- st_buffer(centroids, buffer_m)
  
  road_lengths <- sapply(seq_len(nrow(buffers)), function(i) {
    clipped <- st_intersection(roads_sf, buffers[i, ])
    if (nrow(clipped) == 0) return(0)
    sum(st_length(clipped)) / 1000  # meters → km
  })
  
  grid_sf$road_density_km <- as.numeric(road_lengths)
  return(grid_sf)
}

compute_point_count <- function(points_sf, grid_sf, buffer_m = 10000) {
  #' Count points (e.g., power plants) within buffer of each grid cell.
  
  centroids <- st_centroid(grid_sf)
  buffers <- st_buffer(centroids, buffer_m)
  
  counts <- lengths(st_intersects(buffers, points_sf))
  grid_sf$point_count <- counts
  return(grid_sf)
}

compute_area_coverage <- function(polygons_sf, grid_sf) {
  #' Compute % of each grid cell covered by polygons (e.g., industrial area).
  
  coverage <- sapply(seq_len(nrow(grid_sf)), function(i) {
    inter <- st_intersection(polygons_sf, grid_sf[i, ])
    if (nrow(inter) == 0) return(0)
    sum(st_area(inter)) / st_area(grid_sf[i, ])
  })
  
  grid_sf$area_coverage_pct <- as.numeric(coverage) * 100
  return(grid_sf)
}

message("\n=== OSM Collection Complete ===")
message(glue("Files saved in: {osm_dir}"))
message("Feature engineering functions defined for use in 05_build_grid.R")

# =============================================================================
# NOTE FOR LARGE QUERIES:
# If Overpass API times out for all of Poland, split by voivodship:
#
#   voivodeships <- list(
#     mazowieckie  = c(19.3, 51.0, 22.0, 53.0),
#     malopolskie  = c(19.1, 49.2, 21.3, 50.5),
#     # ... etc
#   )
#
#   for (name in names(voivodeships)) {
#     bb <- voivodeships[[name]]
#     roads <- opq(bbox = bb, timeout = 120) |>
#       add_osm_feature("highway", "motorway") |>
#       osmdata_sf()
#     # ... save per region, then merge
#   }
#
# Alternative: download Poland OSM extract from Geofabrik:
#   https://download.geofabrik.de/europe/poland.html
#   Then read with sf::st_read("poland-latest-free.shp/...")
# =============================================================================