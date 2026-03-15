# =============================================================================
# 08_unsupervised_ml.R — Unsupervised Spatial Machine Learning
# =============================================================================
# Methods: DBSCAN/HDBSCAN for hotspot clustering, spatial agglomeration
#          measurement, comparison of spatial distributions
# Input:   Predicted pollution surface + TROPOMI grid data
# Output:  Cluster maps, agglomeration indices, distribution tests
# =============================================================================

source("00_setup.R")
source("10_utils.R")

library(sf)
library(tidyverse)
library(spdep)
library(dbscan)       # DBSCAN / HDBSCAN clustering
library(cluster)      # silhouette, PAM
library(spatstat)     # point pattern analysis
library(glue)

# --- 1. Load data ------------------------------------------------------------

grid <- st_read(file.path(CONFIG$dir_processed, "analysis_grid.gpkg"), quiet = TRUE)
tropomi <- read_csv(file.path(CONFIG$dir_processed, "tropomi_grid_values.csv"),
                    show_col_types = FALSE)

# Merge annual mean NO2 into grid
no2_annual <- tropomi |>
  filter(pollutant == "NO2") |>
  group_by(cell_id) |>
  summarise(no2_mean = mean(value, na.rm = TRUE), .groups = "drop")

grid <- grid |> left_join(no2_annual, by = "cell_id")

# --- 2. Pollution hotspot detection (DBSCAN on high-NO2 cells) ---------------

message("=== Pollution hotspot clustering ===")

identify_hotspots <- function(grid_sf, value_col = "no2_mean",
                              threshold_quantile = 0.75,
                              eps = 15000, minPts = 5) {
  #' Identify spatial clusters of high-pollution grid cells using DBSCAN.
  #' 1. Filter cells above a pollution threshold
  #' 2. Cluster the high-pollution cells spatially
  #' Returns grid_sf with cluster labels.
  
  threshold <- quantile(grid_sf[[value_col]], threshold_quantile, na.rm = TRUE)
  
  hot_cells <- grid_sf |>
    filter(.data[[value_col]] >= threshold) |>
    filter(!is.na(.data[[value_col]]))
  
  message(glue("  Cells above Q{threshold_quantile*100}: {nrow(hot_cells)}"))
  
  # Get centroids for clustering
  coords <- st_coordinates(st_centroid(hot_cells))
  
  # DBSCAN clustering
  db <- dbscan(coords, eps = eps, minPts = minPts)
  
  hot_cells$cluster <- db$cluster
  hot_cells$is_noise <- db$cluster == 0
  
  n_clusters <- length(unique(db$cluster[db$cluster > 0]))
  message(glue("  Found {n_clusters} hotspot clusters ({sum(db$cluster == 0)} noise points)"))
  
  return(hot_cells)
}

# hotspots <- identify_hotspots(grid, eps = 15000, minPts = 5)

# Plot hotspot clusters
plot_hotspots <- function(grid_sf, hotspots_sf) {
  ggplot() +
    geom_sf(data = grid_sf, fill = "grey95", color = NA) +
    geom_sf(data = hotspots_sf |> filter(!is_noise),
            aes(fill = factor(cluster)), color = NA, alpha = 0.7) +
    geom_sf(data = hotspots_sf |> filter(is_noise),
            fill = "grey50", color = NA, alpha = 0.3) +
    scale_fill_brewer(palette = "Set1", name = "Cluster") +
    labs(title = "NO₂ Pollution Hotspot Clusters (DBSCAN)",
         subtitle = "Top 25% grid cells, clustered spatially") +
    theme_minimal()
}

# --- 3. Multi-variable spatial clustering ------------------------------------

message("=== Multi-variable spatial clustering ===")

spatial_clustering <- function(grid_sf, features, n_clusters = 5,
                               include_coords = TRUE) {
  #' K-medoids clustering on pollution + infrastructure features.
  #' Optionally includes coordinates to enforce spatial contiguity.
  
  data <- grid_sf |>
    st_drop_geometry() |>
    select(all_of(features)) |>
    drop_na()
  
  if (include_coords) {
    coords <- st_coordinates(st_centroid(grid_sf |> filter(!is.na(!!sym(features[1])))))
    data$x <- coords[, 1]
    data$y <- coords[, 2]
  }
  
  # Scale features
  data_scaled <- scale(data)
  
  # PAM (Partitioning Around Medoids) — more robust than k-means
  pam_result <- pam(data_scaled, k = n_clusters)
  
  grid_sf$spatial_cluster <- NA
  valid_idx <- which(!is.na(grid_sf[[features[1]]]))
  grid_sf$spatial_cluster[valid_idx] <- pam_result$clustering
  
  # Silhouette score
  sil <- silhouette(pam_result)
  avg_sil <- mean(sil[, "sil_width"])
  message(glue("  PAM k={n_clusters}: avg silhouette = {round(avg_sil, 3)}"))
  
  return(list(grid = grid_sf, pam = pam_result, silhouette = avg_sil))
}

# Optimal k selection
find_optimal_k <- function(grid_sf, features, k_range = 2:10) {
  #' Test multiple k values and return silhouette scores.
  
  scores <- sapply(k_range, function(k) {
    result <- spatial_clustering(grid_sf, features, n_clusters = k)
    result$silhouette
  })
  
  tibble(k = k_range, silhouette = scores)
}

# cluster_features <- c("no2_mean", "total_road_km", "industrial_pct")
# optimal_k <- find_optimal_k(grid, cluster_features)
# best_k <- optimal_k$k[which.max(optimal_k$silhouette)]
# final_clusters <- spatial_clustering(grid, cluster_features, n_clusters = best_k)

# --- 4. Measuring spatial agglomeration (Moran's I, LISA) --------------------

message("=== Spatial agglomeration measurement ===")

measure_agglomeration <- function(grid_sf, value_col, k = 8) {
  #' Compute global and local Moran's I for a variable.
  #' Returns global test + LISA cluster map data.
  
  clean_grid <- grid_sf |> filter(!is.na(.data[[value_col]]))
  coords <- st_coordinates(st_centroid(clean_grid))
  
  # Spatial weights (k-nearest neighbors)
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W")
  
  # Global Moran's I
  global_moran <- moran.test(clean_grid[[value_col]], lw)
  message(glue("  Global Moran's I = {round(global_moran$estimate[1], 4)}, ",
               "p = {format.pval(global_moran$p.value)}"))
  
  # Local Moran's I (LISA) — identifies local clusters & outliers
  local_moran <- localmoran(clean_grid[[value_col]], lw)
  
  clean_grid$lisa_i <- local_moran[, 1]
  clean_grid$lisa_p <- local_moran[, 5]
  
  # Classify LISA clusters (HH, HL, LH, LL)
  val_scaled <- scale(clean_grid[[value_col]])[, 1]
  lag_scaled <- lag.listw(lw, val_scaled)
  
  clean_grid$lisa_class <- case_when(
    local_moran[, 5] > 0.05          ~ "Not significant",
    val_scaled > 0 & lag_scaled > 0   ~ "High-High",
    val_scaled < 0 & lag_scaled < 0   ~ "Low-Low",
    val_scaled > 0 & lag_scaled < 0   ~ "High-Low",
    val_scaled < 0 & lag_scaled > 0   ~ "Low-High"
  )
  
  return(list(
    global = global_moran,
    local  = clean_grid
  ))
}

# moran_no2 <- measure_agglomeration(grid, "no2_mean")

# LISA cluster map
plot_lisa <- function(lisa_grid, title = "LISA Cluster Map") {
  lisa_colors <- c(
    "High-High"       = "#d7191c",
    "Low-Low"         = "#2c7bb6",
    "High-Low"        = "#fdae61",
    "Low-High"        = "#abd9e9",
    "Not significant" = "grey90"
  )
  
  ggplot(lisa_grid) +
    geom_sf(aes(fill = lisa_class), color = NA) +
    scale_fill_manual(values = lisa_colors, name = "LISA Cluster") +
    labs(title = title) +
    theme_minimal()
}

# --- 5. Comparison of spatial distributions ----------------------------------

message("=== Spatial distribution comparison ===")

compare_distributions <- function(grid_sf, var1, var2) {
  #' Compare spatial distributions of two variables.
  #' Uses: spatial correlation, Lee's L statistic, visual overlay.
  
  clean <- grid_sf |>
    st_drop_geometry() |>
    select(all_of(c(var1, var2))) |>
    drop_na()
  
  # Pearson correlation
  r <- cor(clean[[var1]], clean[[var2]])
  message(glue("  Pearson r({var1}, {var2}) = {round(r, 4)}"))
  
  # Spatial correlation (bivariate Moran's I)
  # coords <- st_coordinates(st_centroid(grid_sf |> drop_na()))
  # nb <- knn2nb(knearneigh(coords, k = 8))
  # lw <- nb2listw(nb, style = "W")
  # bv_moran <- moran_bv(clean[[var1]], clean[[var2]], lw)
  
  return(list(pearson_r = r))
}

# Compare NO2 vs road density, NO2 vs industrial area, etc.
# compare_distributions(grid, "no2_mean", "total_road_km")
# compare_distributions(grid, "no2_mean", "industrial_pct")

# --- 6. Save results ---------------------------------------------------------

message("\n=== Unsupervised ML Complete ===")
message("Key outputs: hotspot clusters, LISA maps, agglomeration stats")

# =============================================================================