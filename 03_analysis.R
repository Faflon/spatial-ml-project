# =============================================================================
# 03_analysis.R — EDA + Supervised ML + Unsupervised ML + Kriging
# =============================================================================
# Run after 02_build_grid.R has produced analysis_grid.gpkg.
# Each section is self-contained and saves outputs independently.
#
# Sections:
#   A. Exploratory Data Analysis
#   B. Supervised ML (RF, GW-RF, ANN)
#   C. Unsupervised ML (DBSCAN, LISA, clustering)
#   D. Kriging (Ordinary, Universal, Regression)
# =============================================================================

source("00_setup.R")

library(sf)
library(terra)
library(tidyverse)
library(spdep)
library(glue)
library(viridis)

# --- Load data ---------------------------------------------------------------

grid <- safe_read_sf(file.path(CONFIG$dir_processed, "analysis_grid.gpkg"))
tropomi <- safe_read_csv(file.path(CONFIG$dir_processed, "tropomi_grid.csv"))
stations <- safe_read_sf(file.path(CONFIG$dir_processed, "stations_gridded.gpkg"))

if (is.null(grid)) stop("Run 02_build_grid.R first to create analysis_grid.gpkg")


# #############################################################################
# A. EXPLORATORY DATA ANALYSIS                                               #
# #############################################################################

run_eda <- function() {

  library(corrplot)
  message("\n========== A. Exploratory Data Analysis ==========\n")

  # --- A1. NO2 map ---
  if ("tropomi_no2" %in% names(grid)) {
    message("Plotting NO2 spatial distribution...")

    p1 <- ggplot(grid |> filter(!is.na(tropomi_no2))) +
      geom_sf(aes(fill = tropomi_no2), color = NA) +
      scale_fill_viridis(option = "inferno",
                         name = expression(NO[2] ~ "(mol/m²)")) +
      labs(title = "Annual Mean Tropospheric NO₂ — Poland",
           subtitle = "Source: Sentinel-5P TROPOMI") +
      theme_minimal()

    ggsave(file.path(CONFIG$dir_figures, "no2_map.png"), p1,
           width = 10, height = 12, dpi = 300)
    message("  Saved: no2_map.png")
  }

  # --- A2. Monthly time series ---
  if (!is.null(tropomi)) {
    message("Plotting monthly time series...")

    monthly <- tropomi |>
      group_by(year_month, pollutant) |>
      summarise(
        mean_val = mean(value, na.rm = TRUE),
        sd_val   = sd(value, na.rm = TRUE),
        .groups  = "drop"
      )

    p2 <- ggplot(monthly, aes(x = year_month, y = mean_val)) +
      geom_line(linewidth = 1, color = "steelblue") +
      geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  alpha = 0.2, fill = "steelblue") +
      facet_wrap(~pollutant, scales = "free_y", ncol = 1) +
      labs(title = "Monthly Mean Pollutant Concentrations",
           x = NULL, y = "Column density (mol/m²)") +
      theme_minimal()

    ggsave(file.path(CONFIG$dir_figures, "monthly_ts.png"), p2,
           width = 10, height = 8, dpi = 300)
    message("  Saved: monthly_ts.png")
  }

  # --- A3. Feature distributions ---
  numeric_cols <- grid |>
    st_drop_geometry() |>
    select(where(is.numeric), -cell_id, -x_coord, -y_coord) |>
    names()

  if (length(numeric_cols) >= 2) {
    message("Plotting feature histograms...")

    grid_long <- grid |>
      st_drop_geometry() |>
      select(all_of(numeric_cols)) |>
      pivot_longer(everything(), names_to = "variable", values_to = "value")

    p3 <- ggplot(grid_long |> filter(!is.na(value)),
                 aes(x = value)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      facet_wrap(~variable, scales = "free", ncol = 3) +
      labs(title = "Feature Distributions") +
      theme_minimal()

    ggsave(file.path(CONFIG$dir_figures, "feature_distributions.png"), p3,
           width = 14, height = 10, dpi = 300)
    message("  Saved: feature_distributions.png")
  }

  # --- A4. Correlation matrix ---
  if (length(numeric_cols) >= 3) {
    message("Computing correlation matrix...")

    cor_data <- grid |>
      st_drop_geometry() |>
      select(all_of(numeric_cols)) |>
      drop_na()

    if (nrow(cor_data) > 30) {
      cor_mat <- cor(cor_data)
      png(file.path(CONFIG$dir_figures, "correlation_matrix.png"),
          width = 800, height = 800, res = 100)
      corrplot(cor_mat, method = "color", type = "lower",
               tl.cex = 0.7, addCoef.col = "black", number.cex = 0.6,
               title = "Feature Correlations", mar = base::c(0, 0, 2, 0))
      dev.off()
      message("  Saved: correlation_matrix.png")
    }
  }

  # --- A5. Spatial autocorrelation (Moran's I) ---
  if ("tropomi_no2" %in% names(grid)) {
    message("Testing spatial autocorrelation...")

    clean <- grid |> filter(!is.na(tropomi_no2))
    coords <- st_coordinates(st_centroid(clean))
    nb <- knearneigh(coords, k = 8) |> knn2nb()
    lw <- nb2listw(nb, style = "W")

    moran <- moran.test(clean$tropomi_no2, lw)
    message(sprintf("  Moran's I = %.4f, p = %s",
                    moran$estimate[1], format.pval(moran$p.value)))
  }

  # --- A6. Summary table ---
  message("\nSummary statistics:")
  summary_df <- grid |>
    st_drop_geometry() |>
    select(all_of(numeric_cols)) |>
    pivot_longer(everything(), names_to = "variable", values_to = "value") |>
    group_by(variable) |>
    summarise(
      n      = sum(!is.na(value)),
      mean   = mean(value, na.rm = TRUE),
      sd     = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      min    = min(value, na.rm = TRUE),
      max    = max(value, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(summary_df, file.path(CONFIG$dir_output, "summary_stats.csv"))
  print(summary_df)

  message("\nEDA complete. Figures in: ", CONFIG$dir_figures)
}


# #############################################################################
# B. SUPERVISED MACHINE LEARNING                                              #
# #############################################################################

run_supervised <- function() {

  library(ranger)
  library(caret)
  library(nnet)

  message("\n========== B. Supervised ML ==========\n")

  # --- B1. Build modelling dataset ---
  # Target: ground-level NO2 from GIOŚ (joined to grid)
  # Features: TROPOMI, OSM, GUS, coordinates

  if (is.null(stations)) {
    stop("No station data. Run collect_gios() and 02_build_grid.R first.")
  }

  # Get station measurements (recent or from bulk parse)
  meas_path <- file.path(CONFIG$dir_raw, "gios", "recent_measurements.csv")
  measurements <- safe_read_csv(meas_path)

  if (is.null(measurements) || nrow(measurements) == 0) {
    stop("No measurement data found. Run collect_gios() first.")
  }

  # Average NO2 per station
  station_no2 <- measurements |>
    filter(param_code == "NO2") |>
    group_by(station_id) |>
    summarise(ground_no2 = mean(value, na.rm = TRUE), .groups = "drop")

  # Join to grid-assigned stations
  model_stations <- stations |>
    st_drop_geometry() |>
    inner_join(station_no2, by = "station_id") |>
    filter(!is.na(cell_id))

  # Join grid features
  grid_features <- grid |>
    st_drop_geometry() |>
    select(-x_coord, -y_coord)

  model_df <- model_stations |>
    inner_join(grid_features, by = "cell_id") |>
    drop_na(ground_no2)

  # Add coordinates
  station_coords <- stations |>
    filter(station_id %in% model_df$station_id)
  coords_mat <- st_coordinates(station_coords)
  model_df$lon <- coords_mat[match(model_df$station_id, station_coords$station_id), 1]
  model_df$lat <- coords_mat[match(model_df$station_id, station_coords$station_id), 2]

  message(glue("Modelling dataset: {nrow(model_df)} station-observations"))

  # Feature columns (exclude IDs and target)
  exclude_cols <- base::c("station_id", "cell_id", "ground_no2",
                           "name", "city", "commune", "province")
  feature_cols <- setdiff(names(model_df), exclude_cols)
  feature_cols <- feature_cols[sapply(model_df[feature_cols], is.numeric)]

  message(glue("Features ({length(feature_cols)}): {paste(feature_cols, collapse = ', ')}"))

  if (length(feature_cols) < 2 || nrow(model_df) < 20) {
    stop("Not enough features or data points for modelling.")
  }

  # --- B2. Spatial train/test split ---
  set.seed(42)
  n_blocks <- 5
  km <- kmeans(model_df[, base::c("lon", "lat")], centers = n_blocks)
  model_df$block <- km$cluster
  test_block <- sample(unique(model_df$block), 1)

  train_df <- model_df |> filter(block != test_block)
  test_df  <- model_df |> filter(block == test_block)

  message(glue("Train: {nrow(train_df)}, Test: {nrow(test_df)} (spatial block CV)"))

  formula_str <- paste("ground_no2 ~", paste(feature_cols, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # --- B3. Random Forest ---
  message("\nTraining Random Forest...")

  rf <- ranger(
    formula    = formula_obj,
    data       = train_df,
    num.trees  = 500,
    importance = "impurity",
    seed       = 42
  )

  message(sprintf("  OOB R² = %.4f, RMSE = %.4f",
                  rf$r.squared, sqrt(rf$prediction.error)))

  # Test set evaluation
  rf_pred <- predict(rf, data = test_df)$predictions
  rf_metrics <- calc_metrics(test_df$ground_no2, rf_pred)
  message(sprintf("  Test R² = %.4f, RMSE = %.4f, MAE = %.4f",
                  rf_metrics$R2, rf_metrics$RMSE, rf_metrics$MAE))

  # Variable importance
  imp_df <- tibble(
    variable   = names(rf$variable.importance),
    importance = rf$variable.importance
  ) |> arrange(desc(importance))

  p_imp <- ggplot(head(imp_df, 15),
                  aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance (Random Forest)", x = NULL) +
    theme_minimal()

  ggsave(file.path(CONFIG$dir_figures, "rf_importance.png"), p_imp,
         width = 8, height = 6, dpi = 300)

  # --- B4. Geographically Weighted RF ---
  message("\nTraining GW-RF...")

  gwrf <- NULL
  if (!requireNamespace("SpatialML", quietly = TRUE)) {
    message("  SpatialML not installed — skipping GW-RF")
    message("  To install: remotes::install_github('SpatialML/SpatialML')")
  } else {
    library(SpatialML)
    tryCatch({
      train_coords <- as.matrix(train_df[, base::c("lon", "lat")])

      gwrf <- grf(
        formula = formula_obj,
        dframe  = as.data.frame(train_df),
        bw      = NULL,
        kernel  = "adaptive",
        coords  = train_coords,
        ntree   = 200
      )

      message(sprintf("  GW-RF bandwidth = %s", gwrf$Bandwidth))
    }, error = function(e) {
      warning(paste0("  GW-RF failed: ", e$message))
      gwrf <<- NULL
    })
  }  # end SpatialML check

  # --- B5. Neural Network ---
  message("\nTraining ANN...")

  tryCatch({
    preproc <- preProcess(train_df[feature_cols], method = base::c("center", "scale"))
    train_scaled <- predict(preproc, train_df)
    test_scaled  <- predict(preproc, test_df)

    ctrl <- trainControl(method = "cv", number = 5)

    ann <- train(
      formula_obj,
      data      = train_scaled,
      method    = "nnet",
      trControl = ctrl,
      tuneGrid  = expand.grid(size = base::c(8, 16), decay = base::c(0.01, 0.001)),
      linout    = TRUE,
      maxit     = 500,
      trace     = FALSE
    )

    ann_pred <- predict(ann, newdata = test_scaled)
    ann_metrics <- calc_metrics(test_df$ground_no2, ann_pred)
    message(sprintf("  ANN Test R² = %.4f, RMSE = %.4f",
                    ann_metrics$R2, ann_metrics$RMSE))
  }, error = function(e) {
    warning(paste0("  ANN failed: ", e$message))
    ann_metrics <- NULL
  })

  # --- B6. Model comparison ---
  comparison <- bind_rows(
    RF  = rf_metrics,
    ANN = if (exists("ann_metrics") && !is.null(ann_metrics)) ann_metrics else NULL,
    .id = "model"
  )

  write_csv(comparison, file.path(CONFIG$dir_output, "model_comparison.csv"))
  message("\nModel comparison:")
  print(comparison)

  # --- B7. Predict on full grid ---
  message("\nPredicting on full grid...")

  grid_data <- grid |> st_drop_geometry()
  # Check which features exist in grid
  missing <- setdiff(feature_cols, names(grid_data))
  if (length(missing) > 0) {
    message(glue("  Warning: missing features in grid: {paste(missing, collapse=', ')}"))
    message("  Skipping full grid prediction.")
  } else {
    grid$predicted_no2 <- predict(rf, data = grid_data)$predictions
    st_write(grid, file.path(CONFIG$dir_output, "predicted_no2.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message("  Saved: predicted_no2.gpkg")
  }

  # Return objects for use in kriging
  base::list(rf = rf, model_df = model_df, train_df = train_df,
             test_df = test_df, feature_cols = feature_cols)
}


# #############################################################################
# C. UNSUPERVISED MACHINE LEARNING                                            #
# #############################################################################

run_unsupervised <- function() {

  library(dbscan)
  library(cluster)

  message("\n========== C. Unsupervised ML ==========\n")

  if (!"tropomi_no2" %in% names(grid)) {
    message("No TROPOMI data in grid — skipping unsupervised analysis")
    return(NULL)
  }

  # --- C1. Hotspot detection (DBSCAN) ---
  message("Hotspot clustering (DBSCAN)...")

  clean <- grid |> filter(!is.na(tropomi_no2))
  threshold <- quantile(clean$tropomi_no2, 0.75, na.rm = TRUE)
  hot_cells <- clean |> filter(tropomi_no2 >= threshold)

  coords <- st_coordinates(st_centroid(hot_cells))
  db <- dbscan(coords, eps = 15000, minPts = 5)

  hot_cells$cluster <- db$cluster
  n_clusters <- length(unique(db$cluster[db$cluster > 0]))
  n_noise <- sum(db$cluster == 0)
  message(sprintf("  Found %d hotspot clusters (%d noise points)", n_clusters, n_noise))

  p_hot <- ggplot() +
    geom_sf(data = clean, fill = "grey95", color = NA) +
    geom_sf(data = hot_cells |> filter(cluster > 0),
            aes(fill = factor(cluster)), color = NA, alpha = 0.7) +
    scale_fill_brewer(palette = "Set1", name = "Cluster") +
    labs(title = "NO₂ Pollution Hotspot Clusters (DBSCAN)") +
    theme_minimal()

  ggsave(file.path(CONFIG$dir_figures, "hotspots_dbscan.png"), p_hot,
         width = 10, height = 12, dpi = 300)

  # --- C2. LISA (Local Moran's I) ---
  message("\nLISA analysis...")

  coords_all <- st_coordinates(st_centroid(clean))
  nb <- knearneigh(coords_all, k = 8) |> knn2nb()
  lw <- nb2listw(nb, style = "W")

  local_m <- localmoran(clean$tropomi_no2, lw)

  val_z <- scale(clean$tropomi_no2)[, 1]
  lag_z <- lag.listw(lw, val_z)

  clean$lisa_class <- case_when(
    local_m[, 5] > 0.05        ~ "Not significant",
    val_z > 0 & lag_z > 0      ~ "High-High",
    val_z < 0 & lag_z < 0      ~ "Low-Low",
    val_z > 0 & lag_z < 0      ~ "High-Low",
    val_z < 0 & lag_z > 0      ~ "Low-High"
  )

  lisa_colors <- base::c(
    "High-High"       = "#d7191c",
    "Low-Low"         = "#2c7bb6",
    "High-Low"        = "#fdae61",
    "Low-High"        = "#abd9e9",
    "Not significant" = "grey90"
  )

  p_lisa <- ggplot(clean) +
    geom_sf(aes(fill = lisa_class), color = NA) +
    scale_fill_manual(values = lisa_colors, name = "LISA") +
    labs(title = "LISA Cluster Map — Tropospheric NO₂") +
    theme_minimal()

  ggsave(file.path(CONFIG$dir_figures, "lisa_no2.png"), p_lisa,
         width = 10, height = 12, dpi = 300)
  message("  Saved: lisa_no2.png")

  # --- C3. Multi-variable clustering (PAM) ---
  cluster_cols <- intersect(
    base::c("tropomi_no2", "total_road_km", "industrial_pct", "railway_km"),
    names(grid)
  )

  if (length(cluster_cols) >= 2) {
    message("\nMulti-variable spatial clustering...")

    clust_data <- clean |>
      st_drop_geometry() |>
      select(cell_id, all_of(cluster_cols)) |>
      drop_na()

    if (nrow(clust_data) >= 50) {
      data_scaled <- scale(clust_data |> select(-cell_id))

      # Optimal k via silhouette
      sil_scores <- sapply(2:8, function(k) {
        pam_fit <- pam(data_scaled, k = k)
        mean(silhouette(pam_fit)[, "sil_width"])
      })

      best_k <- (2:8)[which.max(sil_scores)]
      message(sprintf("  Best k = %d (silhouette = %.3f)", best_k, max(sil_scores)))

      pam_final <- pam(data_scaled, k = best_k)
      clust_data$cluster <- pam_final$clustering

      clean_clustered <- clean |>
        left_join(clust_data |> select(cell_id, cluster), by = "cell_id")

      p_clust <- ggplot(clean_clustered |> filter(!is.na(cluster))) +
        geom_sf(aes(fill = factor(cluster)), color = NA) +
        scale_fill_brewer(palette = "Set2", name = "Cluster") +
        labs(title = sprintf("Spatial Typology (PAM k=%d)", best_k)) +
        theme_minimal()

      ggsave(file.path(CONFIG$dir_figures, "spatial_clusters.png"), p_clust,
             width = 10, height = 12, dpi = 300)
    }
  }

  message("\nUnsupervised ML complete.")
}


# #############################################################################
# D. KRIGING                                                                  #
# #############################################################################

run_kriging <- function(supervised_results = NULL) {

  library(gstat)
  library(automap)

  message("\n========== D. Kriging ==========\n")

  if (is.null(stations)) {
    stop("No station data for kriging.")
  }

  # Load station measurements
  meas <- safe_read_csv(file.path(CONFIG$dir_raw, "gios", "recent_measurements.csv"))
  if (is.null(meas)) stop("No measurement data.")

  # Build station point layer with mean NO2
  station_no2 <- meas |>
    filter(param_code == "NO2") |>
    group_by(station_id) |>
    summarise(no2 = mean(value, na.rm = TRUE), .groups = "drop")

  station_pts <- stations |>
    inner_join(station_no2, by = "station_id") |>
    filter(!is.na(no2))

  message(glue("Kriging with {nrow(station_pts)} station points"))

  # --- D1. Variogram ---
  message("\nFitting variogram...")

  vario <- autofitVariogram(no2 ~ 1, input_data = station_pts)
  message(sprintf("  Model: %s | Nugget=%.1f Sill=%.1f Range=%.0f m",
                  vario$var_model$model[2],
                  vario$var_model$psill[1],
                  sum(vario$var_model$psill),
                  vario$var_model$range[2]))

  png(file.path(CONFIG$dir_figures, "variogram.png"), width = 600, height = 400)
  plot(vario, main = "Semivariogram — Ground-Level NO₂")
  dev.off()
  message("  Saved: variogram.png")

  # --- D2. Ordinary Kriging ---
  message("\nOrdinary Kriging...")

  pred_pts <- st_centroid(grid)

  ok_result <- krige(
    no2 ~ 1,
    locations = station_pts,
    newdata   = pred_pts,
    model     = vario$var_model
  )

  grid$ok_pred <- ok_result$var1.pred
  grid$ok_se   <- sqrt(ok_result$var1.var)

  p_ok <- ggplot(grid |> filter(!is.na(ok_pred))) +
    geom_sf(aes(fill = ok_pred), color = NA) +
    scale_fill_viridis(option = "inferno", name = "NO₂ (µg/m³)") +
    labs(title = "Ordinary Kriging — Ground-Level NO₂") +
    theme_minimal()

  ggsave(file.path(CONFIG$dir_figures, "kriging_ok.png"), p_ok,
         width = 10, height = 12, dpi = 300)

  p_se <- ggplot(grid |> filter(!is.na(ok_se))) +
    geom_sf(aes(fill = ok_se), color = NA) +
    scale_fill_viridis(option = "magma", name = "Std Error") +
    labs(title = "Kriging Uncertainty") +
    theme_minimal()

  ggsave(file.path(CONFIG$dir_figures, "kriging_uncertainty.png"), p_se,
         width = 10, height = 12, dpi = 300)

  # --- D3. Cross-validation ---
  message("\nKriging cross-validation...")

  cv_ok <- krige.cv(no2 ~ 1, locations = station_pts,
                    model = vario$var_model, nfold = 5)

  ok_cv_metrics <- calc_metrics(cv_ok$observed, cv_ok$observed - cv_ok$residual)
  message(sprintf("  OK CV: R²=%.4f RMSE=%.4f MAE=%.4f",
                  ok_cv_metrics$R2, ok_cv_metrics$RMSE, ok_cv_metrics$MAE))

  # --- D4. Universal Kriging ---
  message("\nUniversal Kriging (with coordinate trend)...")

  tryCatch({
    coords_st <- st_coordinates(station_pts)
    station_pts$x <- coords_st[, 1]
    station_pts$y <- coords_st[, 2]

    coords_gr <- st_coordinates(pred_pts)
    pred_pts$x <- coords_gr[, 1]
    pred_pts$y <- coords_gr[, 2]

    vario_uk <- autofitVariogram(no2 ~ x + y, input_data = station_pts)

    uk_result <- krige(
      no2 ~ x + y,
      locations = station_pts,
      newdata   = pred_pts,
      model     = vario_uk$var_model
    )

    grid$uk_pred <- uk_result$var1.pred
    message("  Universal Kriging: done")

    # UK cross-validation
    cv_uk <- krige.cv(no2 ~ x + y, locations = station_pts,
                      model = vario_uk$var_model, nfold = 5)
    uk_cv_metrics <- calc_metrics(cv_uk$observed, cv_uk$observed - cv_uk$residual)
    message(sprintf("  UK CV: R²=%.4f RMSE=%.4f", uk_cv_metrics$R2, uk_cv_metrics$RMSE))
  }, error = function(e) {
    warning(paste0("  UK failed: ", e$message))
  })

  # --- D5. Regression Kriging (RF + kriging of residuals) ---
  if (!is.null(supervised_results)) {
    message("\nRegression Kriging...")

    tryCatch({
      rf <- supervised_results$rf
      feature_cols <- supervised_results$feature_cols

      # RF predictions at station locations
      st_data <- station_pts |> st_drop_geometry()
      available_features <- intersect(feature_cols, names(st_data))

      if (length(available_features) > 0) {
        st_data_joined <- station_pts |>
          st_join(grid |> select(cell_id, all_of(available_features))) |>
          st_drop_geometry() |>
          drop_na()

        # This is a placeholder — full implementation requires matching
        # station data to grid features properly
        message("  Regression Kriging: requires feature alignment (TODO)")
      }
    }, error = function(e) {
      warning(paste0("  Regression Kriging failed: ", e$message))
    })
  }

  # --- D6. Save ---
  st_write(grid, file.path(CONFIG$dir_output, "kriging_results.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)

  kriging_comparison <- bind_rows(
    `Ordinary Kriging`  = ok_cv_metrics,
    `Universal Kriging` = if (exists("uk_cv_metrics")) uk_cv_metrics else NULL,
    .id = "method"
  )

  write_csv(kriging_comparison,
            file.path(CONFIG$dir_output, "kriging_comparison.csv"))

  message("\nKriging comparison:")
  print(kriging_comparison)
  message("\nKriging complete.")
}


# #############################################################################
# RUN ALL                                                                     #
# #############################################################################
# Uncomment sections as data becomes available:

# run_eda()
# sup_results <- run_supervised()
# run_unsupervised()
# run_kriging(supervised_results = sup_results)

message("\n=== 03_analysis.R ready ===")
message("Uncomment run_*() calls above to execute each section.")

# =============================================================================
