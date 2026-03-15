# =============================================================================
# 07_supervised_ml.R — Supervised Spatial Machine Learning
# =============================================================================
# Models: Random Forest, Geographically Weighted RF, ANN
# Target: Ground-level NO2 concentration (from GIOŚ stations)
# Features: TROPOMI columns, OSM infrastructure, GUS demographics, coords
# Output: Model objects, predictions, importance plots, accuracy metrics
# =============================================================================

source("00_setup.R")
source("10_utils.R")

library(sf)
library(tidyverse)
library(ranger)       # fast RF
library(caret)        # model training framework
library(SpatialML)    # geographically weighted RF (grf)
library(spdep)        # spatial weights
library(nnet)         # basic ANN (single hidden layer)
library(glue)

# --- 1. Prepare modelling dataset --------------------------------------------

message("=== Preparing modelling data ===")

grid <- st_read(file.path(CONFIG$dir_processed, "analysis_grid.gpkg"), quiet = TRUE)
stations <- st_read(file.path(CONFIG$dir_processed, "stations_with_grid.gpkg"),
                    quiet = TRUE)

# The modelling dataset = grid cells that contain a GIOŚ station (supervised)
# Target: monthly mean ground-level NO2 from GIOŚ
# Features: TROPOMI NO2, road density, industrial %, population, coords, etc.

# --- Placeholder: construct model_df after all data is integrated ---
# model_df should look like:
#   cell_id | year_month | ground_no2 | tropomi_no2 | tropomi_so2 | tropomi_co |
#   total_road_km | motorway_km | industrial_pct | population_density |
#   lon | lat | ...

# model_df <- stations_monthly |>
#   left_join(tropomi_values, by = c("cell_id", "year_month")) |>
#   left_join(grid |> st_drop_geometry(), by = "cell_id") |>
#   mutate(
#     lon = st_coordinates(geometry)[, 1],
#     lat = st_coordinates(geometry)[, 2]
#   ) |>
#   st_drop_geometry() |>
#   drop_na()

# --- 2. Train/test split (spatial) -------------------------------------------

message("=== Spatial train/test split ===")

create_spatial_split <- function(model_df, test_fraction = 0.2, seed = 42) {
  #' Spatial block cross-validation: split by station location to avoid
  
  #' spatial leakage (nearby stations in both train and test).
  #' Uses k-means on coordinates to create spatial blocks.
  
  set.seed(seed)
  
  coords <- model_df |> select(lon, lat) |> distinct()
  n_blocks <- round(1 / test_fraction)
  
  # Cluster stations into spatial blocks
  km <- kmeans(coords, centers = n_blocks)
  coords$block <- km$cluster
  
  # Assign one block as test
  test_block <- sample(unique(coords$block), 1)
  
  model_df <- model_df |>
    left_join(coords |> select(lon, lat, block), by = c("lon", "lat"))
  
  list(
    train = model_df |> filter(block != test_block),
    test  = model_df |> filter(block == test_block),
    block_assignment = coords
  )
}

# split <- create_spatial_split(model_df)
# train_df <- split$train
# test_df  <- split$test

# --- 3. Model A: Standard Random Forest -------------------------------------

message("=== Model A: Random Forest ===")

train_rf <- function(train_df, target = "ground_no2",
                     features = NULL, ntree = 500) {
  #' Train a standard Random Forest using ranger.
  
  if (is.null(features)) {
    features <- setdiff(names(train_df),
                        c(target, "cell_id", "year_month", "block"))
  }
  
  formula <- as.formula(paste(target, "~", paste(features, collapse = " + ")))
  
  model <- ranger(
    formula    = formula,
    data       = train_df,
    num.trees  = ntree,
    importance = "impurity",
    seed       = 42
  )
  
  message(glue("  RF OOB R² = {round(model$r.squared, 4)}"))
  message(glue("  RF OOB RMSE = {round(sqrt(model$prediction.error), 4)}"))
  
  return(model)
}

evaluate_model <- function(model, test_df, target = "ground_no2") {
  #' Evaluate model on test set. Returns metrics data.frame.
  
  preds <- predict(model, data = test_df)$predictions
  actual <- test_df[[target]]
  
  metrics <- tibble(
    R2   = cor(preds, actual)^2,
    RMSE = sqrt(mean((preds - actual)^2)),
    MAE  = mean(abs(preds - actual)),
    MAPE = mean(abs((preds - actual) / actual)) * 100
  )
  
  message(glue("  Test R² = {round(metrics$R2, 4)}, ",
               "RMSE = {round(metrics$RMSE, 4)}, ",
               "MAE = {round(metrics$MAE, 4)}"))
  
  return(list(predictions = preds, metrics = metrics))
}

# rf_model <- train_rf(train_df)
# rf_eval  <- evaluate_model(rf_model, test_df)

# Variable importance plot
plot_importance <- function(model, top_n = 15) {
  imp <- tibble(
    variable   = names(model$variable.importance),
    importance = model$variable.importance
  ) |>
    arrange(desc(importance)) |>
    head(top_n)
  
  ggplot(imp, aes(x = reorder(variable, importance), y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance (Random Forest)",
         x = NULL, y = "Importance (impurity)") +
    theme_minimal()
}

# --- 4. Model B: Geographically Weighted Random Forest -----------------------

message("=== Model B: Geographically Weighted RF ===")

train_gwrf <- function(train_df, test_df, target = "ground_no2",
                       features = NULL) {
  #' Train a Geographically Weighted Random Forest using SpatialML::grf.
  #' Fits local RF models with spatially varying importance.
  
  if (is.null(features)) {
    features <- setdiff(names(train_df),
                        c(target, "cell_id", "year_month", "block", "lon", "lat"))
  }
  
  # grf requires coordinate columns
  train_coords <- train_df |> select(lon, lat) |> as.matrix()
  test_coords  <- test_df |> select(lon, lat) |> as.matrix()
  
  formula <- as.formula(paste(target, "~", paste(features, collapse = " + ")))
  
  # Fit GW-RF
  gwrf <- grf(
    formula  = formula,
    dframe   = as.data.frame(train_df),
    bw       = NULL,           # automatic bandwidth selection
    kernel   = "adaptive",
    coords   = train_coords,
    ntree    = 200,
    mtry     = NULL,           # auto (sqrt(p))
    forests  = TRUE            # store local forests for prediction
  )
  
  message(glue("  GW-RF fitted with bandwidth = {gwrf$Bandwidth}"))
  
  # Predict on test set
  # gwrf_preds <- predict.grf(gwrf, test_df, x.var.name = "lon", y.var.name = "lat")
  
  return(gwrf)
}

# gwrf_model <- train_gwrf(train_df, test_df)

# Local importance map — shows which variables matter WHERE
plot_local_importance <- function(gwrf_model, grid_sf, variable_name) {
  #' Map the local importance of a variable across space.
  #' Shows spatial non-stationarity in predictor effects.
  
  # gwrf_model$Local.Variable.Importance contains spatially varying importance
  # Join to grid for mapping
  message("  Plotting local importance for: ", variable_name)
}

# --- 5. Model C: Artificial Neural Network -----------------------------------

message("=== Model C: Neural Network ===")

train_ann <- function(train_df, test_df, target = "ground_no2",
                      features = NULL, hidden_units = c(16, 8)) {
  #' Train an ANN using caret + nnet (or keras if available).
  #' Inputs are scaled to [0,1].
  
  if (is.null(features)) {
    features <- setdiff(names(train_df),
                        c(target, "cell_id", "year_month", "block"))
  }
  
  formula <- as.formula(paste(target, "~", paste(features, collapse = " + ")))
  
  # Preprocessing: center & scale
  preproc <- preProcess(train_df[features], method = c("center", "scale"))
  train_scaled <- predict(preproc, train_df)
  test_scaled  <- predict(preproc, test_df)
  
  # Train with caret
  ctrl <- trainControl(method = "cv", number = 5)
  
  model <- train(
    formula,
    data      = train_scaled,
    method    = "nnet",
    trControl = ctrl,
    tuneGrid  = expand.grid(size = hidden_units[1], decay = c(0.01, 0.001)),
    linout    = TRUE,        # regression (not classification)
    maxit     = 500,
    trace     = FALSE
  )
  
  message(glue("  ANN best RMSE (CV) = {round(min(model$results$RMSE), 4)}"))
  
  return(list(model = model, preproc = preproc))
}

# ann_result <- train_ann(train_df, test_df)

# --- 6. Model comparison -----------------------------------------------------

compare_models <- function(results_list) {
  #' Compare all models on test set metrics.
  
  bind_rows(results_list, .id = "model") |>
    arrange(desc(R2))
}

# comparison <- compare_models(list(
#   RF    = rf_eval$metrics,
#   GWRF  = gwrf_eval$metrics,
#   ANN   = ann_eval$metrics
# ))
# write_csv(comparison, file.path(CONFIG$dir_output, "model_comparison.csv"))

# --- 7. Predict on full grid (for mapping) -----------------------------------

predict_full_grid <- function(model, grid_sf) {
  #' Apply trained model to all grid cells (not just station cells).
  #' Returns grid_sf with predicted NO2 column.
  
  grid_data <- grid_sf |> st_drop_geometry()
  preds <- predict(model, data = grid_data)$predictions
  grid_sf$predicted_no2 <- preds
  return(grid_sf)
}

# grid_predicted <- predict_full_grid(rf_model, grid)
# st_write(grid_predicted, file.path(CONFIG$dir_output, "predicted_no2.gpkg"))

message("\n=== Supervised ML Complete ===")

# =============================================================================