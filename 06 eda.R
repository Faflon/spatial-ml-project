# =============================================================================
# 06_eda.R — Exploratory Data Analysis & Visualization
# =============================================================================
# Input:  data/processed/analysis_grid.gpkg
#         data/processed/tropomi_grid_values.csv
#         data/processed/stations_with_grid.gpkg
# Output: figures/ (maps, distributions, correlation plots)
# =============================================================================

source("00_setup.R")
source("10_utils.R")

library(sf)
library(terra)
library(tidyverse)
library(tmap)
library(viridis)
library(corrplot)
library(spdep)
library(glue)

# --- 1. Load integrated data -------------------------------------------------

grid <- st_read(file.path(CONFIG$dir_processed, "analysis_grid.gpkg"), quiet = TRUE)
tropomi <- read_csv(file.path(CONFIG$dir_processed, "tropomi_grid_values.csv"),
                    show_col_types = FALSE)
stations <- st_read(file.path(CONFIG$dir_processed, "stations_with_grid.gpkg"),
                    quiet = TRUE)

# --- 2. TROPOMI spatial distribution maps ------------------------------------

message("=== Mapping TROPOMI NO2 spatial distribution ===")

# Annual mean NO2 per cell
no2_annual <- tropomi |>
  filter(pollutant == "NO2") |>
  group_by(cell_id) |>
  summarise(no2_mean = mean(value, na.rm = TRUE), .groups = "drop")

grid_no2 <- grid |> left_join(no2_annual, by = "cell_id")

# Static map
p_no2_map <- ggplot(grid_no2) +
  geom_sf(aes(fill = no2_mean), color = NA) +
  scale_fill_viridis(option = "inferno", name = expression(NO[2]~"(mol/m²)")) +
  labs(title = "Annual Mean Tropospheric NO₂ — Poland 2024",
       subtitle = "Source: Sentinel-5P TROPOMI") +
  theme_minimal()

ggsave(file.path(CONFIG$dir_figures, "no2_annual_map.png"), p_no2_map,
       width = 10, height = 12, dpi = 300)

# Interactive map (for exploration)
# tmap_mode("view")
# tm_shape(grid_no2) + tm_fill("no2_mean", palette = "inferno", alpha = 0.8)

# --- 3. Monthly time series --------------------------------------------------

message("=== Monthly temporal patterns ===")

monthly_national <- tropomi |>
  group_by(year_month, pollutant) |>
  summarise(
    mean_val   = mean(value, na.rm = TRUE),
    median_val = median(value, na.rm = TRUE),
    sd_val     = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

p_monthly <- ggplot(monthly_national, aes(x = year_month, y = mean_val)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
              alpha = 0.2, fill = "steelblue") +
  facet_wrap(~pollutant, scales = "free_y", ncol = 1) +
  labs(title = "Monthly Mean Pollutant Concentrations — Poland 2024",
       x = "Month", y = "Column density (mol/m²)") +
  theme_minimal()

ggsave(file.path(CONFIG$dir_figures, "monthly_timeseries.png"), p_monthly,
       width = 10, height = 8, dpi = 300)

# --- 4. GIOŚ station data distribution --------------------------------------

message("=== Ground station data distributions ===")

# Histogram of ground-level NO2 values
# (placeholder — adapt after loading actual measurement data)
# p_station_hist <- ggplot(station_monthly, aes(x = mean_value)) +
#   geom_histogram(bins = 50, fill = "coral", alpha = 0.7) +
#   facet_wrap(~param_code, scales = "free") +
#   labs(title = "Distribution of Monthly Mean Pollutant Concentrations",
#        x = "Concentration (µg/m³)", y = "Count") +
#   theme_minimal()

# --- 5. Satellite vs. ground correlation ------------------------------------

message("=== Satellite–Ground correlation ===")

# Join TROPOMI cell values to GIOŚ station locations
# (key validation: how well does satellite NO2 predict ground NO2?)

# station_vs_satellite <- station_monthly |>
#   filter(param_code == "NO2") |>
#   left_join(tropomi |> filter(pollutant == "NO2"),
#             by = c("cell_id", "year_month")) |>
#   rename(ground_no2 = mean_value, satellite_no2 = value)
#
# p_scatter <- ggplot(station_vs_satellite, aes(x = satellite_no2, y = ground_no2)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(method = "lm", color = "red") +
#   labs(title = "Satellite vs. Ground-Level NO₂",
#        x = "TROPOMI tropospheric NO₂ (mol/m²)",
#        y = "Ground NO₂ (µg/m³)") +
#   annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
#            label = paste0("R² = ", round(cor(...)^2, 3))) +
#   theme_minimal()

# --- 6. Feature correlations -------------------------------------------------

message("=== Feature correlation matrix ===")

# Numeric columns from the grid (OSM features, GUS, TROPOMI)
# grid_numeric <- grid |>
#   st_drop_geometry() |>
#   select(where(is.numeric), -cell_id)
#
# cor_matrix <- cor(grid_numeric, use = "pairwise.complete.obs")
# png(file.path(CONFIG$dir_figures, "correlation_matrix.png"),
#     width = 800, height = 800)
# corrplot(cor_matrix, method = "color", type = "lower",
#          tl.cex = 0.7, addCoef.col = "black", number.cex = 0.6)
# dev.off()

# --- 7. Spatial autocorrelation (Global Moran's I) --------------------------

message("=== Spatial autocorrelation ===")

# Test whether NO2 values are spatially clustered (expected: yes)
# grid_no2_clean <- grid_no2 |> filter(!is.na(no2_mean))
# coords <- st_coordinates(st_centroid(grid_no2_clean))
# nb <- knn2nb(knearneigh(coords, k = 8))
# lw <- nb2listw(nb, style = "W")
# moran_result <- moran.test(grid_no2_clean$no2_mean, lw)
# message(glue("Moran's I = {round(moran_result$estimate[1], 4)}, ",
#              "p = {format.pval(moran_result$p.value)}"))

# --- 8. Summary statistics table ---------------------------------------------

message("=== Summary statistics ===")

# summary_table <- grid |>
#   st_drop_geometry() |>
#   select(-cell_id) |>
#   pivot_longer(everything(), names_to = "variable", values_to = "value") |>
#   group_by(variable) |>
#   summarise(
#     n    = sum(!is.na(value)),
#     mean = mean(value, na.rm = TRUE),
#     sd   = sd(value, na.rm = TRUE),
#     min  = min(value, na.rm = TRUE),
#     q25  = quantile(value, 0.25, na.rm = TRUE),
#     median = median(value, na.rm = TRUE),
#     q75  = quantile(value, 0.75, na.rm = TRUE),
#     max  = max(value, na.rm = TRUE)
#   )
# write_csv(summary_table, file.path(CONFIG$dir_output, "summary_statistics.csv"))

message("\n=== EDA Complete ===")
message(glue("Figures saved to: {CONFIG$dir_figures}"))

# =============================================================================