# Spatial ML Project — Air Pollution Modeling in Poland

**Course:** Spatial Machine Learning in R (WNE UW, 2026)  
**Authors:** Ondřej Marvan & Adam Jaworski  
**Supervisor:** prof. Katarzyna Kopczewska

## Project Overview

Predicting ground-level air pollution concentrations across Poland using satellite-derived atmospheric data (Sentinel-5P TROPOMI) combined with ground station measurements (GIOŚ), demographic features (GUS BDL), and spatial infrastructure data (OpenStreetMap). The project applies supervised and unsupervised spatial machine learning methods including geographically weighted random forests, spatial clustering, and kriging.

## Data Sources

| Source | What | Resolution | Access |
|--------|------|-----------|--------|
| Sentinel-5P TROPOMI | NO₂, SO₂, CO atmospheric columns | ~5.5km, daily | Google Earth Engine / Copernicus |
| GIOŚ | Ground-level NO₂, SO₂, PM10, PM2.5, CO, O₃ | ~250 stations, hourly | REST API / bulk CSV |
| GUS BDL | Population, vehicles, industry stats | Gmina level | REST API |
| OpenStreetMap | Roads, railways, industrial areas, power plants | Vector features | `osmdata` package |

## Repository Structure

```
├── 00_setup.R            # Packages, config, helpers
├── 01_collect_data.R     # ALL data collection (TROPOMI, GIOŚ, GUS, OSM)
├── 02_build_grid.R       # 5km grid + data integration
├── 03_analysis.R         # EDA + supervised ML + unsupervised ML + kriging
├── 04_report.Rmd         # Final RPubs report
├── data/
│   ├── raw/              # Downloaded source data
│   ├── processed/        # Cleaned & integrated datasets
│   └── output/           # Final results (predictions, comparisons)
├── figures/              # Maps and plots
└── README.md
```

## Workflow

```
00_setup  →  01_collect_data  →  02_build_grid  →  03_analysis  →  04_report.Rmd  →  RPubs
```

Each script in `01` and `03` is organized into **independently runnable sections** — uncomment the function calls you need. This lets you work incrementally as data becomes available.

## Methods

- **Supervised ML:** Random Forest, Geographically Weighted RF, ANN
- **Unsupervised ML:** DBSCAN hotspot clustering, LISA (Local Moran's I), PAM spatial typology
- **Geostatistics:** Ordinary Kriging, Universal Kriging, Regression Kriging
- **Spatial analysis:** Spatial weight matrices, Moran's I, variogram analysis
