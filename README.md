# Spatial ML Project — Air Pollution Modeling in Poland
Course: Spatial Machine Learning in R (WNE UW, 2026)
Authors: Ondřej Marvan & Adam Jaworski
Supervisor: prof. Katarzyna Kopczewska
Project Overview
Predicting ground-level air pollution concentrations across Poland using satellite-derived atmospheric data (Sentinel-5P TROPOMI) combined with ground station measurements (GIOŚ), demographic features (GUS BDL), and spatial infrastructure data (OpenStreetMap). The project applies supervised and unsupervised spatial machine learning methods including geographically weighted random forests, spatial clustering, and kriging.
Data Sources
SourceWhatResolutionAccessSentinel-5P TROPOMINO₂, SO₂, CO atmospheric columns~5.5km, dailyGoogle Earth Engine / CopernicusGIOŚGround-level NO₂, SO₂, PM10, PM2.5, CO, O₃~250 stations, hourlyREST API / bulk CSVGUS BDLPopulation, vehicles, industry statsGmina levelREST API / bdl packageOpenStreetMapRoads, railways, industrial areas, power plantsVector featuresosmdata packageERA5 (optional)Wind, temperature, boundary layer height0.25°, hourlyecmwfr / GEE
Repository Structure
├── 00_setup.R              # Package installation & project config
├── 01_collect_tropomi.R    # Sentinel-5P data via Google Earth Engine
├── 02_collect_gios.R       # GIOŚ ground station data via API
├── 03_collect_gus.R        # GUS BDL demographic data
├── 04_collect_osm.R        # OpenStreetMap spatial features
├── 05_build_grid.R         # Grid creation & data integration
├── 06_eda.R                # Exploratory data analysis & visualization
├── 07_supervised_ml.R      # RF, GW-RF, ANN models
├── 08_unsupervised_ml.R    # Clustering, agglomeration, LISA
├── 09_kriging.R            # Ordinary, Universal & Regression Kriging
├── 10_utils.R              # Shared helper functions & themes
├── 11_report.Rmd           # Final RPubs report (sources all scripts)
├── data/
│   ├── raw/                # Downloaded source data
│   ├── processed/          # Cleaned & integrated datasets
│   └── output/             # Final analysis outputs
├── figures/                # Maps and plots
└── README.md
Workflow
Run scripts in order (00 → 05 for data, 06 → 09 for analysis). The report 11_report.Rmd sources the analysis scripts and knits to HTML for RPubs publication.
00_setup → 01-04 (parallel data collection) → 05_build_grid → 06_eda →
→ 07_supervised_ml → 08_unsupervised_ml → 09_kriging → 11_report.Rmd → RPubs
Methods (planned)

Supervised ML: Random Forest, Geographically Weighted RF, ANN, CNN
Unsupervised ML: Spatial clustering (pollution hotspots), agglomeration measurement, spatial distribution comparison
Geostatistics: Kriging for surface interpolation
Spatial analysis: Spatial weight matrices, spatial autocorrelation (Moran's I)
