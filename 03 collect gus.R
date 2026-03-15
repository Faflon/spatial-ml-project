# =============================================================================
# 03_collect_gus.R — GUS Bank Danych Lokalnych (BDL) Demographic Data
# =============================================================================
# Data source: GUS BDL API (https://bdl.stat.gov.pl/api/v1/)
# What:        Population, urbanization, industry, transport at gmina level
# Output:      Gmina-level data.frame with features for the ML model
# =============================================================================

source("00_setup.R")

library(httr2)
library(jsonlite)
library(tidyverse)
library(glue)

# --- API Configuration -------------------------------------------------------

BDL_BASE <- "https://bdl.stat.gov.pl/api/v1"

# Optional: register at https://bdl.stat.gov.pl to get an API key for higher
# rate limits. Without a key, you get 10 requests/15 seconds.
BDL_API_KEY <- Sys.getenv("BDL_API_KEY", unset = NA)

# --- Helper: paginated GET ---------------------------------------------------

bdl_get <- function(endpoint, params = list(), delay = 1.5) {
  url <- paste0(BDL_BASE, endpoint)
  
  req <- request(url) |>
    req_url_query(!!!params, format = "json") |>
    req_headers("Accept" = "application/json")
  
  # Add API key if available
  if (!is.na(BDL_API_KEY)) {
    req <- req |> req_headers("X-ClientId" = BDL_API_KEY)
  }
  
  resp <- tryCatch({
    req |>
      req_retry(max_tries = 3, backoff = ~3) |>
      req_perform()
  }, error = function(e) {
    warning(glue("BDL API error: {e$message}"))
    return(NULL)
  })
  
  Sys.sleep(delay)
  if (is.null(resp)) return(NULL)
  resp_body_json(resp, simplifyVector = TRUE)
}

bdl_get_all_pages <- function(endpoint, params = list(), page_size = 100) {
  #' Fetch all pages from a paginated BDL endpoint.
  
  all_results <- list()
  page <- 0
  
  repeat {
    params_page <- c(params, list(`page-size` = page_size, page = page))
    result <- bdl_get(endpoint, params_page)
    
    if (is.null(result) || length(result$results) == 0) break
    
    all_results <- c(all_results, list(result$results))
    page <- page + 1
    
    total_pages <- ceiling((result$totalRecords %||% 0) / page_size)
    message(glue("  Page {page}/{total_pages}"))
    
    if (page >= total_pages) break
  }
  
  bind_rows(all_results)
}

# --- 1. Search for relevant BDL variables ------------------------------------
# BDL organizes data hierarchically: Subject > Group > Subgroup > Variable
# We need to find variable IDs for the features we want.

message("=== Searching BDL for relevant variables ===\n")

# Key search terms for our project
search_terms <- c(
  "ludność",            # population
  "gęstość zaludnienia",# population density
  "pojazdy samochodowe",# motor vehicles
  "przemysł",           # industry
  "emisja",             # emissions
  "powierzchnia"        # area
)

found_variables <- list()
for (term in search_terms) {
  message(glue("Searching: '{term}'..."))
  result <- bdl_get("/variables/search", params = list(name = term))
  if (!is.null(result) && length(result$results) > 0) {
    found_variables[[term]] <- result$results |>
      as_tibble() |>
      select(any_of(c("id", "n1", "n2", "n3", "subjectId",
                      "measureUnitName", "level")))
  }
}

# Print found variables for manual selection
all_found <- bind_rows(found_variables, .id = "search_term")
message(glue("\nFound {nrow(all_found)} variables. Review and select IDs below."))

# --- 2. Pre-selected variable IDs -------------------------------------------
# These are common BDL variable IDs — VERIFY them via the search above
# or via the BDL web interface at bdl.stat.gov.pl
#
# You need to look up actual variable IDs from BDL for your year of interest.
# Below are PLACEHOLDER IDs — replace with actual ones after searching.

VARIABLE_IDS <- list(
  population_total   = NULL,  # e.g., "72305" — ludność ogółem
  population_density = NULL,  # e.g., "60559" — gęstość zaludnienia (os/km²)
  area_km2           = NULL,  # e.g., "3510"  — powierzchnia ogółem
  registered_cars    = NULL,  # e.g., "75488" — pojazdy samochodowe zarejestrowane
  urbanization       = NULL   # e.g., derived from urban/rural population
)

message("\n!!! IMPORTANT !!!")
message("Replace NULL values in VARIABLE_IDS with actual BDL variable IDs.")
message("Browse variables at: https://bdl.stat.gov.pl/bdl/metadane")
message("Or use the search results printed above.")

# --- 3. Download data for selected variables ---------------------------------

download_bdl_variable <- function(var_id, var_name, year = 2024) {
  #' Download data for one BDL variable at gmina level for a given year.
  #' Level 6 = gminy (municipalities)
  
  if (is.null(var_id)) {
    message(glue("  Skipping {var_name} (no variable ID set)"))
    return(NULL)
  }
  
  message(glue("  Downloading {var_name} (var_id={var_id})..."))
  
  data <- bdl_get_all_pages(
    glue("/data/by-variable/{var_id}"),
    params = list(year = year, `unit-level` = 6)  # level 6 = gmina
  )
  
  if (is.null(data) || nrow(data) == 0) {
    warning(glue("  No data returned for {var_name}"))
    return(NULL)
  }
  
  data |>
    as_tibble() |>
    mutate(variable_name = var_name)
}

# Download all variables
year <- as.integer(format(as.Date(CONFIG$date_start), "%Y"))

gus_data_list <- list()
for (var_name in names(VARIABLE_IDS)) {
  var_id <- VARIABLE_IDS[[var_name]]
  gus_data_list[[var_name]] <- download_bdl_variable(var_id, var_name, year)
}

# --- 4. Alternative: use the `bdl` R package (simpler) ----------------------

message("\n=== Alternative: using `bdl` R package ===")
message("The `bdl` package wraps BDL API with convenient R functions:\n")
message('  library(bdl)')
message('  # Search for variables')
message('  search_variables("ludność")')
message('  search_variables("gęstość")')
message('  # Get data at gmina level (level = 6)')
message('  pop <- get_data_by_variable(varId = "72305", year = 2024, unitLevel = 6)')
message("")

# --- 5. Get gmina boundaries (for spatial join) -----------------------------

message("=== Gmina boundaries ===")
message("Download administrative boundaries from:")
message("  https://www.geoportal.gov.pl/ → Dane do pobrania → PRG")
message("  Or use: https://gadm.org/download_country.html → Poland → Level 2 (gminy)")
message("")
message("In R:")
message('  # From GADM:')
message('  library(geodata)')
message('  poland_gminy <- gadm(country = "POL", level = 2, path = "data/raw/gus")')
message('  # Convert to sf:')
message('  poland_gminy_sf <- st_as_sf(poland_gminy)')

# --- 6. Save -----------------------------------------------------------------

gus_combined <- bind_rows(gus_data_list)
if (nrow(gus_combined) > 0) {
  write_csv(gus_combined, file.path(CONFIG$dir_raw, "gus", "bdl_variables.csv"))
  message(glue("\nSaved {nrow(gus_combined)} rows to data/raw/gus/bdl_variables.csv"))
}

# Print summary of search results for Adam to review
if (nrow(all_found) > 0) {
  write_csv(all_found, file.path(CONFIG$dir_raw, "gus", "bdl_variable_search.csv"))
  message("Variable search results saved to data/raw/gus/bdl_variable_search.csv")
  message("Review this file and update VARIABLE_IDS in the script.")
}

message("\n=== GUS BDL Collection Complete ===")

# =============================================================================