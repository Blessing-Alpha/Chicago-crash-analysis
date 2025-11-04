# ---- CLEAN START ----
rm(list = ls())   # Remove all previous data
graphics.off()    # Close all plots
cat("\014")       # Clear console
gc()             # Free up memory

# ---- Packages ----
pkgs <- c("sf","spNetwork","dplyr","units", "lubridate")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(inst)) install.packages(inst)

library(sf); library(spNetwork); library(dplyr); library(units); library(lubridate)

# ---- Paths & layers ----
gpkg <- "C:/Users/MSI/OneDrive - Oklahoma A and M System/Data_thesis_Chicago/Processed_Data/exports/chicago_crash_study.gpkg"

roads   <- st_read(gpkg, "roads", quiet = TRUE) |>
  st_cast("MULTILINESTRING", warn = FALSE)
crashes <- st_read(gpkg, "crashes_2020_2023", quiet = TRUE) |>
  st_cast("POINT", warn = FALSE)
city    <- st_read(gpkg, "city", quiet = TRUE)

# ---- Project to meters ----
roads_m   <- st_transform(roads, 26971)
crashes_m <- st_transform(crashes, 26971)
city_m    <- st_transform(city, 26971)

# Extract year from crash date (adjust column name as needed)
if("CrashDT" %in% names(crashes_m)) {
  crashes_m$year <- year(crashes_m$CrashDT)
} else if("CRASH_DATE" %in% names(crashes_m)) {
  crashes_m$year <- year(crashes_m$CRASH_DATE)
} else {
  stop("Could not find date column for year extraction")
}

# Clip roads to city
roads_m <- st_intersection(roads_m, st_union(city_m))

# ---- Severity weights function ----
nz <- function(x) ifelse(is.na(x), 0, x)
calculate_severity <- function(crashes_df) {
  if(all(c("TotalFatals","AInjuries","BInjuries","CInjuries") %in% names(crashes_df))){
    return(5*nz(crashes_df$TotalFatals) +
             3*nz(crashes_df$AInjuries) +
             1*nz(crashes_df$BInjuries) +
             0.5*nz(crashes_df$CInjuries))
  } else {
    return(rep(1.0, nrow(crashes_df)))
  }
}

# ---- Sampling ----
samp_dist <- 50
samples <- lines_points_along(roads_m, dist = samp_dist)

# ---- NKDE function ----
run_nkde_safe <- function(lines, events, w, samples, bw,
                          grid = c(10,10), max_depth = 12) {
  out <- try(nkde(
    lines = lines, events = events, w = w, samples = samples,
    kernel_name = "quartic", bw = bw,
    method = "discontinuous", div = "bw",
    digits = 2, tol = 0.1, max_depth = max_depth,
    grid_shape = grid, verbose = TRUE
  ), silent = TRUE)
  
  if (inherits(out, "try-error")) {
    message("Retrying NKDE with finer grid (12x12) …")
    out <- nkde(
      lines = lines, events = events, w = w, samples = samples,
      kernel_name = "quartic", bw = bw,
      method = "discontinuous", div = "bw",
      digits = 2, tol = 0.1, max_depth = max_depth,
      grid_shape = c(12,12), verbose = TRUE
    )
  }
  out
}

# ---- Bandwidths ----
bws <- c(`bw213m`=213.36, `bw305m`=304.8, `bw457m`=457.2)

# ---- Get unique years ----
years <- sort(unique(crashes_m$year))
message("Years found: ", paste(years, collapse = ", "))

# ---- Run NKDE for each year ----
all_samples <- list()

for(yr in years) {
  message("\n=== Processing year: ", yr, " ===")
  
  # Filter crashes for this year
  crashes_yr <- crashes_m[crashes_m$year == yr, ]
  message("Crashes in ", yr, ": ", nrow(crashes_yr))
  
  if(nrow(crashes_yr) == 0) {
    message("No crashes for year ", yr, ", skipping...")
    next
  }
  
  # Calculate severity weights for this year
  crashes_yr$sev_wt <- calculate_severity(crashes_yr)
  
  # Initialize samples for this year
  samples_yr <- samples
  samples_yr$year <- yr
  
  # Run NKDE for each bandwidth
  for (nm in names(bws)) {
    bw <- bws[[nm]]
    message("Computing NKDE ", nm, " for ", yr, " …")
    
    samples_yr[[paste0("I_u_", nm)]] <- run_nkde_safe(
      lines = roads_m, events = crashes_yr, w = rep(1, nrow(crashes_yr)),
      samples = samples, bw = bw
    )
    
    samples_yr[[paste0("I_w_", nm)]] <- run_nkde_safe(
      lines = roads_m, events = crashes_yr, w = crashes_yr$sev_wt,
      samples = samples, bw = bw
    )
  }
  
  all_samples[[as.character(yr)]] <- samples_yr
}

# ---- Combine all years and save ----
samples_all_years <- do.call(rbind, all_samples)

# Save points with all years
st_write(samples_all_years, gpkg, layer = "nkde_points_by_year", delete_layer = TRUE, quiet = TRUE)

# Also save individual years for easier mapping
for(yr in years) {
  yr_data <- samples_all_years[samples_all_years$year == yr, ]
  if(nrow(yr_data) > 0) {
    st_write(yr_data, gpkg, layer = paste0("nkde_points_", yr), delete_layer = TRUE, quiet = TRUE)
  }
}

# ---- Aggregate to roads by year ----
roads_all_years <- list()

for(yr in years) {
  yr_samples <- samples_all_years[samples_all_years$year == yr, ]
  
  if(nrow(yr_samples) > 0) {
    nearest_idx <- sf::st_nearest_feature(yr_samples, roads_m)
    yr_samples$road_id <- nearest_idx
    
    vals <- yr_samples |>
      st_drop_geometry() |>
      dplyr::group_by(road_id) |>
      dplyr::summarise(dplyr::across(starts_with("I_"), mean, na.rm = TRUE), .groups = "drop")
    
    roads_yr <- roads_m
    roads_yr$road_id <- seq_len(nrow(roads_m))
    roads_yr$year <- yr
    roads_yr <- dplyr::left_join(roads_yr, vals, by = "road_id")
    
    roads_all_years[[as.character(yr)]] <- roads_yr
  }
}

# Combine and save road data
if(length(roads_all_years) > 0) {
  roads_combined <- do.call(rbind, roads_all_years)
  st_write(roads_combined, gpkg, layer = "roads_nkde_by_year", delete_layer = TRUE, quiet = TRUE)
  
  # Save individual years
  for(yr in years) {
    if(as.character(yr) %in% names(roads_all_years)) {
      st_write(roads_all_years[[as.character(yr)]], gpkg, 
               layer = paste0("roads_nkde_", yr), delete_layer = TRUE, quiet = TRUE)
    }
  }
}

message("✓ NKDE finished for years: ", paste(years, collapse = ", "))
message("✓ Layers saved:")
message("  - nkde_points_by_year (all years combined)")
message("  - nkde_points_YYYY (individual years)")
message("  - roads_nkde_by_year (all years combined)")
message("  - roads_nkde_YYYY (individual years)")

# Check what layers are in the GeoPackage
st_layers("C:/Users/MSI/OneDrive - Oklahoma A and M System/Data_thesis_Chicago/Processed_Data/exports/chicago_crash_study.gpkg")