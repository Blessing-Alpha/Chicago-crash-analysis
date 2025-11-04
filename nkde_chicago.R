# ---- CLEAN START ----
rm(list = ls())   # Remove all previous data
graphics.off()    # Close all plots
cat("\014")       # Clear console
gc()             # Free up memory

# ---- Packages ----
pkgs <- c("sf","spNetwork","dplyr","units")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(inst)) install.packages(inst)

library(sf); library(spNetwork); library(dplyr); library(units)

# ---- Paths & layers (matches your successful export) ----
gpkg <- "C:/Users/MSI/OneDrive - Oklahoma A and M System/Data_thesis_Chicago/Processed_Data/exports/chicago_crash_study.gpkg"

roads   <- st_read(gpkg, "roads", quiet = TRUE) |>
  st_cast("MULTILINESTRING", warn = FALSE)
crashes <- st_read(gpkg, "crashes_2020_2023", quiet = TRUE) |>
  st_cast("POINT", warn = FALSE)
city    <- st_read(gpkg, "city", quiet = TRUE)

# ---- Project to meters (EPSG:26971 = NAD83 / Illinois East (m)) ----
roads_m   <- st_transform(roads, 26971)
crashes_m <- st_transform(crashes, 26971)
city_m    <- st_transform(city, 26971)

# Optional: clip roads to the city polygon (speeds up NKDE)
roads_m <- st_intersection(roads_m, st_union(city_m))

# ---- Severity weights (fallback=1 if fields missing) ----
nz <- function(x) ifelse(is.na(x), 0, x)
if(all(c("TotalFatals","AInjuries","BInjuries","CInjuries") %in% names(crashes_m))){
  crashes_m$sev_wt <- 5*nz(crashes_m$TotalFatals) +
    3*nz(crashes_m$AInjuries)   +
    1*nz(crashes_m$BInjuries)   +
    0.5*nz(crashes_m$CInjuries)
} else {
  crashes_m$sev_wt <- 1.0
}

# --- sampling: every 50 m instead of 25 m ---
samp_dist <- 50
samples <- lines_points_along(roads_m, dist = samp_dist)

run_nkde_safe <- function(lines, events, w, samples, bw,
                          grid = c(10,10), max_depth = 12) {
  # first try with 10x10 tiling
  out <- try(nkde(
    lines = lines, events = events, w = w, samples = samples,
    kernel_name = "quartic", bw = bw,
    method = "discontinuous", div = "bw",
    digits = 2, tol = 0.1, max_depth = max_depth,
    grid_shape = grid,           # <-- finer tiling
    verbose = TRUE
  ), silent = TRUE)
  
  if (inherits(out, "try-error")) {
    message("Retrying NKDE with finer grid (12x12) …")
    out <- nkde(
      lines = lines, events = events, w = w, samples = samples,
      kernel_name = "quartic", bw = bw,
      method = "discontinuous", div = "bw",
      digits = 2, tol = 0.1, max_depth = max_depth,
      grid_shape = c(12,12),    # <-- even finer
      verbose = TRUE
    )
  }
  out
}

# bandwidths (~700ft, 1000ft, 1500ft)
bws <- c(`bw213m`=213.36, `bw305m`=304.8, `bw457m`=457.2)

samples_out <- samples
for (nm in names(bws)) {
  bw <- bws[[nm]]
  message("Computing NKDE ", nm, " …")
  
  samples_out[[paste0("I_u_", nm)]] <- run_nkde_safe(
    lines = roads_m, events = crashes_m, w = rep(1, nrow(crashes_m)),
    samples = samples, bw = bw
  )
  
  samples_out[[paste0("I_w_", nm)]] <- run_nkde_safe(
    lines = roads_m, events = crashes_m, w = crashes_m$sev_wt,
    samples = samples, bw = bw
  )
}

# ---- SAVE RESULTS TO GEOPACKAGE FOR MAPPING ----
st_write(samples_out, gpkg, layer = "nkde_points", delete_layer = TRUE, quiet = TRUE)

# Aggregate back to roads for easy mapping
nearest_idx <- sf::st_nearest_feature(samples_out, roads_m)
samples_out$road_id <- nearest_idx
vals <- samples_out |>
  st_drop_geometry() |>
  dplyr::group_by(road_id) |>
  dplyr::summarise(dplyr::across(starts_with("I_"), mean, na.rm = TRUE), .groups = "drop")

roads_nkde <- roads_m
roads_nkde$road_id <- seq_len(nrow(roads_m))
roads_nkde <- dplyr::left_join(roads_nkde, vals, by = "road_id")
st_write(roads_nkde, gpkg, layer = "roads_nkde", delete_layer = TRUE, quiet = TRUE)

message("✓ NKDE finished! Results saved to GeoPackage for mapping.")
message("✓ Layers created: 'nkde_points' (points) and 'roads_nkde' (lines)") 