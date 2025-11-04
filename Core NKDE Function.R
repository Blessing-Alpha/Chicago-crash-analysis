run_nkde_safe <- function(lines, events, w, samples, bw, grid = c(10,10)) {
  # Primary calculation with error handling
  # Includes automatic retry with refined grid
}
# --- Packages (install once if needed) -----------------------------------------
# install.packages(c("sf","spNetwork","dplyr","units"))
suppressPackageStartupMessages({
  library(sf)
  library(spNetwork)
  library(dplyr)
  library(units)
})

# --- Paths & layer names -------------------------------------------------------
gpkg <- "C:/Users/MSI/OneDrive - Oklahoma A and M System/Data_thesis_Chicago/Processed_Data/exports/chicago_crash_study.gpkg"
roads_name   <- "roads"
crashes_name <- "crashes_2020_2023"
city_name    <- "city"

# --- Read layers ( data are StatePlane feet) -------------------------------
message("Reading layers from GPKG…")
roads   <- st_read(gpkg, roads_name, quiet = TRUE)
crashes <- st_read(gpkg, crashes_name, quiet = TRUE)
city    <- st_read(gpkg, city_name, quiet = TRUE)

# Ensure clean geometry types, drop Z/M if present
roads   <- st_zm(roads)   |> st_cast("MULTILINESTRING", warn = FALSE)
crashes <- st_zm(crashes) |> st_cast("POINT", warn = FALSE)
city    <- st_zm(city)

# --- Reproject to meters: EPSG:26971 (NAD83 / Illinois East) -------------------
message("Reprojecting to EPSG:26971 (meters)…")
roads_m   <- st_transform(roads, 26971)
crashes_m <- st_transform(crashes, 26971)
city_m    <- st_transform(city, 26971)

# Optional: clip roads to city to speed up NKDE a bit
# (skip if your 'roads' are already Chicago-only)
suppressWarnings({
  city_m <- st_make_valid(city_m)
})
roads_m <- st_intersection(roads_m, city_m)

# --- Build a severity weight (adjust coefficients to your needs) ---------------
sev_cols <- c("TotalFatals","AInjuries","BInjuries","CInjuries")
have_cols <- sev_cols %in% names(crashes_m)

if (all(have_cols)) {
  nz <- function(x) ifelse(is.na(x), 0, x)
  crashes_m <- crashes_m |>
    mutate(
      sev_wt = 5*nz(.data$TotalFatals) +
        3*nz(.data$AInjuries)   +
        1*nz(.data$BInjuries)   +
        0.5*nz(.data$CInjuries)
    )
  message("✓ Severity weights built from TotalFatals/A/B/CInjuries.")
} else {
  crashes_m$sev_wt <- 1.0
  message("⚠ Severity fields not found; using unweighted (all ones).")
}

# --- Bandwidths (meters) -------------------------------------------------------
bws <- c(`bw213m`=213.36, `bw305m`=304.8, `bw457m`=457.2)

# --- Small helper to run NKDE and attach fields to a roads object --------------
run_nkde_set <- function(roads_sf, events_sf, label_prefix) {
  out <- roads_sf
  for (nm in names(bws)) {
    bw <- bws[[nm]]
    
    message(sprintf("  NKDE %s — bw=%s m (unweighted)…", nm, bw))
    nk_u <- nkde(
      lines       = roads_sf,
      events      = events_sf,
      w           = rep(1, nrow(events_sf)),
      kernel_name = "quartic",
      bw          = bw,
      method      = "discontinuous",
      div         = "bw",
      agg         = 10,       # ↑ increase for speed, ↓ for precision (5–20 is common)
      digits      = 3,
      tol         = 0.1,
      max_depth   = 8,
      grid_shape  = "regular",
      verbose     = FALSE
    )
    out[[paste0(label_prefix, "_", nm)]] <- as.numeric(nk_u)
    
    message(sprintf("  NKDE %s — bw=%s m (severity-weighted)…", nm, bw))
    nk_w <- nkde(
      lines       = roads_sf,
      events      = events_sf,
      w           = events_sf$sev_wt,
      kernel_name = "quartic",
      bw          = bw,
      method      = "discontinuous",
      div         = "bw",
      agg         = 10,
      digits      = 3,
      tol         = 0.1,
      max_depth   = 8,
      grid_shape  = "regular",
      verbose     = FALSE
    )
    out[[paste0(label_prefix, "w_", nm)]] <- as.numeric(nk_w)
  }
  out
}

# --- 1) Combined 2020–2023 NKDE ------------------------------------------------
message("Running NKDE for ALL years (2020–2023)…")
roads_all <- run_nkde_set(roads_m, crashes_m, "NKDE")

# Write combined result
st_write(roads_all, gpkg, layer = "roads_nkde_all", delete_layer = TRUE, quiet = TRUE)
message("✓ Wrote: roads_nkde_all")

# --- 2) Per-year NKDE layers ---------------------------------------------------
# Use the CrashDT field already exported from ArcGIS
if (!("CrashDT" %in% names(crashes_m))) {
  stop("CrashDT field not found in crashes; export included it previously in ArcGIS.")
}
crashes_m$Year <- as.integer(format(crashes_m$CrashDT, "%Y"))

years <- c(2020L, 2021L, 2022L, 2023L)

for (yy in years) {
  message(sprintf("\nRunning NKDE for year %d…", yy))
  ev_y <- crashes_m %>% filter(Year == yy)
  if (nrow(ev_y) < 50) {
    message(sprintf("  (Skipped %d; too few points)", yy))
    next
  }
  rd_y <- run_nkde_set(roads_m, ev_y, paste0("NKDE", yy))
  st_write(rd_y, gpkg, layer = paste0("roads_nkde_", yy), delete_layer = TRUE, quiet = TRUE)
  message(sprintf("✓ Wrote: roads_nkde_%d", yy))
}

message("\nAll NKDE layers are in your GeoPackage. Suggested viz:")
message("  - roads_nkde_all: fields NKDE_bw213m, NKDEw_bw213m, …")
message("  - roads_nkde_2020..2023: NKDE2020_bw213m, NKDE2020w_bw213m, … etc.")

