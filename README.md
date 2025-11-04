# Chicago Crash Network Analysis

Spatial analysis of traffic crashes in Chicago using Network Kernel Density Estimation (NKDE) to identify high-risk road segments.

## Overview
This project analyzes Chicago traffic crash data (2020-2023) to compute crash density along road networks, incorporating injury severity weights and multiple spatial bandwidths.

## Features
- Network Kernel Density Estimation (NKDE)
- Temporal analysis (annual + combined)
- Severity-weighted crash density
- Multiple spatial scales (213m, 305m, 457m)
- Geospatial processing and visualization

## Data Sources
- Chicago road network
- Crash reports (2020-2023)
- City boundaries

## Methods
- Discontinuous NKDE with quartic kernel
- Severity weighting based on injury types
- Coordinate reprojection (EPSG:26971)
- Geometric operations and spatial joins

## Output
GeoPackage layers with crash density values for:
- Combined 2020-2023 period
- Individual years (2020, 2021, 2022, 2023)
- Both unweighted and severity-weighted densities
