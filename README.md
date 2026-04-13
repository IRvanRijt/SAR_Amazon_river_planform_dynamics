# Assessment of River Planform Dynamics in the Amazon Basin using Sentinel-1 SAR Data (2017-2025)  — Research Code Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## 1. Introduction

This repository contains all code used for data acquisition, preprocessing, analysis, and validation in the MSc thesis research on Amazon sub-basin river planform dynamics. The research was conducted at **Wageningen University & Research (WUR)**, 2026.

The repository is organised into three languages reflecting the different stages of the workflow:

- **JavaScript (Google Earth Engine)** — Sentinel-1 data selection, preprocessing and binary water mask generation
- **Python** — Postprocessing of water masks, water change class generation and planform river metric extraction
- **R** — Accuracy assessment of water change classifications using Olofsson's validation method

All scripts are provided as-is for reproducibility and transparicy.

---

## 2. Repository Structure & File Descriptions

### JavaScript — Google Earth Engine (GEE)

**`GEE_Sentinel1_WaterMask.js`**

Retrieves and processes Sentinel-1 SAR data to produce quarterly binary water masks for the Amazon Basin.

- **Input:** Sentinel-1 GRD imagery, NASA SRTM DEM, ESA WorldCover land cover map, Amazon Basin shapefile
- **Steps:**
  1. Selects and filters Sentinel-1 scenes by spatial extent and temporal window
  2. Applies spatial and temporal filtering to reduce speckle noise
  3. Masks incidence angle noise at scene edges
  4. Generates mean composite images per quarter
  5. Classifies water using an Otsu threshold
  6. Applies terrain artefact masking using the SRTM DEM
  7. Applies land cover artefact masking using ESA WorldCover
  8. Filters isolated pixels using connected pixel analysis
- **Output:** Quarterly binary water masks (water = 1, non-water = 0)

---

### Python

**`majority_filter.py`**

Reduces classification noise by applying a temporal majority filter across three consecutive quarterly water masks.

- **Input:** Quarterly binary water masks (t−1, t, t+1)
- **Steps:** Retains a pixel as water only if at least 2 of the 3 images classify it as water
- **Output:** Majority-filtered binary water masks

---

**`gap_filling.py`**

Improves river channel connectivity by closing small gaps in the majority-filtered water masks.

- **Input:** Majority-filtered binary water masks
- **Steps:** Applies morphological closing using a 5×5 pixel window to fill gaps within the water body
- **Output:** Gap-closed binary water masks

---

**`subbasin_clip.py`**

Clips the full-extent water masks to individual sub-basin boundaries.

- **Input:** Majority-filtered and gap-closed binary water masks, sub-basin boundary shapefiles
- **Steps:** Clips each quarterly water mask to each sub-basin polygon
- **Output:** Per-sub-basin quarterly binary water masks

---

**`water_change_classes.py`**

Generates water change classification maps by comparing consecutive quarterly water masks.

- **Input:** Majority-filtered and gap-closed binary water masks (t and t+1)
- **Steps:** Assigns one of four change classes per pixel based on bitwise comparison between consecutive quarters:

  | Class | t | t+1 | Description |
  |-------|---|-----|-------------|
  | 0 | Non-water | Non-water | Stable non-water |
  | 1 | Water | Water | Stable water |
  | 2 | Non-water | Water | Water gain |
  | 3 | Water | Non-water | Water loss |

- **Output:** Raster maps with four water change classes per consecutive quarter pair

---

**`river_metrics_extraction.py`**

Extracts six planform river metrics per quarter for each sub-basin from the clipped binary water masks.

- **Input:** Sub-basin clipped binary water masks
- **Steps:**
  1. Reprojects water masks to equal-area CRS (EPSG:6933)
  2. Skeletonizes the water mask using medial-axis thinning
  3. Builds a width-weighted skeleton graph (edge cost = step × (1/width)²)
  4. Extracts mainstem centerlines using Dijkstra's algorithm between manually defined upstream and downstream seed points
  5. Smooths centerlines using a Savitzky-Golay filter
  6. Calculates mean channel width, sinuosity, meander migration rate, stable water area, water gain and water loss
- **Output:** Per-sub-basin CSV files with quarterly metric time series and PNG time-series plots

---

### R

**`validation_olofsson.R`**

Performs accuracy assessment of the water change classification maps following the area-weighted validation method of Olofsson et al. (2014).

- **Input:** Water change class rasters, CSV file containing mapped and reference labels for validation sample points
- **Steps:** Calculates overall accuracy (OA), user's accuracy (UA), producer's accuracy (PA) and associated standard errors at a 95% confidence interval
- **Output:** CSV files and PNG figures with validation metrics and confidence intervals

---

## 3. Usage

### Requirements

| Language | Version |
|----------|---------|
| Python | 3.13 |
| R | 4.4.2 |
| Google Earth Engine account | — |

**Python packages** — install via `pip install -r requirements.txt`:

```
numpy
pandas
rasterio
geopandas
scipy
scikit-image
networkx
shapely
matplotlib
pyproj
```

**R packages:**

```r
install.packages(c("tidyverse", "terra", "ggplot2"))
```

---

### Running the code

**Google Earth Engine (JavaScript)**

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
2. Copy the contents of `GEE_Sentinel1_WaterMask.js` into a new script
3. Upload your Amazon Basin shapefile as a GEE asset and update the asset path in **Section 2** of the script
4. The Otsu threshold value can be adjusted in **Section 2** — the default value used during development is provided as a starting point
5. Run the script and export results to Google Drive

**Python scripts**

1. At the top of each script, update the input and output folder paths to match your local directory structure — these are clearly marked in **Section 2** of each file
2. Processing parameters are set to the values used during development and are found in **Section 2** of each script — adjust if needed
3. For `river_metrics_extraction.py`, add upstream and downstream seed point coordinates for each sub-basin in the `SUBBASIN_SEEDS` dictionary before running
4. Run scripts in the following order:
```
   majority_filter.py → gap_filling.py → subbasin_clip.py
   water_change_classes.py
   river_metrics_extraction.py
```

**R script**

1. At the top of `validation_olofsson.R`, update the input paths for the water change class rasters and validation sample CSV
2. Run the full script — outputs are saved to the defined output folder

---

## 4. Citation

If you use this code in your own research, please cite the associated paper:

> **[Author(s)]** (*2026*). *[Full paper title]*. *[Journal name]*, [Volume]([Issue]), [Pages]. [https://doi.org/XXXXXXXX](https://doi.org/XXXXXXXX)

To cite the code repository directly:

> **IRvanRijt** (2026). *Assessment of River Planform Dynamics in the Amazon Basin using Sentinel-1 SAR Data (2017-2025)  — Research Code Repository
* (Version 1.0). Zenodo. [https://doi.org/10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

---

## 5. License

This repository is licensed under the **MIT License**.
