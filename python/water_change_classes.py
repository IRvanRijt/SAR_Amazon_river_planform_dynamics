# Author: Ivar van Rijt
# Date: 27-FEB-2026 (latest update)
# Purpose: Extracting quatitative river metrics from binary rivermasks

#---------------------------------------------------------------------------#
# 1. IMPORTING REQUIRED PACKAGES
#---------------------------------------------------------------------------#

import os
import glob
import numpy as np
import rasterio

#---------------------------------------------------------------------------#
# 2. DEFINING INPUT AND OUTPUT PATHS
#---------------------------------------------------------------------------#

# Input: gap-closed, majority-filtered quarterly water masks - change this path to your own input and output path
input_folder = r"C:\Users\"

# Output directory path
output_folder = r"C:\Users\"

# Create output directory if it does not exist
os.makedirs(output_folder, exist_ok=True)

#---------------------------------------------------------------------------#
# 3. LOADING AND SORTING INPUT RASTER FILES
#---------------------------------------------------------------------------#

# Collect all processed quarterly rasters

files = glob.glob(os.path.join(input_folder, "*_majority_filtered_gap_closed.tif"))

# Ensure that files were found
if len(files) == 0:
    raise FileNotFoundError("No raster files found. Check filename pattern.")

def extract_year_quarter(path):
    """
    Extracts YYYYQ# from filename:
    20m_2017Q1_majority_filtered_gap_closed.tif
         ↑
    """
    return os.path.basename(path).split("_")[1]

# Sort files chronologically based on year and quarter
files_sorted = sorted(files, key=lambda x: extract_year_quarter(x))

print(f"\nFound {len(files_sorted)} raster files:\n")
for f in files_sorted:
    print(" -", os.path.basename(f))

#---------------------------------------------------------------------------#
# 4. DEFINING WATER CHANGE CLASSIFICATION FUNCTION
#---------------------------------------------------------------------------#

# Function for pixel comparison between two consecutive images and assigns class
def compute_change_raster(prev_arr, curr_arr):
    """
    Change classes:
        0 = stable non-water (00)
        1 = stable water (11)
        2 = water gain (01)
        3 = water loss (10)
    """
    # Initialize output array
    change = np.zeros(prev_arr.shape, dtype=np.uint8)

    # Boolean masks for water pixels
    prev_water = prev_arr == 1
    curr_water = curr_arr == 1
    
    change[(prev_water & curr_water)] = 1 # Stable water
    change[(~prev_water & curr_water)] = 3 # Water gain
    change[(prev_water & ~curr_water)] = 4 # Water loss

    return change

#---------------------------------------------------------------------------#
# 5. RUNNING CHANGE RASTERS
#---------------------------------------------------------------------------#

# Loop through consecutive raster pairs
for i in range(len(files_sorted) - 1):

    prev_file = files_sorted[i]
    curr_file = files_sorted[i + 1]

    prev_name = extract_year_quarter(prev_file)
    curr_name = extract_year_quarter(curr_file)

    print(f"\nProcessing: {prev_name} → {curr_name}")

    # Open previous and current rasters
    with rasterio.open(prev_file) as src_prev, \
         rasterio.open(curr_file) as src_curr:

        # Safety checks
        if src_prev.width != src_curr.width or src_prev.height != src_curr.height:
            raise ValueError("Raster dimensions do not match!")

        if src_prev.transform != src_curr.transform:
            raise ValueError("Raster geotransforms do not match!")

        meta = src_prev.meta.copy()
        meta.update({
            "count": 1,
            "dtype": rasterio.uint8,
            "compress": "lzw"
        })

        out_name = f"{prev_name}_to_{curr_name}_change.tif"
        out_path = os.path.join(output_folder, out_name)

        # Block-wise processing
        with rasterio.open(out_path, "w", **meta) as dst:

            total_blocks = len(list(src_prev.block_windows(1)))
            block_counter = 0

            for ji, window in src_prev.block_windows(1):

                prev_block = src_prev.read(1, window=window)
                curr_block = src_curr.read(1, window=window)

                change_block = compute_change_raster(prev_block, curr_block)

                dst.write(change_block, 1, window=window)

    print(f"Saved: {out_name}")

print("\n✅ All change layers created successfully!")