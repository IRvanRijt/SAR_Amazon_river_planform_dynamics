# Author: Ivar van Rijt
# Date: 22-FEB-2026 (latest update)
# Purpose: Majority filter for temporal stability in binary watermasks

#---------------------------------------------------------------------------#
# 1. IMPORTING REQUIRED PACKAGES
#---------------------------------------------------------------------------#

import rasterio
import numpy as np
import os

#---------------------------------------------------------------------------#
# 2. DEFINING INPUT AND OUTPUT PATHS
#---------------------------------------------------------------------------#

# Input: three consecutive quarterly Sentinel-1 water masks these represent t-1, t, and t+1 - change this path to your own input and output path
q1_path = r"C:/Users/.tif"
q2_path = r"C:/Users/.tif"
q3_path = r"C:/Users/.tif"

# Resulting majority filtered water mask
output_path = r"C:/Users/.tif"

#---------------------------------------------------------------------------#
# 3. RUNNING MAJORITY FILTER (CENTERED 3-IMAGE)
#---------------------------------------------------------------------------#

# Open raster files
src1 = rasterio.open(q1_path)
src2 = rasterio.open(q2_path)
src3 = rasterio.open(q3_path)

# Copy spatial metadata from first raster and update binary output, remove nodata and apply LWZ
meta = src1.meta.copy()
meta.update(dtype="uint8", nodata=None, compress="LZW")  # Binary and LWZ compression

# Remove existing output file if it exists

if os.path.exists(output_path):
    os.remove(output_path)

# Block-wise processing of majority filter

with rasterio.open(output_path, "w", **meta) as dst:

    for ji, window in dst.block_windows(1):
        # Read the three quarters
        q1 = src1.read(1, window=window)
        q2 = src2.read(1, window=window)
        q3 = src3.read(1, window=window)

        # Majority filter: 1 if at least 2 of 3 quarters are water
        majority = ((q1 + q2 + q3) >= 2).astype(np.uint8)

        # Write block
        dst.write(majority, 1, window=window)

print("Done! Output saved to:")
print(output_path)

# #---------------------------------------------------------------------------#
# # 4. RUNNING MAJORITY FILTER (2-IMAGE BOUNDARIES)
# #---------------------------------------------------------------------------#

# # This section is used only for boundary quarters (e.g., 2017Q1 or 2025Q4),
# # where a centered 3-image window is not available

# # 4.1 --- Input files ---

# # Two consecutive quarterly Sentinel-1 water masks these represent t, and t+1 OR t and t-1
# img1_path = r"C:/Users/rober/Documents/MGI_WUR/Thesis_Dynamics_Amazon_River/Data/20m_s1_WaterMasks/20m_s1_Merged_Watermasks/s1_2025Q3_Watermask_20m.tif"
# img2_path = r"C:/Users/rober/Documents/MGI_WUR/Thesis_Dynamics_Amazon_River/Data/20m_s1_WaterMasks/20m_s1_Merged_Watermasks/s1_2025Q4_Watermask_20m.tif"

# # 4.2 --- Output path ---

# output_path = r"C:/Users/rober/Documents/MGI_WUR/Thesis_Dynamics_Amazon_River/Data/20m_s1_WaterMasks/20m_s1_Merged_Watermask_Majority_filtered/20m_2025Q4_majority_filtered.tif"

# # 4.3 ---- Open raster files ----

# src1 = rasterio.open(img1_path)
# src2 = rasterio.open(img2_path)

# # Copy spatial metadata from first raster and update binary output, remove nodata and apply LWZ
# meta = src1.meta.copy()
# meta.update(dtype="uint8", nodata=None, compress="LZW")

# # 4.4 --- Remove existing output file if it exists ---

# if os.path.exists(output_path):
#     os.remove(output_path)

# # 4.5 --- Block-wise processing of AND filter ---

# with rasterio.open(output_path, "w", **meta) as dst:

#     for ji, window in dst.block_windows(1):
#         img1 = src1.read(1, window=window)
#         img2 = src2.read(1, window=window)

#         # 2-image AND rule:
#         # Water only if BOTH quarters classify as water
#         result = ((img1 == 1) & (img2 == 1)).astype(np.uint8)

#         dst.write(result, 1, window=window)

# print("Done! Output saved to:")
# print(output_path)