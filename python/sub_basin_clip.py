# Author: Ivar van Rijt
# Date: 27-FEB-2026 (latest update)
# Purpose: Clipping sub-basin shapes out of full basin watermasks

#---------------------------------------------------------------------------#
# 1. IMPORTING REQUIRED PACKAGES
#---------------------------------------------------------------------------#

import os
import glob
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import numpy as np

#---------------------------------------------------------------------------#
# 2. DEFINING INPUT AND OUTPUT PATHS
#---------------------------------------------------------------------------#
# Input: desired raster Geotiff, clipping border shapefile and outpit DIR - change this path to your own input and output path
input_raster_dir = r"C:\Users.tif"
shapefile_path = r"C:\Users\.shp"
output_dir = r"C:\Users\"

# Create output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

#---------------------------------------------------------------------------#
# 3. CLIP WATERMASKS TO SUB-BASIN SHAPEFILE
#---------------------------------------------------------------------------#

# Load shapefile
gdf = gpd.read_file(shapefile_path)

# Ensure geometry is valid
gdf = gdf[gdf.geometry.notnull()]
geometries = gdf.geometry.values

# Get list of GeoTIFFs
tif_files = glob.glob(os.path.join(input_raster_dir, "*.tif"))
print(f"Found {len(tif_files)} GeoTIFF files.")

# Clip each raster
for tif_path in tif_files:
    filename = os.path.basename(tif_path)
    output_path = os.path.join(output_dir, filename)

    with rasterio.open(tif_path) as src:

        # Reproject shapefile if CRS does not match raster
        if gdf.crs != src.crs:
            gdf_proj = gdf.to_crs(src.crs)
            geometries = gdf_proj.geometry.values
        else:
            gdf_proj = gdf
            geometries = gdf.geometry.values

        # Clip raster
        clipped_image, clipped_transform = mask(
            src,
            geometries,
            crop=True,
            nodata=0
        )

        # Convert to uint8
        clipped_image = clipped_image.astype(np.uint8)

        # Update metadata
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": clipped_image.shape[1],
            "width": clipped_image.shape[2],
            "transform": clipped_transform,
            "dtype": rasterio.uint8,
            "compress": "lzw",
            "tiled": True
        })

        # Write clipped raster
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(clipped_image)

    print(f"Clipped and saved: {filename}")

print("All rasters successfully clipped.")