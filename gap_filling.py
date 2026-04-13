# Author: Ivar van Rijt
# Date: 23-FEB-2026 (latest update)
# Purpose: Filling gaps caused by mis-classification or anthropogenic structures to preserve river connectivity

#---------------------------------------------------------------------------#
# 1. IMPORTING REQUIRED PACKAGES
#---------------------------------------------------------------------------#

import os
import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.ndimage import binary_closing

#---------------------------------------------------------------------------#
# 2. DEFINING INPUT AND OUTPUT PATHS
#---------------------------------------------------------------------------#

# Input: majority filtered binary water mask - change this path to your own input and output path
input_path = r"C:\Users\.tif"

# Output directory path
output_folder = r"C:\Users\"
os.makedirs(output_folder, exist_ok=True)

# Output filename (adds suffix to original name)
output_name = os.path.splitext(os.path.basename(input_path))[0] + "_gap_closed.tif"
output_path = os.path.join(output_folder, output_name)

#---------------------------------------------------------------------------#
# 3. PARAMETER SETTINGS
#---------------------------------------------------------------------------#

structure = np.ones((5, 5), dtype=np.uint8) # 5x5 kernel 
buffer = 2 # Prevents edge artifacts when processing blocks independently
block_size = 2048 # Block size for memory-efficient processing

#---------------------------------------------------------------------------#
# 4. MORPHOLOGICAL GAP FILLING (BINARY CLOSING)
#---------------------------------------------------------------------------#

# Binary closing fills gaps to preserve river morphology
with rasterio.open(input_path) as src:
    profile = src.profile # Copy spatial metadata from input raster
    profile.update(dtype="uint8", compress="LZW", count=1) # Binary and LWZ compression

    # Iterate over raster in blocks defined in parameter settings
    with rasterio.open(output_path, "w", **profile) as dst:
        for row in range(0, src.height, block_size):
            for col in range(0, src.width, block_size):

                # Create buffered processing window to make sure operation does not create artificial seams at block boundaries
                row_start = max(row - buffer, 0)
                col_start = max(col - buffer, 0)

                row_end = min(row + block_size + buffer, src.height)
                col_end = min(col + block_size + buffer, src.width)

                window = Window(
                    col_start,
                    row_start,
                    col_end - col_start,
                    row_end - row_start
                )

                # Read buffered window
                arr = src.read(1, window=window)
                arr = (arr == 1)

                # Apply binary closing
                closed = binary_closing(arr, structure=structure)

                # Remove buffer before writing
                r0 = row - row_start
                c0 = col - col_start
                r1 = r0 + min(block_size, src.height - row)
                c1 = c0 + min(block_size, src.width - col)

                out = closed[r0:r1, c0:c1].astype(np.uint8)

                # Write processed block to output raster
                dst.write(out, 1, window=Window(col, row, out.shape[1], out.shape[0]))

print("Finished! Output saved to:")
print(output_path)
