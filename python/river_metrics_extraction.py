# Author: Ivar van Rijt
# Date: 07_APR_2026 (latest update)
# Purpose: Extracting planform river metrics from quarterly binary water masks

#---------------------------------------------------------------------------#
# 1. IMPORTING REQUIRED PACKAGES
#---------------------------------------------------------------------------#

import os
import re
import glob
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio import warp
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.signal import savgol_filter
from skimage.morphology import skeletonize
from skimage.transform import resize as sk_resize
import networkx as nx
import geopandas as gpd
from shapely.geometry import LineString
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------#
# 2. USER SETTINGS
#---------------------------------------------------------------------------#

# Active subbasin — change this path to your own input and output path
SUBBASIN_DIR = r"C:\Users\"
OUTPUT_BASE  = r"C:\Users\"
TARGET_CRS   = "EPSG:6933"

# Processing parameters
ROLLING_STD_WINDOW      = 4       # quarters for rolling-std band in plots
CENTERLINE_MIN_LENGTH_M = 500     # minimum accepted centerline length (m)
SKELETON_FLOOD_LIMIT    = 5e7     # water pixel count above which coarse skeleton is used
WIDTH_WEIGHT_ALPHA      = 2.0     # edge cost = step x (1/width)^ALPHA
SMOOTHING_WINDOW        = 11      # Savitzky-Golay window (must be odd)
SMOOTHING_POLYORDER     = 3
GPKG_NAME               = "centerlines.gpkg"

# Seed-point snap settings
SEED_SNAP_RADIUS_M    = 5_000    # search radius around each seed (metres)
SEED_WIDTH_PERCENTILE = 50       # width percentile filter within snap radius

# Per-subbasin seed points (WGS-84 lon, lat)
SUBBASIN_SEEDS: dict[str, dict[str, tuple[float, float]]] = {
    "20m s1 abacaxis watermasks": {
        "upstream":   (-58.88445, -7.17635),
        "downstream": (-57.6809,  -3.314),
    },
    "20m s1 amazon floodplain watermasks": {
        "upstream":   (-73.43591, -4.41964),
        "downstream": (-49.58656,  0.128),
    },
    "20m s1 coastal basin north watermasks": {
        "upstream":   (-52.49545,  1.57448),
        "downstream": (-50.15293,  1.33445),
    },
    "20m s1 coastal basin south watermasks": {
        "upstream":   (-44.8317,   -3.54592),
        "downstream": (-44.3665,   -2.4293),
    },
    "20m s1 japura caqueta watermasks": {
        "upstream":   (-72.39661,  -0.61277),
        "downstream": (-65.74981,  -1.79088),
    },
    "20m s1 javari watermasks": {
        "upstream":   (-73.16463,  -6.27476),
        "downstream": (-70.20234,  -4.3377),
    },
    "20m s1 jurua watermasks": {
        "upstream":   (-72.78285,  -8.95089),
        "downstream": (-65.75846,  -2.63143),
    },
    "20m s1 madeira watermasks": {
        "upstream":   (-64.97192, -14.44881),
        "downstream": (-59.05153,  -3.81108),
    },
    "20m s1 manacupuru watermasks": {
        "upstream":   (-61.88866,  -2.84401),
        "downstream": (-60.67828,  -3.26771),
    },
    "20m s1 maranon watermasks": {
        "upstream":   (-78.449647, -5.174448),
        "downstream": (-73.45234,  -4.44286),
    },
    "20m s1 napo watermasks": {
        "upstream":   (-77.79397,  -1.04418),
        "downstream": (-72.91281,  -3.27125),
    },
    "20m s1 negro watermasks": {
        "upstream":   (-69.14524,   2.34714),
        "downstream": (-59.94319,  -3.15788),
    },
    "20m s1 parana madeirinha watermasks": {
        "upstream":   (-61.54187,  -5.28935),
        "downstream": (-59.17294,  -3.58029),
    },
    "20m s1 purus watermasks": {
        "upstream":   (-71.78749, -10.74917),
        "downstream": (-62.0128,   -4.4029),
    },
    "20m s1 putumayo ica watermasks": {
        "upstream":   (-76.59906,   0.60431),
        "downstream": (-67.96167,  -3.16033),
    },
    "20m s1 rio araua watermasks": {
        "upstream":   (-65.09402,  -5.33898),
        "downstream": (-63.13691,  -4.07572),
    },
    "20m s1 rio capim watermasks": {
        "upstream":   (-48.89916,  -3.82854),
        "downstream": (-48.13569,  -1.54738),
    },
    "20m s1 rio curuatingua watermasks": {
        "upstream":   (-54.47292,  -3.30717),
        "downstream": (-54.18215,  -2.63959),
    },
    "20m s1 rio jutai watermasks": {
        "upstream":   (-66.97367,  -2.88006),
        "downstream": (-69.43965,  -5.8107),
    },
    "20m s1 rio nanay watermasks": {
        "upstream":   (-74.5081,   -3.47985),
        "downstream": (-73.2429,   -3.6957),
    },
    "20m s1 rio pacaja watermasks": {
        "upstream":   (-50.28001,  -3.28514),
        "downstream": (-50.89991,  -1.99207),
    },
    "20m s1 rio paru watermasks": {
        "upstream":   (-53.15939,   0.47565),
        "downstream": (-52.1086,   -1.1774),
    },
    "20m s1 rio piorini watermasks": {
        "upstream":   (-63.947,    -2.86918),
        "downstream": (-63.2972,   -3.5179),
    },
    "20m s1 rio uatuma watermasks": {
        "upstream":   (-59.57402,  -0.88815),
        "downstream": (-57.79436,  -2.54538),
    },
    "20m s1 tapajos watermasks": {
        "upstream":   (-54.68759, -14.35608),
        "downstream": (-54.71985,  -2.4095),
    },
    "20m s1 tocantins watermasks": {
        "upstream":   (-48.29385, -13.83188),
        "downstream": (-49.468499, -2.299397),
    },
    "20m s1 trombetas watermasks": {
        "upstream":   (-56.93731,   0.11871),
        "downstream": (-55.88931,  -1.76624),
    },
    "20m s1 ucayali watermasks": {
        "upstream":   (-74.19527,  -9.63775),
        "downstream": (-73.49513,  -4.52191),
    },
    "20m s1 xingu watermasks": {
        "upstream":   (-53.35859, -14.02664),
        "downstream": (-52.21716,  -1.87276),
    },
}

#---------------------------------------------------------------------------#
# 3. HELPER FUNCTIONS
#---------------------------------------------------------------------------#

def extract_quarter(fname: str) -> str | None:
    """Extract YYYYQN label from filename."""
    m = re.search(r"(\d{4})Q([1-4])", fname)
    return f"{m.group(1)}Q{m.group(2)}" if m else None


def _basin_key(path: str) -> str:
    """Normalise folder name for lookup in SUBBASIN_SEEDS."""
    return os.path.basename(path).lower().replace("_", " ").replace("-", " ")


def skeletonize_coarse(mask: np.ndarray, factor: int = 4) -> tuple[np.ndarray, int]:
    """Downsample mask, skeletonize, return small skeleton and scale factor."""
    small      = mask[::factor, ::factor]
    skel_small = skeletonize(small)
    return skel_small.astype(bool), factor


def smooth_centerline(line: LineString) -> LineString:
    """Apply Savitzky-Golay filter to remove skeleton staircase artefacts."""
    coords = np.array(line.coords)
    if len(coords) < SMOOTHING_WINDOW:
        return line
    smooth = savgol_filter(
        coords, window_length=SMOOTHING_WINDOW,
        polyorder=SMOOTHING_POLYORDER, axis=0,
    )
    return LineString(smooth)


def distance_transform_tiled(
    mask: np.ndarray,
    tile_size: int = 4096,
    overlap: int   = 512,
) -> np.ndarray:
    """
    Euclidean distance transform processed in overlapping tiles.
    Raises MemoryError early if the output array would exceed 2 GB.
    """
    H, W = mask.shape
    size_gb = H * W * 4 / 1024**3
    if size_gb > 2.0:
        raise MemoryError(
            f"distance_transform_tiled: output would be {size_gb:.1f} GB "
            f"({H}x{W} float32). Use streaming accumulation instead."
        )
    dist = np.zeros((H, W), dtype=np.float32)
    r0 = 0
    while r0 < H:
        r1 = min(r0 + tile_size, H)
        c0 = 0
        while c0 < W:
            c1  = min(c0 + tile_size, W)
            er0 = max(r0 - overlap, 0);  er1 = min(r1 + overlap, H)
            ec0 = max(c0 - overlap, 0);  ec1 = min(c1 + overlap, W)
            tile_dist = ndi.distance_transform_edt(
                mask[er0:er1, ec0:ec1]
            ).astype(np.float32)
            cr0 = r0 - er0;  cr1 = cr0 + (r1 - r0)
            cc0 = c0 - ec0;  cc1 = cc0 + (c1 - c0)
            dist[r0:r1, c0:c1] = tile_dist[cr0:cr1, cc0:cc1]
            c0 = c1
        r0 = r1
    return dist


def build_width_weighted_graph(
    skel: np.ndarray,
    dist_transform: np.ndarray,
) -> tuple[nx.Graph, dict, np.ndarray]:
    """
    Build skeleton graph where edge cost = step × (1/mean_width)^ALPHA.
    Wide channels are cheap paths, narrow tributaries are expensive.
    """
    coords = np.argwhere(skel)
    idx    = {(int(r), int(c)): i for i, (r, c) in enumerate(coords)}

    G = nx.Graph()
    for (r, c), i in idx.items():
        G.add_node(i, rc=(r, c), width=float(dist_transform[r, c]))

    for (r, c), i in idx.items():
        w_i = float(dist_transform[r, c])
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == dc == 0:
                    continue
                nb = (r + dr, c + dc)
                if nb not in idx:
                    continue
                j      = idx[nb]
                w_j    = float(dist_transform[nb[0], nb[1]])
                mean_w = max((w_i + w_j) / 2.0, 1e-6)
                step   = 1.4142 if (dr != 0 and dc != 0) else 1.0
                cost   = step * (1.0 / mean_w) ** WIDTH_WEIGHT_ALPHA
                if not G.has_edge(i, j):
                    G.add_edge(i, j, weight=cost)

    return G, idx, coords


def _lonlat_to_pixel(lon: float, lat: float, transform) -> tuple[int, int]:
    """Project WGS-84 lon/lat to raster pixel (row, col) in TARGET_CRS."""
    from pyproj import Transformer

    if not (-180 <= lon <= 180):
        warnings.warn(f"Longitude {lon} outside -180..180 — may be swapped.")
    if lon > 0:
        warnings.warn(f"Longitude {lon} is positive — Amazon seeds should be negative.")
    if lat < -25 or lat > 10:
        warnings.warn(f"Latitude {lat} outside expected Amazon range (-25..10).")

    tr       = Transformer.from_crs("EPSG:4326", TARGET_CRS, always_xy=True)
    x, y     = tr.transform(lon, lat)
    col, row = ~transform * (x, y)
    row_i, col_i = int(round(row)), int(round(col))

    print(f"    seed lon={lon:.4f} lat={lat:.4f} "
          f"-> {TARGET_CRS} ({x:.0f}, {y:.0f}) "
          f"-> pixel (row={row_i}, col={col_i})")
    return row_i, col_i


def extract_centerline(
    water_mask      : np.ndarray,
    transform,
    pixel_size      : float,
    seed_upstream   : tuple[float, float],
    seed_downstream : tuple[float, float],
    buffer          : int = 10,
) -> LineString | None:
    """
    Extract a width-weighted Dijkstra centerline between two seed points.
    Returns a smoothed LineString or None if extraction fails.
    Failures are recorded as NaN — no fallback methods are used.

    All heavy arrays (distance transform, skeleton, graph) are kept in
    crop-space to avoid full-image allocations on large basins.
    Basins exceeding 2 GB or 50 M water pixels use a 4x coarse skeleton.
    """
    rows, cols = np.nonzero(water_mask)
    if len(rows) < 10:
        warnings.warn("Water mask too small — quarter recorded as NaN.")
        return None

    # Crop tightly around water pixels
    r0 = max(rows.min() - buffer, 0)
    r1 = min(rows.max() + buffer, water_mask.shape[0])
    c0 = max(cols.min() - buffer, 0)
    c1 = min(cols.max() + buffer, water_mask.shape[1])

    water_crop          = water_mask[r0:r1, c0:c1]
    n_water             = np.count_nonzero(water_crop)
    crop_h, crop_w      = water_crop.shape
    crop_gb             = crop_h * crop_w * 4 / 1024**3
    force_coarse        = (n_water > SKELETON_FLOOD_LIMIT) or (crop_gb > 2.0)

    if force_coarse:
        coarse_factor = 4
        print(f"    large crop ({crop_h}x{crop_w}, {n_water/1e6:.1f}M water px,"
              f" {crop_gb:.1f} GB) — using coarse mode (factor={coarse_factor})")
        water_coarse  = water_crop[::coarse_factor, ::coarse_factor]
        struct_c      = ndi.generate_binary_structure(2, 1)
        dilate_coarse = ndi.binary_dilation(water_coarse, structure=struct_c,
                                            iterations=1)
        skel_small, _ = skeletonize_coarse(dilate_coarse, factor=1)
        dist_small    = distance_transform_tiled(water_coarse)
        skel          = skel_small
        dist_crop_use = dist_small
        coord_scale   = coarse_factor
    else:
        coarse_factor = 1
        coord_scale   = 1
        # 1-pixel dilation using cross structuring element closes narrow pixel gaps
        struct       = ndi.generate_binary_structure(2, 1)
        water_dilate = ndi.binary_dilation(water_crop, structure=struct, iterations=1)
        try:
            dist_crop_use = distance_transform_tiled(water_crop)
        except MemoryError:
            warnings.warn("Distance transform too large — quarter recorded as NaN.")
            return None
        try:
            skel = skeletonize(water_dilate)
        except MemoryError:
            skel_small, _ = skeletonize_coarse(water_dilate, factor=4)
            skel          = skel_small
            dist_crop_use = distance_transform_tiled(water_crop[::4, ::4])
            coord_scale   = 4

    if np.count_nonzero(skel) < 2:
        warnings.warn("Skeleton too small — quarter recorded as NaN.")
        return None

    G, idx, coords = build_width_weighted_graph(skel, dist_crop_use)
    if len(G.nodes) == 0:
        warnings.warn("Empty skeleton graph — quarter recorded as NaN.")
        return None

    all_nodes = list(G.nodes)
    all_rc    = np.array([G.nodes[n]["rc"] for n in all_nodes])

    def snap_in_crop(lon, lat):
        """Snap WGS-84 seed to nearest wide skeleton node in crop/coarse space."""
        full_row, full_col = _lonlat_to_pixel(lon, lat, transform)
        crop_row = (full_row - r0) / coord_scale
        crop_col = (full_col - c0) / coord_scale
        seed_pt  = np.array([[crop_row, crop_col]], dtype=np.float32)
        tree     = cKDTree(all_rc)

        for multiplier in (1, 2, 5):
            radius_px = (SEED_SNAP_RADIUS_M * multiplier) / pixel_size
            within    = np.array(tree.query_ball_point(seed_pt[0], radius_px),
                                 dtype=np.intp)
            if len(within) > 0:
                if multiplier > 1:
                    print(f"    ({lon:.4f}, {lat:.4f}): expanded to "
                          f"{SEED_SNAP_RADIUS_M * multiplier / 1000:.0f} km, "
                          f"found {len(within)} candidates")
                break
        else:
            closest_dist, _ = tree.query(seed_pt[0], k=1)
            closest_km      = float(closest_dist) * pixel_size / 1000
            warnings.warn(
                f"Seed ({lon:.4f}, {lat:.4f}): no nodes found at 5x radius. "
                f"Closest node is {closest_km:.1f} km away. "
                f"Quarter recorded as NaN."
            )
            return None

        widths_within = np.array([G.nodes[all_nodes[i]]["width"] for i in within])
        thresh        = np.percentile(widths_within, SEED_WIDTH_PERCENTILE)
        wide_within   = within[widths_within >= thresh]
        if len(wide_within) == 0:
            wide_within = within

        dists_wide = np.linalg.norm(all_rc[wide_within] - seed_pt, axis=1)
        best       = wide_within[int(np.argmin(dists_wide))]
        snapped_rc = G.nodes[all_nodes[best]]["rc"]
        print(f"    snapped -> crop pixel (row={snapped_rc[0]}, col={snapped_rc[1]}), "
              f"width={G.nodes[all_nodes[best]]['width']:.1f} px")
        return all_nodes[best]

    s = snap_in_crop(*seed_upstream)
    e = snap_in_crop(*seed_downstream)

    if s is None or e is None:
        return None

    if s == e:
        warnings.warn("Seeds snapped to same node — quarter recorded as NaN.")
        return None

    # Bridge disconnected components if the gap is on water; NaN if dry land
    if nx.has_path(G, s, e):
        work_graph = G
    else:
        comp_s  = set(nx.node_connected_component(G, s))
        comp_e  = set(nx.node_connected_component(G, e))
        nodes_s = list(comp_s)
        nodes_e = list(comp_e)
        rc_s    = np.array([G.nodes[n]["rc"] for n in nodes_s])
        rc_e    = np.array([G.nodes[n]["rc"] for n in nodes_e])

        tree_e          = cKDTree(rc_e)
        dists_b, idxs_b = tree_e.query(rc_s, k=1)
        best_s_i  = int(np.argmin(dists_b))
        best_e_i  = int(idxs_b[best_s_i])
        bridge_s  = nodes_s[best_s_i]
        bridge_e  = nodes_e[best_e_i]
        gap_px    = float(dists_b[best_s_i])
        gap_m     = gap_px * pixel_size

        print(f"    seeds in disconnected components — gap {gap_m/1000:.2f} km")

        water_check = (water_crop[::coord_scale, ::coord_scale]
                       if coord_scale > 1 else water_crop)
        rs0, cs0 = G.nodes[bridge_s]["rc"]
        re0, ce0 = G.nodes[bridge_e]["rc"]
        n_steps  = max(int(gap_px * 2), 2)
        rs_arr   = np.clip(np.round(np.linspace(rs0, re0, n_steps)).astype(int),
                           0, water_check.shape[0] - 1)
        cs_arr   = np.clip(np.round(np.linspace(cs0, ce0, n_steps)).astype(int),
                           0, water_check.shape[1] - 1)

        # Gaps <= 3 pixels are rasterisation artefacts — bridge unconditionally
        if gap_px <= 3.0:
            accept_bridge = True
            print(f"    artefact gap ({gap_px:.1f} px = {gap_m:.0f} m) — bridged")
        else:
            accept_bridge = bool(water_check[rs_arr, cs_arr].all())

        if accept_bridge:
            mean_w      = max(
                (dist_crop_use[rs0, cs0] + dist_crop_use[re0, ce0]) / 2.0, 1e-6)
            bridge_cost = gap_px * (1.0 / mean_w) ** WIDTH_WEIGHT_ALPHA
            work_graph  = G.copy()
            work_graph.add_edge(bridge_s, bridge_e, weight=bridge_cost)
            print(f"    water bridge accepted ({gap_m/1000:.2f} km)")
        else:
            warnings.warn(
                f"Gap of {gap_m/1000:.1f} km crosses dry land — "
                f"quarter recorded as NaN."
            )
            return None

    path = nx.shortest_path(work_graph, s, e, weight="weight")
    if len(path) < 2:
        warnings.warn("Path collapsed to single node — quarter recorded as NaN.")
        return None

    # Convert skeleton-space coords back to projected metres
    px = [work_graph.nodes[n]["rc"] for n in path]
    xy = [rasterio.transform.xy(
              transform,
              int(r) * coord_scale + r0,
              int(c) * coord_scale + c0,
              offset="center")
          for r, c in px]
    line = LineString(xy)

    if line.length < CENTERLINE_MIN_LENGTH_M:
        warnings.warn(
            f"Centerline too short ({line.length:.0f} m) — quarter recorded as NaN.")
        return None

    return smooth_centerline(line)


#---------------------------------------------------------------------------#
# 4. METRIC CALCULATION FUNCTIONS
#---------------------------------------------------------------------------#

def mean_width(water: np.ndarray, pixel_size: float) -> float:
    """
    Mean channel width = 2 x mean distance-transform value across water pixels.
    Computed tile-by-tile to avoid full-image float32 allocation on large basins.
    """
    H, W      = water.shape
    tile_size = 4096
    overlap   = 512
    total_sum = 0.0
    total_cnt = 0

    r0 = 0
    while r0 < H:
        r1 = min(r0 + tile_size, H)
        c0 = 0
        while c0 < W:
            c1  = min(c0 + tile_size, W)
            er0 = max(r0 - overlap, 0);  er1 = min(r1 + overlap, H)
            ec0 = max(c0 - overlap, 0);  ec1 = min(c1 + overlap, W)
            tile_dist      = ndi.distance_transform_edt(
                water[er0:er1, ec0:ec1]).astype(np.float32)
            cr0 = r0 - er0;  cr1 = cr0 + (r1 - r0)
            cc0 = c0 - ec0;  cc1 = cc0 + (c1 - c0)
            interior_dist  = tile_dist[cr0:cr1, cc0:cc1]
            interior_water = water[r0:r1, c0:c1]
            vals           = interior_dist[interior_water]
            total_sum     += float(vals.sum())
            total_cnt     += int(vals.size)
            c0 = c1
        r0 = r1

    if total_cnt == 0:
        return np.nan
    return float(2.0 * (total_sum / total_cnt) * pixel_size)


def sinuosity(centerline: LineString | None) -> float:
    """Ratio of centerline arc length to straight-line endpoint distance."""
    if centerline is None or len(centerline.coords) < 2:
        return np.nan
    straight = LineString([centerline.coords[0], centerline.coords[-1]])
    sl = straight.length
    return float(centerline.length / sl) if sl > 1e-6 else np.nan


def migration_rate(cl1: LineString, cl2: LineString) -> float:
    """Symmetric mean nearest-neighbour distance between two centerlines."""
    c1, c2 = np.array(cl1.coords), np.array(cl2.coords)
    return float((cKDTree(c2).query(c1)[0].mean()
                + cKDTree(c1).query(c2)[0].mean()) / 2.0)


#---------------------------------------------------------------------------#
# 5. OUTPUT SETUP
#---------------------------------------------------------------------------#

basin   = os.path.basename(SUBBASIN_DIR)
OUT_DIR = os.path.join(OUTPUT_BASE, f"{basin}_River_Metrics")
os.makedirs(os.path.join(OUT_DIR, "centerlines"), exist_ok=True)
GPKG_PATH = os.path.join(OUT_DIR, "centerlines", GPKG_NAME)

# Match active subbasin folder to seed dictionary entry
_key   = _basin_key(SUBBASIN_DIR)
_seeds = next(
    (v for k, v in SUBBASIN_SEEDS.items()
     if k.lower().replace("_", " ") == _key), None
)
SEED_UP   = _seeds["upstream"]   if _seeds else None
SEED_DOWN = _seeds["downstream"] if _seeds else None

if SEED_UP is None:
    print(
        f"\n  WARNING: No seed points found for '{basin}'.\n"
        f"  All centerline-dependent metrics will be NaN.\n"
        f"  Add coordinates to SUBBASIN_SEEDS in Section 2.\n"
    )
else:
    print(f"\n  Seeds loaded — upstream {SEED_UP}  downstream {SEED_DOWN}\n")

#---------------------------------------------------------------------------#
# 6. MAIN PROCESSING LOOP
#---------------------------------------------------------------------------#

records          = []
centerlines      = {}
prev_water       = None
prev_quarter     = None
pixel_area_m2    = None
shared_transform = None
shared_w = shared_h = None
nan_quarters: list[str] = []

tif_files = sorted(glob.glob(os.path.join(SUBBASIN_DIR, "*.tif")))
if not tif_files:
    raise FileNotFoundError(f"No .tif files in {SUBBASIN_DIR}")

for tif in tif_files:
    quarter = extract_quarter(tif)
    if quarter is None:
        print(f"  Skipping (no quarter label): {os.path.basename(tif)}")
        continue

    print(f"\nProcessing {quarter} ...")

    # Load and reproject water mask to TARGET_CRS
    with rasterio.open(tif) as src:
        water   = src.read(1).astype(bool)
        profile = src.profile
        bounds  = src.bounds

    transform, w, h = warp.calculate_default_transform(
        profile["crs"], TARGET_CRS,
        profile["width"], profile["height"], *bounds,
    )

    if shared_transform is None:
        shared_transform = transform
        shared_w, shared_h = w, h
    elif transform != shared_transform or w != shared_w or h != shared_h:
        warnings.warn(f"{quarter}: raster grid differs from first file.")

    water_r = np.zeros((h, w), dtype=np.uint8)
    warp.reproject(
        water.astype(np.uint8), water_r,
        src_transform=profile["transform"], src_crs=profile["crs"],
        dst_transform=transform, dst_crs=TARGET_CRS,
        resampling=warp.Resampling.nearest,
    )
    water_r = water_r.astype(bool)

    pixel_size = (abs(transform.a) + abs(transform.e)) / 2.0
    if pixel_area_m2 is None:
        pixel_area_m2 = pixel_size ** 2

    # Extract centerline using seed points only — NaN on any failure
    centerline = None
    if SEED_UP is not None and SEED_DOWN is not None:
        centerline = extract_centerline(
            water_r, transform, pixel_size,
            seed_upstream=SEED_UP, seed_downstream=SEED_DOWN,
        )

    if centerline is not None:
        centerlines[quarter] = centerline
        print(f"  centerline OK  length={centerline.length/1000:.2f} km")
        gpd.GeoDataFrame(
            {"geometry": [centerline], "quarter": [quarter]},
            crs=TARGET_CRS,
        ).to_file(GPKG_PATH, layer=quarter, driver="GPKG")
    else:
        nan_quarters.append(quarter)
        print(f"  centerline NaN")

    # Calculate metrics for this quarter
    record = {
        "quarter"            : quarter,
        "mean_width_m"       : mean_width(water_r, pixel_size),
        "sinuosity"          : sinuosity(centerline),
        "centerline_length_m": centerline.length if centerline else np.nan,
    }

    if prev_water is not None:
        record["stable_water_m2"] = (
            np.count_nonzero( prev_water &  water_r) * pixel_area_m2)
        record["water_gain_m2"]   = (
            np.count_nonzero(~prev_water &  water_r) * pixel_area_m2)
        record["water_loss_m2"]   = (
            np.count_nonzero( prev_water & ~water_r) * pixel_area_m2)
        if prev_quarter in centerlines and quarter in centerlines:
            record["migration_rate_m"] = migration_rate(
                centerlines[prev_quarter], centerlines[quarter])
        else:
            record["migration_rate_m"] = np.nan
    else:
        record.update(dict(
            stable_water_m2=np.nan, water_gain_m2=np.nan,
            water_loss_m2=np.nan,   migration_rate_m=np.nan,
        ))

    records.append(record)
    prev_water   = water_r
    prev_quarter = quarter

#---------------------------------------------------------------------------#
# 7. EXPORT RESULTS
#---------------------------------------------------------------------------#

df = pd.DataFrame(records).sort_values("quarter").reset_index(drop=True)

# Full time-series CSV
df.to_csv(os.path.join(OUT_DIR, "river_metrics_timeseries.csv"), index=False)

# Overall means CSV
df.drop(columns=["quarter"], errors="ignore") \
  .mean(numeric_only=True).to_frame("mean").T \
  .to_csv(os.path.join(OUT_DIR, "river_metrics_overall_means.csv"), index=False)

print(f"\nCSVs saved -> {OUT_DIR}")

#---------------------------------------------------------------------------#
# 8. VISUALISATION
#---------------------------------------------------------------------------#

def plot_metric(df, col, ylabel, fname, structural_nan_quarters=None):
    """
    Plot metric time series with linear trend and rolling std band.
    Grey bands = structurally missing (no predecessor quarter).
    Red bands  = extraction failed.
    """
    if structural_nan_quarters is None:
        structural_nan_quarters = set()

    x     = np.arange(len(df))
    y     = df[col].values.astype(float)
    valid = ~np.isnan(y)
    if valid.sum() < 2:
        print(f"  Skipping plot '{fname}' — fewer than 2 valid data points.")
        return

    trend       = np.poly1d(np.polyfit(x[valid], y[valid], 1))(x)
    rolling_std = pd.Series(y).rolling(ROLLING_STD_WINDOW, center=True).std()

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(x, y, "o-", label="Data", linewidth=1.5)
    ax.plot(x, trend, "--", label="Trend", linewidth=1.2, color="tomato")
    ax.fill_between(x, y - rolling_std, y + rolling_std,
                    alpha=0.20, label=f"±1 std (rolling {ROLLING_STD_WINDOW}Q)")

    # Vertical lines at year boundaries
    years = df["quarter"].str[:4]
    for i in range(1, len(years)):
        if years.iloc[i] != years.iloc[i - 1]:
            ax.axvline(i - 0.5, color="lightgrey", linestyle="--", lw=0.8)

    # NaN shading
    for i, (q, is_nan) in enumerate(zip(df["quarter"], np.isnan(y))):
        if is_nan:
            color = "grey" if q in structural_nan_quarters else "red"
            alpha = 0.12   if q in structural_nan_quarters else 0.08
            ax.axvspan(i - 0.4, i + 0.4, color=color, alpha=alpha)

    ax.set_xticks(x)
    ax.set_xticklabels(df["quarter"], rotation=45, ha="right")
    ax.set_xlabel("Quarter")
    ax.set_ylabel(ylabel)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=300)
    plt.close(fig)


# First quarter is structurally NaN for difference metrics — shade grey
first_quarter = {df["quarter"].iloc[0]} if len(df) > 0 else set()

for col, label, name, structural in [
    ("mean_width_m",       "Mean Width (m)",         "mean_width.png",      set()),
    ("sinuosity",          "Sinuosity (-)",           "sinuosity.png",       set()),
    ("migration_rate_m",   "Migration Rate (m/Q)",    "migration_rate.png",  first_quarter),
    ("stable_water_m2",    "Stable Water Area (m²)",  "stable_water.png",    first_quarter),
    ("water_gain_m2",      "Water Gain (m²)",         "water_gain.png",      first_quarter),
    ("water_loss_m2",      "Water Loss (m²)",         "water_loss.png",      first_quarter),
]:
    plot_metric(df, col, label, name, structural_nan_quarters=structural)

#---------------------------------------------------------------------------#
# 9. RUN SUMMARY
#---------------------------------------------------------------------------#

total_q = len(records)
n_ok    = len(centerlines)
n_nan   = len(nan_quarters)
nan_pct = 100 * n_nan / total_q if total_q > 0 else 0

print("\n" + "=" * 54)
print("  RUN SUMMARY")
print("=" * 54)
print(f"  Subbasin : {basin}")
print(f"  Quarters processed         : {total_q}")
print(f"  Centerlines extracted (OK) : {n_ok}")
print(f"  Centerlines failed (NaN)   : {n_nan}  ({nan_pct:.0f}%)")

if nan_quarters:
    print()
    print("  NaN quarters (check warnings above for reason):")
    for q in nan_quarters:
        print(f"    {q}")
else:
    print()
    print("  All quarters extracted successfully.")

print("=" * 54)
print(f"\nAll results -> {OUT_DIR}")
print("Done.")