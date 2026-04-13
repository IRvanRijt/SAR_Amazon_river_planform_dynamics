// Author: Ivar van Rijt
// Date: 20-FEB-2026 (latest update)
// Purpose: Sentinel-1 VH water mask extraction (Amazon Basin)

//---------------------------------------------------------------------------//
// 1. PREPROCESSING
//---------------------------------------------------------------------------//

// 1.1 --- Define parameters ---

// Quarterly date ranges for Sentinel-1 temporal selection (2017–2025).
// These ranges serve as a reference lookup when selecting a specific
// quarter for the analysis. Only one chosen period is actively used in
// this script, the rest are included for convenience.

// Pick required timestamp:
// 2017 Q1: '2017-01-01', '2017-03-31'
// 2017 Q2: '2017-04-01', '2017-06-30'
// 2017 Q3: '2017-07-01', '2017-09-30'
// 2017 Q4: '2017-10-01', '2017-12-31'

// 2018 Q1: '2018-01-01', '2018-03-31'
// 2018 Q2: '2018-04-01', '2018-06-30'
// 2018 Q3: '2018-07-01', '2018-09-30'
// 2018 Q4: '2018-10-01', '2018-12-31'

// 2019 Q1: '2019-01-01', '2019-03-31'
// 2019 Q2: '2019-04-01', '2019-06-30'
// 2019 Q3: '2019-07-01', '2019-09-30'
// 2019 Q4: '2019-10-01', '2019-12-31'

// 2020 Q1: '2020-01-01', '2020-03-31'
// 2020 Q2: '2020-04-01', '2020-06-30'
// 2020 Q3: '2020-07-01', '2020-09-30'
// 2020 Q4: '2020-10-01', '2020-12-31'

// 2021 Q1: '2021-01-01', '2021-03-31'
// 2021 Q2: '2021-04-01', '2021-06-30'
// 2021 Q3: '2021-07-01', '2021-09-30'
// 2021 Q4: '2021-10-01', '2021-12-31'

// 2022 Q1: '2022-01-01', '2022-03-31'
// 2022 Q2: '2022-04-01', '2022-06-30'
// 2022 Q3: '2022-07-01', '2022-09-30'
// 2022 Q4: '2022-10-01', '2022-12-31'

// 2023 Q1: '2023-01-01', '2023-03-31'
// 2023 Q2: '2023-04-01', '2023-06-30'
// 2023 Q3: '2023-07-01', '2023-09-30'
// 2023 Q4: '2023-10-01', '2023-12-31'

// 2024 Q1: '2024-01-01', '2024-03-31'
// 2024 Q2: '2024-04-01', '2024-06-30'
// 2024 Q3: '2024-07-01', '2024-09-30'
// 2024 Q4: '2024-10-01', '2024-12-31'

// 2025 Q1: '2025-01-01', '2025-03-31'
// 2025 Q2: '2025-04-01', '2025-06-30'
// 2025 Q3: '2025-07-01', '2025-09-30'
// 2025 Q4: '2025-10-01', '2025-12-31'

// Study area geometry
var geometry = AB1_Amazon_Basin.geometry();

// 1.2 --- Load Sentinel-1 ImageCollection ---

// Two collections are created: ASCENDING and DESCENDING orbit passes, 
// Both are later merged to improve temporal density.

var asc = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate('2017-01-01', '2017-03-31')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filterBounds(geometry);

var desc = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate('2017-01-01', '2017-03-31')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filterBounds(geometry);

var sentinel1 = asc.merge(desc);

// 1.3 --- Mask Image borders ---

// Masks out low-incidence-angle edges and scene-start/end artefacts,
// where the incidence angle drops below ~30° or exceeds ~45°.

function maskBorderNew(image) {
  var angle = image.select('angle');
  return image.updateMask(
    createSceneStartEndMask(image, 500)
      .and(angle.gt(30.73993).and(angle.lt(45.03993)))
  );
}

// 1.4 --- Function to handle errors at beginning and end of images ---

// Detects borders by comparing image footprint geometry with orbit slice metadata.
function createSceneStartEndMask(image, bufferMeters) {
  var geometry = image.geometry();
  var coordinates = ee.Array(ee.List(geometry.coordinates().get(0)));
  var size = coordinates.length().get([0]);

  // Determine bounding coordinates (min/max eges)
  var max = coordinates.reduce(ee.Reducer.max(), [0]);
  var min = coordinates.reduce(ee.Reducer.min(), [0]);

  // Extract extreme points for masking
  var xMax = ee.Geometry.Point(coordinates.mask(
    max.slice(1, 0, 1).repeat(0, size).eq(coordinates.slice(1, 0, 1))
  ).slice(0, 0, 1).project([1]).toList());

  var xMin = ee.Geometry.Point(coordinates.mask(
    min.slice(1, 0, 1).repeat(0, size).eq(coordinates.slice(1, 0, 1))
  ).slice(0, 0, 1).project([1]).toList());

  var yMax = ee.Geometry.Point(coordinates.mask(
    max.slice(1, 1).repeat(0, size).eq(coordinates.slice(1, 1))
  ).slice(0, 0, 1).project([1]).toList());

  var yMin = ee.Geometry.Point(coordinates.mask(
    min.slice(1, 1).repeat(0, size).eq(coordinates.slice(1, 1))
  ).slice(0, 0, 1).project([1]).toList());

  // Create geomerty based slice
  var totalSlices = image.getNumber('totalSlices');
  var features = ee.FeatureCollection([
    ee.Feature(ee.Geometry.LineString([xMax, yMax]), {sliceNumber: 1, orbit: 'DESCENDING'}),
    ee.Feature(ee.Geometry.LineString([xMin, yMin]), {sliceNumber: totalSlices, orbit: 'DESCENDING'}),
    ee.Feature(ee.Geometry.LineString([xMax, yMin]), {sliceNumber: 1, orbit: 'ASCENDING'}),
    ee.Feature(ee.Geometry.LineString([xMin, yMax]), {sliceNumber: totalSlices, orbit: 'ASCENDING'})
  ])
  .filter(ee.Filter.and(
    ee.Filter.eq('sliceNumber', image.getNumber('sliceNumber')),
    ee.Filter.eq('orbit', image.getString('orbitProperties_pass'))
  ));

  // Mask expansion and edge removal
  var buffered = features.geometry().buffer(30000);
  var coords = ee.List(geometry.coordinates().get(0));

  // Convert footprint ring to segement FeatureCollection.
  var segments = ee.FeatureCollection(
    coords.zip(coords.slice(1).cat(coords.slice(0, 1)))
      .map(function (segment) {
        return ee.Feature(ee.Geometry.LineString(ee.List(segment)));
      })
  );

  // Retain only footprint segments inside buffered boundary
  var filteredSegments = segments.filter(
    ee.Filter.isContained('.geo', buffered)
  );

  // Mask regions close to the boundary
  return filteredSegments.distance(bufferMeters).not().unmask(1);
}

// 1.5 --- Apply border mask ---
var s1_preprocces_view = sentinel1.map(maskBorderNew);

// 1.6 --- Create mean composite to reduce speckle ---
var composite = s1_preprocces_view.select('VH').mean().clip(geometry);

// 1.7 --- Visualize preprocessing results ---
Map.centerObject(geometry, 5);
Map.addLayer(composite, {bands: 'VH', min: -25, max: -5}, 'VH Composite', false);

//---------------------------------------------------------------------------//
// 2. CLASSIFICATION
//---------------------------------------------------------------------------//

// 2.1 --- Define processing settings ---
var band = 'VH';
var scale = 10;

// 2.2 --- Load image and region ---
var vh = composite.select(band);
var region = geometry;

// 2.3 --- Create binary water mask using fixed threshold ---

// Threshold options:
//   • -24.624 : Histogram classification (fewer false positives)
//   • -17.373 : Otsu classification (more permissive, captures braiding and shoaling waters)

var vhSmoothed = vh.focal_mean({
  radius: 30,
  units: 'meters'
});

var water = vhSmoothed.lt(-17.373).rename('Raw Water Mask');

Map.addLayer(
  water.updateMask(water),
  {palette: ['00BFFF']},
  'Raw Water Mask (preprocessing)',
  false
);

// 2.4 --- Remove layover/shadow with DEM-based elevation filtering ---

// Clip DEM to AOI
var demAOI = DEM.clip(geometry);

// Keep only elevations below 1500 m
var lowElevationMask = demAOI.lt(1500)
  .resample('bilinear')
  .reproject({
    crs: vh.projection(),
    scale: 10
  });

// Apply elevation mask to water classification
var waterLowElevation = water.updateMask(lowElevationMask);

// 2.5 --- Additional DEM-based slope filtering with ESA water override ---

// Load ESA WorldCover 2021 for slope override (2.5) and land-cover filtering (2.6)
var worldCover2021 = Land_Cover
  .filter(ee.Filter.eq('system:index', '2021'))
  .first()
  .select(0)
  .reproject({
    crs: vh.projection(),
    scale: 10
  });

// Class 80 = Permanent water bodies
var esaWaterMaskRaw = worldCover2021.eq(80);

// Apply 5-pixel (50 m) buffer to reduce omission errors
var esaWaterMaskBuffered = esaWaterMaskRaw
  .focal_max({
    radius: 5,
    units: 'pixels'
  });

// Use buffered mask for slope override
var esaWaterMask = esaWaterMaskBuffered;

// Compute slope (degrees) from DEM
var slope = ee.Terrain.slope(demAOI)
  .reproject({
    crs: vh.projection(),
    scale: 10
  });

// Define slope threshold (degrees)
var maxSlopeDeg = 20;

// Standard low-slope mask
var lowSlopeMask = slope.lt(maxSlopeDeg);

// Slope mask, keep if slope < threshold OR if ESA says "water"
var slopeWithESAOverride = lowSlopeMask.or(esaWaterMask);

// Apply combined mask
var waterLowElevationLowSlope = waterLowElevation.updateMask(
  slopeWithESAOverride
);

// 2.6 --- Remove false positives using ESA WorldCover (classes 20, 30, 40, 50 and 60) ---

// Classes to REMOVE:
// 20 = Shrubland
// 30 = Grassland
// 40 = Cropland
// 50 = Built-up
// 60 = Bare / sparse vegetation
var nonWaterLandCover = worldCover2021
  .eq(20)
  .or(worldCover2021.eq(30))
  .or(worldCover2021.eq(40))
  .or(worldCover2021.eq(50))
  .or(worldCover2021.eq(60));

// Keep water pixels NOT overlapping those classes
var waterLandCoverFiltered = waterLowElevationLowSlope.updateMask(
  nonWaterLandCover.not()
);

// 2.7 --- Masking classification noise ---

// to remove isolated misclassified pixels, connect fragmented river segments and retain only large, continuous water patches (rivers)

// Filter to keep only continuous rivers and connect small gaps in the mask
var smoothed = waterLandCoverFiltered
  .focal_max(1)
  .focal_min(1);

// Count connected pixels (8-connected)
var patchSize = smoothed.connectedPixelCount(1024, true);
var minPatchSize = 1024;

// Smooth final river edges 
var largeWater = smoothed
  .updateMask(patchSize.gte(minPatchSize))
  .focal_max(1)
  .focal_min(1)
  .selfMask();

Map.addLayer(largeWater, {palette: ['0000FF']}, 'Filtered River Mask', false);

//---------------------------------------------------------------------------//
// 3. EXPORT
//---------------------------------------------------------------------------//

// 3.1 --- Set desired export scale ---
var exportScale = 20;

// 3.2 --- Compute bounding box of study area ---
var bounds = geometry.bounds();
var coords = ee.List(bounds.coordinates().get(0));

var xs = coords.map(function(pt) { return ee.Number(ee.List(pt).get(0)); });
var ys = coords.map(function(pt) { return ee.Number(ee.List(pt).get(1)); });

var xmin = ee.Number(xs.reduce(ee.Reducer.min()));
var xmax = ee.Number(xs.reduce(ee.Reducer.max()));
var ymin = ee.Number(ys.reduce(ee.Reducer.min()));
var ymax = ee.Number(ys.reduce(ee.Reducer.max()));

// 3.3 --- Create 4 x 4 grid (16 tiles) ---
var nCols = 4;
var nRows = 4;

var xStep = xmax.subtract(xmin).divide(nCols);
var yStep = ymax.subtract(ymin).divide(nRows);

// Reproject final image before export
var largeWater10m = largeWater.reproject({
  crs: 'EPSG:4326',
  scale: 10
});

// 3.4 --- Export function ---
function exportTile(col, row) {

  var x0 = xmin.add(xStep.multiply(col));
  var x1 = x0.add(xStep);

  var y0 = ymin.add(yStep.multiply(row));
  var y1 = y0.add(yStep);

  var tileRect = ee.Geometry.Rectangle([x0, y0, x1, y1]);

  // Intersect with actual basin geometry
  var tileGeom = geometry.intersection(tileRect, 1);

  // Skip empty tiles (important!)
  var area = tileGeom.area(1);
  
  Export.image.toDrive({
    image: largeWater10m.toUint8(),
    description: 's1_2017Q1_WaterMask_tile_' + col + '_' + row,
    folder: 'GEE_Exports_Amazon_Basin_20m_Final_Images_16tiles',
    fileNamePrefix: 's1_2017Q1_WaterMask_tile_' + col + '_' + row,
    region: tileGeom,
    scale: exportScale,
    crs: 'EPSG:4326',
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

// 3.5 --- Loop over 16 tiles ---
for (var col = 0; col < nCols; col++) {
  for (var row = 0; row < nRows; row++) {
    exportTile(col, row);
  }
}