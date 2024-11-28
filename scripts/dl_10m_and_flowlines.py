from pathlib import Path
import json

from pygeohydro import WBD
from pynhd import NHD
import py3dep
import rasterio


# Create base directory for data
def ensure_dir(directory):
    Path(directory).mkdir(parents=True, exist_ok=True)

def huc_boundary(huc12, layer='huc12'):
    wbd = WBD(layer=layer)
    boundary = wbd.byids(layer, huc12)
    geom = boundary.loc[0, 'geometry']
    return geom


# load huc 12s:
with open("../input_data/huc12s.json") as f:
    huc12s = json.load(f)

for huc12 in huc12s:
    # Create directory structure for outputs for that huc
    print(huc12)
    base_dir = f"../data/"
    ensure_dir(base_dir)

    boundary = huc_boundary(huc12)

    print('downloading dem')
    dem_10m = py3dep.get_dem(boundary, resolution=10)
    # ensure no ocean values by filtering out values below 0
    dem_10m = dem_10m.where(dem_10m > 0)

    dem_10m_path = f"{base_dir}/{huc12}-dem_10m.tif"
    # convert to 3310
    dem_10m = dem_10m.rio.reproject("EPSG:3310", resampling=rasterio.enums.Resampling.bilinear)
    dem_10m.rio.to_raster(dem_10m_path)


    print('download flowline')
    nhd = NHD('flowline_hr')
    flowlines = nhd.bygeom(boundary)
    if not flowlines.empty:
        flowlines = flowlines.to_crs("EPSG:3310")
        flowlines.to_file(f"{base_dir}/{huc12}-flowlines.shp")
