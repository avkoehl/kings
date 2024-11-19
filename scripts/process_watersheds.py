import os
from pathlib import Path

from pygeohydro import WBD
from pynhd import NHD
import rioxarray as rxr
import whitebox
import geopandas as gpd


# Create base directory for data
def ensure_dir(directory):
    Path(directory).mkdir(parents=True, exist_ok=True)

def huc_boundary(hucid, layer='huc12'):
    wbd = WBD(layer=layer)
    boundary = wbd.byids(layer, huc12)
    geom = boundary.loc[0, 'geometry']
    return geom

def process_watershed(huc12):
    # Create directory structure
    base_dir = f"data/{huc12}"
    ensure_dir(base_dir)

    boundary = huc_boundary(huc12)

    dem_10m = py3dep.get_dem(boundary, resolution=10)
    dem_10m_path = f"{base_dir}/dem_10m.tif"
    dem_10m.rio.to_raster(dem_10m_path)
    
    # Run isobasins
    wbt.isobasins(
        dem_10m_path,
        f"{base_dir}/isobasins.tif",
        size=2000  # Adjust size parameter as needed
    )
    
    # Convert isobasins raster to vector
    wbt.raster_to_vector_polygons(
        f"{base_dir}/isobasins.tif",
        f"{base_dir}/isobasins.shp"
    )
    
    # Read isobasins as geodataframe
    isobasins = gpd.read_file(f"{base_dir}/isobasins.shp")
    
    # Process each isobasin
    for idx, isobasin in isobasins.iterrows():
        basin_id = str(idx).zfill(3)
        basin_geom = isobasin.geometry
        
        # Try to get 1m DEM if available
        # check if available
        try:
            dem_1m = py3dep.get_dem(boundary, resolution=10)
            dem_1m = get_3dep_dem(basin_geom, resolution=1)
            dem_1m.rio.to_raster(f"{base_dir}/{basin_id}-dem-1m.tif")
        except Exception as e:
            print(f"Could not get 1m DEM for basin {basin_id}: {e}")
            # Fall back to 10m DEM
            dem_10m.rio.clip([basin_geom]).rio.to_raster(f"{base_dir}/{basin_id}-dem-10m.tif")
        
        # Get NHD flowlines
        nhd = NHD('flowline')
        flowlines = nhd.bygeom(basin_geom)
        if not flowlines.empty:
            flowlines.to_file(f"{base_dir}/{basin_id}-flowlines.shp")

# Initialize WhiteboxTools
wbt = whitebox.WhiteboxTools()
wbt.set_workingdir("../workingdir/")

# Process all HUC12s
huc12s = [
    '180101070402',
    '180101070401',
    '180101070201',
    '180101070204',
    '180101070207'
]

for huc12 in huc12s:
    print(f"Processing HUC12: {huc12}")
    process_watershed(huc12)
    print(f"Completed processing HUC12: {huc12}")
