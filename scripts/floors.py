"""
iterate through each region (hucID, catchmentID)
load flowlines, load dem
run valleyx workflow:
    - align flowlines (10m - 1m)
    - delineate reaches
    - find walls
    - label floors
save floors, wallpoints, flowlines
"""
import os
import shutil

import pandas as pd
import rioxarray as rxr
import geopandas as gpd
from valleyx import extract_valleys
from valleyx import ValleyConfig
from valleyx.utils import setup_wbt

input_dir = "../data/catchments/"
output_dir = "../data/floors/"

# if output_dir doesn't exist make it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

params = {
    'hand_threshold': 10,
    'spacing': 20,
    'minsize': 200,
    'window': 5,
    'sigma': 1.5,
    'line_spacing': 3,
    'line_width': 100,
    'line_max_width': 300,
    'point_spacing': 1,
    'min_hand_jump': 10,
    'ratio': 2.5,
    'min_peak_prominence': 20,
    'min_distance': 5,
    'num_cells': 6,
    'slope_threshold': 13,
    'foundation_slope': 6,
    'buffer': 2,
    'min_points': 10,
    'percentile': .90,
    'max_floor_slope': 15
}
config = ValleyConfig.from_dict(params)

files = {}
for file in os.listdir(input_dir):
    prefix = file.split('-')[0]
    if prefix not in files:
        files[prefix] = [None, None]
    if file.lower().endswith('_dem.tif'):
        files[prefix][0] = file
    if file.lower().endswith('.shp'):
        files[prefix][1] = file

# drop any that already have a .tif file in output_dir
for file in os.listdir(output_dir):
    prefix = file.split('-')[0]
    if file.lower().endswith('.tif'):
        files.pop(prefix, None)

# joblib parallel the rest

count = 0
for cID,(dem_file,flow_file) in files.items():
    print(f"processing {cID}, {count} / {len(files.keys())}")
    # make working dir
    working_dir = os.path.join(output_dir, f"{cID}-working_dir/")
    os.makedirs(working_dir)
    wbt = setup_wbt(working_dir, verbose=False, max_procs=1)

    dem = rxr.open_rasterio(os.path.join(input_dir, dem_file), masked=True).squeeze()
    flowlines = gpd.read_file(os.path.join(input_dir, flow_file))

    try:
        results = extract_valleys(dem, flowlines, wbt, config)
        results['floor'].rio.to_raster(os.path.join(output_dir, f"{cID}-floor.tif"))
        results['flowlines'].to_file(os.path.join(output_dir, f"{cID}-flowlines.shp"))
        results['wallpoints'].to_file(os.path.join(output_dir, f"{cID}-wallpoints.shp"))
        shutil.rmtree(working_dir)
        count +=1
    except:
        print('failed', cID)
        shutil.rmtree(working_dir)
        count +=1
