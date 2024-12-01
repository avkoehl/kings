from pathlib import Path 
import os
import shutil

import pandas as pd
import py3dep
import rasterio
import geopandas as gpd
import rioxarray as rxr
from rioxarray.merge import merge_arrays
from shapely.ops import linemerge
from shapely.geometry import MultiLineString, LineString

def clean_flowline(geom):
    if isinstance(geom, LineString):
        return geom
    merged = linemerge(geom)
    if isinstance(merged, LineString):
        return merged
    else:
        return geom #multilinestring

ODIR = "../data/catchments/"
regions = gpd.read_file("../data/all_regions.shp")
regions.crs = "EPSG:3310"
flowlines = gpd.read_file("../data/all_flowlines.shp")

if not os.path.exists(ODIR):
   os.makedirs(ODIR)

failed_items = []
for i,(index,row) in enumerate(regions.iterrows()):
    percent = round(i / len(regions) * 100, 2)

    full_catchmentID = row['hucID'] + '_' + str(row['cID'])

    if Path(f"{ODIR}/{full_catchmentID}-1m_dem.tif").exists():
        #print(f"already downloaded {full_catchmentID}")
        continue

    print(f"processing region {i}, {percent}%, {row['hucID']} {row['cID']}")

    row = row.copy()
    flow = flowlines.clip(row['geometry'], keep_geom_type=True, sort=True)
    flow['geometry'] = flow['geometry'].apply(clean_flowline)
    # now for cases where that failed
    flow = flow.explode()
    flow = flow.loc[flow.length > 200]

    flow.to_file(f"{ODIR}/{full_catchmentID}-flowlines.shp")

    row['geometry'] = row['geometry'].buffer(20)
    reprojected = gpd.GeoSeries(row['geometry'], crs=regions.crs).to_crs("EPSG:4326")
    geom = reprojected.item()

    try:
        dem = py3dep.get_dem(geom, crs="EPSG:4326", resolution=1)
        dem = dem.rio.reproject("EPSG:3310", 
            resampling=rasterio.enums.Resampling.bilinear)
        dem.rio.to_raster(f"{ODIR}/{full_catchmentID}-1m_dem.tif")
    except Exception as e:
        print(f"failed to download {full_catchmentID} exception: {e}")
        failed_items.append(full_catchmentID)

# for each failed one, split into smaller parts
# download each part and mosaic
regions['fullID'] = regions['hucID'] + '_' + regions['cID'].astype(str)
for item in failed_items:
    print(f"processing {item}")
    row = regions.loc[regions['fullID'] == item].copy()
    row['geometry'] = row['geometry'].buffer(20)
    row = row.to_crs("4326")
    geom = row['geometry']
    bbox = geom.bounds
    bbox['mid_x'] = (bbox['minx'] + bbox['maxx']) / 2
    bbox['mid_y'] = (bbox['miny'] + bbox['maxy']) / 2
    
    # Create four smaller bounding boxes
    split_bboxes = pd.DataFrame([
        {'minx': bbox.iloc[0]['minx'], 'miny': bbox.iloc[0]['miny'], 'maxx': bbox.iloc[0]['mid_x'], 'maxy': bbox.iloc[0]['mid_y']},  # Bottom-left
        {'minx': bbox.iloc[0]['mid_x'], 'miny': bbox.iloc[0]['miny'], 'maxx': bbox.iloc[0]['maxx'], 'maxy': bbox.iloc[0]['mid_y']},  # Bottom-right
        {'minx': bbox.iloc[0]['minx'], 'miny': bbox.iloc[0]['mid_y'], 'maxx': bbox.iloc[0]['mid_x'], 'maxy': bbox.iloc[0]['maxy']},  # Top-left
        {'minx': bbox.iloc[0]['mid_x'], 'miny': bbox.iloc[0]['mid_y'], 'maxx': bbox.iloc[0]['maxx'], 'maxy': bbox.iloc[0]['maxy']}   # Top-right
    ])
    
    dems = []
    for ind,box in split_bboxes.iterrows():
        dem = py3dep.get_dem(tuple(box), crs="EPSG:4326", resolution=1)
        dem = dem.rio.reproject("EPSG:3310", 
            resampling=rasterio.enums.Resampling.bilinear)
        dems.append(dem)


    mosaic = merge_arrays(dems)
    region = regions.loc[regions['fullID'] == item, "geometry"]
    clipped = mosaic.rio.clip(region)
    clipped.rio.to_raster(f"{ODIR}/{item}-1m_dem.tif")
