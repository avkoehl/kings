from pathlib import Path 
import os
import shutil

import py3dep
import rasterio
import geopandas as gpd


ODIR = "../data/catchments/"
regions = gpd.read_file("../data/all_regions.shp")
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


    repro = gpd.GeoSeries(row['geometry'], crs=regions.crs).to_crs("EPSG:4326")
    region = repro.item()
    flow = flowlines.clip(row['geometry'], keep_geom_type=True, sort=True)
    flow.to_file(f"{ODIR}/{full_catchmentID}-flowlines.shp")


    try:
        dem = py3dep.get_dem(region, crs="EPSG:4326", resolution=1)
        dem = dem.rio.reproject("EPSG:3310", 
            rescampling=rasterio.enums.Resampling.bilinear)
        dem.rio.to_raster(f"{ODIR}/{full_catchmentID}-1m_dem.tif")
    except Exception as e:
        print(f"failed to download {full_catchmentID} exception: {e}")
        failed_items.append(full_catchmentID)

if len(failed_items):
   for item in failed_items:
      print(item)
