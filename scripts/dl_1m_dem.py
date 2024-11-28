import os
import shutil

import py3dep
import rasterio
import geopandas as gpd


ODIR = "../data/catchments/"
regions = gpd.read_file("../data/all_regions.shp")
regions.crs="EPSG:3310"
flowlines = gpd.read_file("../data/all_flowlines.shp")

if os.path.exists(ODIR):
   shutil.rmtree(ODIR)
os.makedirs(ODIR)


for i,(index,row) in enumerate(regions.iterrows()):
    percent = round(i / len(regions) * 100, 2)
    print(f"processing region {i}, {percent}%, {row['hucID']} {row['cID']}")
    full_catchmentID = row['hucID'] + '_' + str(row['cID'])

    repro = gpd.GeoSeries(row['geometry'], crs=regions.crs).to_crs("EPSG:4326")
    region = repro.item()
    dem = py3dep.get_dem(region, crs="EPSG:4326", resolution=1)
    dem = dem.rio.reproject("EPSG:3310", rescampling=rasterio.enums.Resampling.bilinear)
    flow = flowlines.clip(row['geometry'], keep_geom_type=True, sort=True)

    # save
    flow.to_file(f"{ODIR}/{full_catchmentID}-flowlines.shp")
    dem.rio.to_raster(f"{ODIR}/{full_catchmentID}-1m_dem.tif")
