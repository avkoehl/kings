import rioxarray as rxr
import geopandas as gpd
from valleyx import flow_analysis
from valleyx.utils import setup_wbt
# remove coastline flowlines
# align flowlines to 10m dem
# get outlet points
# compute watersheds
# for those boundaries get the 1m dem

working_dir = "./working_dir/"
wbt =  setup_wbt(working_dir, verbose=True, max_procs=-1)
flowlines = gpd.read_file("../data/180101070401-flowlines.shp")
dem = rxr.open_rasterio("../data/180101070401-dem_10m.tif", masked=True).squeeze()

flowlines = keep_only_streams(flowlines)
flowlines, dataset = flow_analysis(dem, flowlines, wbt)

def keep_only_streams(flowlines):
    return flowlines.loc[flowlines['ftype'] == 460]

def flowline_to_network(flowlines):
    
    return network

def strahler_order(network):

    pass

def intersection_points(network):
    # any intersection where neither stream is first order counts
    pass

def catchments(dem, flowlines):
    # remove coast (keep only stream)
    # align flowline
    # to network
    # strahler order
    # intersection points
    # snap
    # watershed
    # return catchments raster, flowpaths raster, flowlines shapefile
