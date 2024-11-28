import os
import shutil

import rioxarray as rxr
import geopandas as gpd
import networkx as nx
from shapely.geometry import Point
from valleyx import flow_analysis
from valleyx.utils import setup_wbt
from valleyx.terrain.subbasins import label_subbasins_pour_points

def catchments(dem, flowlines, wbt):
    """
    dem, xarray
    flowlines, xarray
    wbt, whitebox instance

    split the huc12 watershed into smaller catchments:

    1. align flowlines to flow dir
    2. process each network within flowlines to:
        - compute strahler order
        - find the mainstem path
        - save all points where a stream order 2 or greater stream meets the mainstem
    3. run watershed algorithm on those pour points
    """
    flowlines = keep_only_streams(flowlines)
    flowlines, dataset = flow_analysis(dem, flowlines, wbt)
    pour_points = subcatchment_breaks(flowlines)
    watersheds = label_subbasin_pour_points(dataset['flow_dir'], pour_points, wbt)
    return watersheds



def keep_only_streams(flowlines):
    return flowlines.loc[flowlines['ftype'] == 460]

def flowlines_to_network(flowlines):
    # brittle if coordinates aren't exact to all decimals
    G = nx.DiGraph()
    for index, line in flowlines.items():
        start = line.coords[0]
        end = line.coords[-1]
        G.add_node(start, streamID=index)
        G.add_node(end, streamID=index)
        G.add_edge(start, end, streamID=index)
    
    return G

def subcatchment_breaks(flowlines, min_order=2):
    G = flowlines_to_network(flowlines)
    for u, v in G.edges():
        G.edges[u, v]['mainstem'] = False

    outlets = [node for node in G.nodes() if G.out_degree(node) == 0]

    # Process each subnetwork
    junction_nodes = []
    for outlet in outlets:
        # Extract subgraph
        upstream = nx.ancestors(G, outlet)
        upstream.add(outlet)
        subG = G.subgraph(upstream).copy()

        calculate_strahler(subG, outlet)
        # Copy edge attributes back to main graph
        for u, v in subG.edges():
            G.edges[u,v].update(subG.edges[u, v])

        mainstem_nodes = find_mainstem(subG, outlet)
        mainstem_edges = []
        for i in range(len(mainstem_nodes)-1):
            u = mainstem_nodes[i]
            v = mainstem_nodes[i+1]
            subG.edges[u,v]['mainstem'] = True 
            G.edges[u,v]['mainstem'] = True 

        # junction nodes
        # start at bottom work up
        for node in mainstem_nodes[::-1]:
            in_edges = list(subG.in_edges(node, data=True))
            for u,v,data in in_edges:
                if data['mainstem']:
                    continue

                if data['strahler'] >= min_order:
                    junction_nodes.append(Point(node))

    pour_points = [Point(o) for o in outlets]
    pour_points = pour_points + junction_nodes
    pour_points = gpd.GeoSeries(pour_points, crs=flowlines.crs)

    orders = []
    for u, v in G.edges():
        orders.append({
            'streamID': G.edges[u,v].get('streamID'),
            'mainstem': G.edges[u,v].get('mainstem'),
            'strahler_order': G.edges[u,v].get('strahler')})

    orders = pd.DataFrame.from_records(orders)
    gdf = flowlines.to_frame('geometry')
    gdf = gdf.merge(
            orders[['streamID', 'strahler_order', 'mainstem']],
            left_index=True,
            right_on='streamID'
    )
    gdf = gdf.set_index('streamID')
    return pour_points

   
def calculate_strahler(G, root_node):
    # Assumes G is a directed graph (DiGraph) with flow direction
    
    def _strahler_recursive(node):
        # Get upstream edges
        in_edges = list(G.in_edges(node))
        
        if not in_edges:  # Leaf node/headwater
            return 1
        
        # Get Strahler numbers of upstream edges
        upstream_orders = []
        for u, v in in_edges:
            if 'strahler' not in G.edges[u, v]:
                upstream_order = _strahler_recursive(u)
                G.edges[u, v]['strahler'] = upstream_order
            upstream_orders.append(G.edges[u, v]['strahler'])
        
        # Calculate Strahler number for current segment
        max_order = max(upstream_orders)
        if upstream_orders.count(max_order) > 1:
            strahler = max_order + 1
        else:
            strahler = max_order
            
        # Set order for edge going downstream from current node
        out_edges = list(G.out_edges(node))
        if out_edges:  # if not outlet
            u, v = out_edges[0]  # should only be one downstream edge
            G.edges[u, v]['strahler'] = strahler
            
        return strahler
    
    return _strahler_recursive(root_node)

def find_mainstem(G, outlet_node):
    """
    Find mainstem by following highest order stream upstream from outlet
    
    Parameters:
    G: NetworkX DiGraph with 'strahler' edge attribute
    outlet_node: starting node
    
    Returns:
    list of nodes representing mainstem path from headwater to outlet
    """
    mainstem = [outlet_node]
    current_node = outlet_node
    
    while True:
        # Get all upstream edges
        upstream_edges = list(G.in_edges(current_node))
        
        # If no more upstream edges, we've reached a headwater
        if not upstream_edges:
            break
            
        # Find the edge with highest Strahler order
        max_order = -1
        next_node = None
        
        for u, v in upstream_edges:
            order = G.edges[u, v]['strahler']
            if order > max_order:
                max_order = order
                next_node = u
                
        # Add to mainstem and continue upstream
        mainstem.append(next_node)
        current_node = next_node
        
    return mainstem[::-1]

if __name__ == "__main__":
    working_dir = "./working_dir/"
    dem_file
    flowlines_file = 


    if os.path.exists(working_dir):
        shutil.rmtree(working_dir)
    os.makedirs(working_dir)
    
    wbt =  setup_wbt(working_dir, verbose=True, max_procs=-1)
    flowlines = gpd.read_file("../data/180101070401-flowlines.shp")
    dem = rxr.open_rasterio("../data/180101070401-dem_10m.tif", masked=True).squeeze()
    watersheds = catchments(dem, flowlines, wbt):
