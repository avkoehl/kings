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
import multiprocessing as mp
from functools import partial

import pandas as pd
import rioxarray as rxr
import geopandas as gpd
from valleyx import extract_valleys
from valleyx import ValleyConfig
from valleyx.utils import setup_wbt
from tqdm import tqdm
from loguru import logger

def process_single_pair(args):
    cID, (dem_file, flow_file), input_dir, output_dir, config = args

    logger.enable('valleyx')
    logger.remove()
    log_file = os.path.join(output_dir, f"{cID}.log")
    logger.add(log_file, level="DEBUG")
    
    # Create a process-specific working directory
    working_dir = os.path.join(output_dir, f"{cID}-working_dir-{os.getpid()}/")
    try:
        os.makedirs(working_dir)
        
        # Initialize WhiteboxTools for this process
        wbt = setup_wbt(working_dir, verbose=False, max_procs=1)
        
        # Load data
        dem = rxr.open_rasterio(os.path.join(input_dir, dem_file), masked=True).squeeze()
        flowlines = gpd.read_file(os.path.join(input_dir, flow_file))
        
        # Process
        results = extract_valleys(dem, flowlines, wbt, config)
        
        # Save results
        results['floor'].rio.to_raster(os.path.join(output_dir, f"{cID}-floor.tif"))
        results['flowlines'].to_file(os.path.join(output_dir, f"{cID}-flowlines.shp"))
        results['wallpoints'].to_file(os.path.join(output_dir, f"{cID}-wallpoints.shp"))
        
        success = True
        
    except Exception as e:
        print(f'Failed {cID}: {str(e)}')
        success = False
        
    finally:
        # Clean up
        if os.path.exists(working_dir):
            shutil.rmtree(working_dir)
        logger.remove()
    
    return cID, success

def parallel_process_pairs(files, input_dir, output_dir, config, n_processes=None):
    """
    Parallel processing of DEM and flowline pairs
    
    Parameters:
    -----------
    files : dict
        Dictionary of file pairs {cID: (dem_file, flow_file)}
    input_dir : str
        Input directory path
    output_dir : str
        Output directory path
    config : dict
        Configuration dictionary
    n_processes : int, optional
        Number of parallel processes. Defaults to CPU count - 1
    """
    # Determine number of processes
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)
    
    # Prepare arguments for each process
    process_args = [
        (cID, files[cID], input_dir, output_dir, config)
        for cID in files.keys()
    ]
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize progress tracking
    total_pairs = len(files)
    processed = 0
    failed = []
    
    # Create process pool and run processing
    with mp.Pool(processes=n_processes) as pool:
        # Use tqdm for progress tracking
        for cID, success in tqdm(
            pool.imap_unordered(process_single_pair, process_args),
            total=total_pairs,
            desc="Processing pairs"
        ):
            processed += 1
            if not success:
                failed.append(cID)
    
    # Print summary
    print(f"\nProcessing complete:")
    print(f"Total processed: {processed}/{total_pairs}")
    if failed:
        print(f"Failed pairs: {len(failed)}")
        print("Failed IDs:", failed)
    else:
        print("All pairs processed successfully")

# Example usage:
if __name__ == "__main__":
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
    
    parallel_process_pairs(
        files,
        input_dir=input_dir,
        output_dir=output_dir,
        config=config,
        n_processes=4 
    )
