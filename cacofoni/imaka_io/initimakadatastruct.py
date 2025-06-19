# FILE: initimakadatastruct.py

# Import packages
import numpy as np
from cacofoni.imaka_io.getiparm import get_param_values, get_single_value
from cacofoni.config import CacofoniConfig


def initimakadatastruct(fparam, 
                        ntimes,
                        silent=False):
    """
    Initialize the imaka data structure for storing AO telemetry over time.

    Parameters
    ----------
    fparam : str
        Path to the imaka parameter file (i.e., imakaparm.txt).
    
    ntimes : int
        Number of time steps (frames) to initialize the structure for.
        Set by the telemetry FITS file. 

    Returns
    -------
    imakadata : list of dict
        A list of data frames (one per time step), each holding loop, WFS, and DM data.
    """
    
    # Load configuration object for max allowed number of WFS
    # Leftover idl logic, not sure if it is needed anymore
    # Only 1 WFS appears to be used 
    config = CacofoniConfig()
    nwfs_max = config.nwfs_max

    # --- A) Parse system parameters ---
    
    nsub = get_single_value(fparam, "sys_parm.nsub") # Number of subapertures per axis (e.g., 12 for 12x12)
    nact = get_single_value(fparam, "sys_parm.nact") # Number of DM actuators (reports 64 but actually 36)
    nwfs = get_single_value(fparam, "sys_parm.nwfs") # Number of WFS cameras actually in use (not max)

    # Pixel width of each WFS camera, assumes all WFS cameras are same size
    npixx_all = get_param_values(fparam, "wfscam_parm.npixx", which_column=1, cast_type=int)
    
    # Sanity check: confirm that total pixel widths match expected shape
    if sum(npixx_all) != nwfs * npixx_all[0]:
        print("Error in NPIXX") # A warning, not a fatal error
    npixx = npixx_all[0] # Take the value for one WFS camera (they're all assumed equal)

    # --- B) Derived quantities ---
    
    nsub_total = nsub * nsub              # Total number of subapertures (e.g., 144)
    npix_per_sub = npixx / nsub           # Pixels per subaperture (assumed square)
    frame_size = int(nsub * npix_per_sub) # Final image size in pixels per WFS (e.g., 120x120)

    # Print configuration summary
    if not silent:
        print("Assumptions from parameter file:")
        print("------------------------------------------------")
        print(f"Number of subaps across     = {nsub}")
        print(f"Number of drivers (not act) = {nact}")
        print(f"Number of WFS               = {nwfs}")
        print(f"Total number of subaps      = {nsub_total}")
        print(f"Image frame size (pixels)   = {frame_size}x{frame_size}") # A little redundant 
        print(f"Pixels per subap (1D)       = {npix_per_sub}")
        print("------------------------------------------------\n")
    
    # --- C) Build full data structure ---
    
    imakadata = []
    
    for _ in range(ntimes):
        frame = {
            
            # --- Loop control state ---
            "loop": {
                "state": 0,   # 0 = open loop; may be toggled during runtime
                "cntr": 0,    # frame counter
            },
            
            # --- WFS camera image data ---
            # Pre-allocate raw and processed pixel arrays for each WFS camera slot
            "wfscam": [
                {
                    "timestamp": 0, # placeholder for time the image was taken
                    "fieldcount": 0, # counter for sequencing or field ID
                    "tsample": 0.0, # time interval or sampling rate
                    "rawpixels": np.zeros((frame_size, frame_size), dtype=np.uint16),
                    "pixels": np.zeros((frame_size, frame_size), dtype=np.float32),
                    "avepixels": np.zeros((frame_size, frame_size), dtype=np.float32),
                }
                for _ in range(nwfs_max)
            ],
            
            # --- WFS centroid measurements ---
            # Each centroid = (x, y) pair for each subaperture
            "wfs": [
                {
                    "raw_centroids": np.zeros(2 * nsub_total, dtype=np.float32),
                    "centroids": np.zeros(2 * nsub_total, dtype=np.float32),
                    "avecentroids": np.zeros(2 * nsub_total, dtype=np.float32),
                }
                for _ in range(nwfs_max)
            ],
            
            # --- DM actuator voltage commands ---
            "dm": {
                "deltav": np.zeros(nact, dtype=np.float32),      # voltage change at this step
                "voltages": np.zeros(nact, dtype=np.float32),    # absolute voltages
                "avevoltages": np.zeros(nact, dtype=np.float32), # average over time
            },
        }
        
        # Append this frame to the full time series
        imakadata.append(frame)

    return imakadata, nwfs