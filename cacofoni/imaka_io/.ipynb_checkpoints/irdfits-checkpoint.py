# FILE: irdfits.py

# Import packages
import numpy as np
from astropy.io import fits
from cacofoni.imaka_io.initimakadatastruct import initimakadatastruct
from cacofoni.imaka_io.getiparm import get_param_values, get_single_value
from cacofoni.config import CacofoniConfig


def irdfits(ftele, 
            fparam,
            exten=None,
            silent=False,):
    """
    Load an imaka telemetry FITS file into a structured array matching IDL format.

    Parameters
    ----------
    ftele : str
        Path to FITS telemetry file
        
    fparam : str
        Path to parameter text file
        
    exten : list[int] or int, optional
        Controls which FITS extensions to load (default is all 8).
        If an integer is passed, only that extension is loaded.
        If a list of length 8 is passed, each index corresponds to:
            0 = Loop control (always loaded internally)
            1 = Raw WFS camera pixels
            2 = Processed WFS camera pixels
            3 = WFS centroids
            4 = DM voltages
            5 = Average WFS camera pixels
            6 = Average WFS centroids
            7 = Average DM voltages

    Returns
    -------
    cb : list of dict
        Time-series list of telemetry frames. Each entry is a dictionary with keys:
        ['loop', 'wfscam', 'wfs', 'dm'] and sub-keys depending on data type.
    """
    
    # Normalize the extension selector
    if exten is None:
        # Default: load all 8 extensions
        exten = [1] * 8
    elif isinstance(exten, int):
        # Load only the specified extension, if requested
        tmp = [0] * 8
        tmp[exten] = 1
        exten = tmp

    # --- Step 1: Open FITS file and load base telemetry info (loop state) ---
    
    with fits.open(ftele) as hdul:
        
        # Load first extenion to get ntimes and nwfs
        h0 = hdul[0].header   # Primary header
        d0 = hdul[0].data     # Loop state data
        ntimes = h0.get('NAXIS2', d0.shape[-1])    # Number of time steps
        nwfs = int(h0['NWFS'])                     # Number of WFS cameras in use, should match imaka parameter file
        
        if not silent:
            print("Assumptions from telemetery FITS file:") # From configuration file
            print("------------------------------------------------")
            print(f"Number of time steps        = {ntimes}")
            print(f"Number of WFS               = {nwfs}")
            print("------------------------------------------------\n")

        # --- Step 2: Allocate empty data structure for storing telemetry ---
        
        if not silent:
            print("Setting up empty data structure with initimakadatastruct...\n")
        
        cb, nwfs_param = initimakadatastruct(ntimes=ntimes, fparam=fparam, silent=silent)
        
        if nwfs != nwfs_param:
            raise (f"Mismatch in WFS count: FITS file = {nwfs}, Param file = {nwfs_param}")
            
        if not silent:
            print("Done setting up empty data structure with initimakadatastruct.\n")

        # --- Step 3: Extension 0 - Loop state ---
        
        if not silent:
            print("Filling empty structure with telemetry data...\n")
            
            print("Running E0: Loop state\n")
        
        # Always load the first extension (index 0)
        d0 = hdul[0].data.T  # Transpose: from (2, ntimes) → (ntimes, 2)
        loop_state_array = np.int16(d0[0, :])
        loop_cntr_array  = np.uint32(d0[1, :])

        for t in range(ntimes):
            cb[t]['loop']['state'] = loop_state_array[t]
            cb[t]['loop']['cntr'] = loop_cntr_array[t]
         
       # --- Extension 1: Raw WFS camera pixels ---
        if exten[1]:
            if not silent:
                print("Running E1: Raw wfscam data\n")

            d1 = hdul[1].data
            h1 = hdul[1].header
            header_vals = [str(v) for v in h1.values()]

            if not any("No wfscam (rawpixels) data saved" in v for v in header_vals):
                if d1 is not None and d1.ndim == 4:
                    if not silent:
                        print("d1 shape:", d1.shape)
                    for t in range(ntimes):
                        for i in range(nwfs):
                            cb[t]['wfscam'][i]['rawpixels'] = d1[:, :, i, t].T
                else:
                    if not silent:
                        print("E1 skipped: data exists but not in 4D shape — treating as empty.\n")
            else:
                if not silent:
                    print("No wfscam (rawpixels) data saved (E1)\n")
        else:
            if not silent:
                print("Skipping E1: Raw wfscam data\n")


        # --- Extension 2: Processed WFS camera pixels ---
        if exten[2]:
            if not silent:
                print("Running E2: Processed wfscam data\n")

            d2 = hdul[2].data
            h2 = hdul[2].header
            header_vals = [str(v) for v in h2.values()]

            if not any("No wfscam (processed pixels) data saved" in v for v in header_vals):
                if d2 is not None and d2.ndim == 4:
                    if not silent:
                        print("d2 shape:", d2.shape)
                    for t in range(ntimes):
                        for i in range(nwfs):
                            cb[t]['wfscam'][i]['pixels'] = d2[:, :, i, t].T
                else:
                    if not silent:
                        print("E2 skipped: data exists but not in 4D shape — treating as empty.\n")
            else:
                if not silent:
                    print("No wfscam (processed pixels) data saved (E2)\n")
        else:
            if not silent:
                print("Skipping E2: Processed wfscam data\n")



        # --- Step 6: Extension 3 - WFS centroids ---
        
        if exten[3]:
            
            if not silent:
                print("Running E3: WFS centroids\n")
            
            d3 = hdul[3].data  # shape: (ntimes, nwfs, 2*nsub)
            for t in range(ntimes):
                for i in range(nwfs):
                    cb[t]['wfs'][i]['centroids'] = d3[t, i, :]
        else:
            if not silent:
                print("Skipping E3: WFS centroids\n")

        # --- Step 7: Extension 4 - DM voltages and delta voltages ---
        
        if exten[4]:
            
            if not silent:
                print("Running E4: DM voltages\n")
            
            d4 = hdul[4].data
            h4 = hdul[4].header

            header_vals = [str(v) for v in h4.values()]
            if not any("No dm data saved" in v for v in header_vals):
                if d4 is not None and d4.ndim == 3:
                    for t in range(ntimes):
                        cb[t]['dm']['deltav']   = d4[t, 0, :]  # time first
                        cb[t]['dm']['voltages'] = d4[t, 1, :]  # time first
                else:
                    if not silent:
                        print("Warning: d4 exists but is not the expected 3D shape.")
            else:
                if not silent:
                    print("No dm data saved")
        else:
            if not silent:
                print("Skipping E4: DM voltages\n")

        # --- Step 8: Extension 5 - Average WFS camera pixels ---
        
        if exten[5]:
            
            if not silent:
                print("Running E5: Average wfscam\n")
            
            d5 = hdul[5].data  # shape: (nwfs, y, x)
            for i in range(nwfs):
                cb[0]['wfscam'][i]['avepixels'] = d5[i, :, :].T
        else:
            if not silent:
                print("Skipping E5: Average wfscam\n")

        # --- Step 9: Extension 6 - Average WFS centroids ---
        
        if exten[6]:
            
            if not silent:
                print("Running E6: Average WFS centroids\n")
            
            d6 = hdul[6].data  # shape: (nwfs, 2*nsub)
            for i in range(nwfs):
                cb[0]['wfs'][i]['avecentroids'] = d6[i, :]
        else:
            if not silent:
                print("Skipping E6: Average WFS centroids\n")

        # --- Step 10: Extension 7 - Average DM voltages ---
        
        if exten[7]:
            
            if not silent:
                print("Running E7: Average DM voltages")
            
            d7 = hdul[7].data
            h7 = hdul[7].header
            if not any("No average dm data saved" in str(v) for v in h7.values()):
                if d7 is not None:
                    cb[0]['dm']['avevoltages'] = np.array(d7).flatten()
                    # Need to flatten to match the idl style
                else:
                    if not silent:
                        print("Warning: d7 is None")
            else:
                if not silent:
                    print("No average dm data saved")
        else:
            if not silent:
                print("Skipping E7: Average DM voltages\n")

    if not silent:
        print("Finished loading FITS file into Python imakadata structure.\n")
    
    return cb



    
            

