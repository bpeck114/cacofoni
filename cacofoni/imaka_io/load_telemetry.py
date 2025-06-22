# FILE: load_telemetry.py

import numpy as np
from astropy.io import fits
from cacofoni.config import CacofoniConfig

def load_telemetry(ftele, 
                   silent=False):
    """
    Simple loader for iMaka telemetry FITS file.
    
    Parameters
    ----------
    ftele : str
        Path to the telemetry FITS file.

    silent : bool
        Suppress print statements.

    Returns
    -------
    data : dict
        Dictionary of numpy arrays with keys:
        'loop', 'centroids', 'dm'
    """

    data = {}
    data['loop'] = {}
    data['wfs'] = {}
    data['dm'] = {}
    
    # Actual nact does not match size of telemetry array
    config = CacofoniConfig
    nact = config.num_actuators
    nsub = config.num_subapertures

    with fits.open(ftele) as hdul:

        if not silent: 
            print(f"Loading Extension 0: Loop state")
            print(f"------------------------------------------------")
        
        d0 = hdul[0].data # shape: (27000, 5) = (ntimes, parameter)
        # parameter = (0) counter, (1) state, (2) clocktime, (3) dTime, (4) WFScam time
        
        if not silent:
            print(f"Shape of Extension             = {d0.shape}")
        
        data['loop']['clocktime'] = d0[:, 2] # shape: (27000,) = (ntimes, clocktime)

        if not silent: 
            for main_key in data:
                for sub_key in data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = '{sub_key}': {data['loop']['clocktime'].shape}")
                    warn_if_zero("loop.clocktime", data['loop']['clocktime'], silent) # checks if array is all zeroes

            print(f"------------------------------------------------\n")

            print(f"Loading Extension 3: Centroids")
            print(f"------------------------------------------------")
        
        d3 = hdul[3].data  # shape: (27000, 1, 288) = (ntimes, nwfs, ncentroids)
        # ncentroids = 144 x-slopes + 144 y-slopes
        # 144 = 12x12 nsub
        
        if not silent: 
            print(f"Shape of Extension             = {d3.shape}")

        data['wfs']['xcentroids'] = d3[:, :, :nsub] # shape: (27000, 1, 144) = (ntimes, nwfs, xcentroids)
        data['wfs']['ycentroids'] = d3[:, :, nsub:] # shape: (27000, 1, 144) = (ntimes, nwfs, ycentroids)
        
        if not silent: 
            for main_key in list(data.keys())[1:]:
                for sub_key in data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = {data[main_key][sub_key].shape}")
                    warn_if_zero("loop.clocktime", data['loop']['clocktime'], silent) # checks if array is all zeroes
            print(f"------------------------------------------------\n")
            
            print(f"Loading Extension 4: DM Voltages")
            print(f"------------------------------------------------")
        
        d4 = hdul[4].data # shape: (27000, 2, 64) = (ntimes, voltage_mes, nact_driver)
        # voltage_mes = deltav, voltage
        # nact = 36 actuators with 36 padded zeroes 
        
        if not silent: 
            print(f"Shape of Extension             = {d4.shape}")
        
        data['dm']['deltav'] = d4[:, 0, :nact] # shape: (27000, 36) = (ntimes, deltav, nact)
        data['dm']['voltage'] = d4[:, 1, :nact] # shape: (27000, 36) = (ntimes, voltage, nact)
        
        if not silent: 
            for main_key in list(data.keys())[2:]:
                for sub_key in data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = {data[main_key][sub_key].shape}")
                    warn_if_zero("loop.clocktime", data['loop']['clocktime'], silent) # checks if array is all zeroes
            print(f"------------------------------------------------\n")
        
    return data

def warn_if_zero(name, 
                 arr, 
                 silent):
    
    if not silent: 
        if np.all(arr == 0):
            print(f"WARNING: The array '{name}' contains all zeroes.")
            
    return 

