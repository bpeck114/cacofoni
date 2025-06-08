# FILE: irdfits.py
# Still needs to be tested 

# Import packages
from astropy.io import fits
import numpy as np
from cacofoni.config import CacophonyConfig

def irdfits(filename, fparm=None, silent=False, exten=None):
    """
    Reads in a FITS file.
    
    Inputs:
    -------
    filename : str
             Path to FITS file.
    
    Optional Inputs:
    ----------------
    fparm : str
          Path to parameter file.
          Default is located in: "cacofoni/data/imakaparm/txt"
          
    silent : bool
           Whether to suppress warning messages.
    
    exten : 
          Which FITS extensions to read 
          (list or a single integer)
          Default: loads all 8 extensions.
    
    Outputs:
    ---------
    
    """
    
    # If not parameter file is chosen, go with the default
    if fparm is None:
        fparm = CacophonyConfig.fparam_path
        print(f"Parameter File: {CacophonyConfig.fparam_path}")
    
    # If "exten" was not specified, load all 8 extensions by default
    if exten is None:
        exten = [1, 1, 1, 1, 1, 1, 1] # Default
    elift isinstance(exten, int): # checking if exten is an integer
        # E.g., exten=4 then builds a list like [0, 0, 0, 0, 1, 0, 0, 0]
        tmp = [0] * 8
        tmp[0] = 1
        tmp[exten] = 1
        exten = tmp 
        
    # Open the FITS file
    with fits.open(fname) as hdul:
        h0 = hdul[0].header # Read the primary header
        ntimes = h0['NAXIS2'] # NAXIS2 = number of time samples in telemetry
        
        # Making a placeholder list of "ntimes" in empty dictionaries
        data = [{} for _ in range(ntimes)]
        dm_data = {} # deformable mirror data
        
        # If exten=4 is specified then:
        if exten[4]:
            d = hdul[4].data # Shape assumed: [n_actuators, 2, ntimes]
            dm_data['deltav'] = d[:, 0, :]
            dm_data['voltages'] = d[:, 1, :]
            
        # If exten=7 is specified then read average DM voltages
        if exten[7]:
            try:
                d = hdul[7].data
                dm_data['avevoltages'] = d
            except Exception as e:
                if not silent:
                    print("No average DM data save:", e)
                    
        return {'dm': dm_data}
                