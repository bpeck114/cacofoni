# FILE: params.py 
# Still needs to be tested and made more efficient

from dataclasses import dataclass, field
import numpy as np
from typing import List
from math import pow
from pathlib import Path

# Need to go back through and understand what is going on
def initimakaparmstruct(fparm):
    """
    Initializes a nested parameter structure from
    the imaka parameter file. 
    
    Inputs:
    -------
    fparm : str
          Path to parameter file. 
    
    Outputs:
    --------
    structure : dict 
              A dictionary representing the full 
              imaka parameter structure. 
    
    """
    
    # may not need to be a bool anymore
    success, nsub_vals = getiparm(fparm, "sys_parm.nsub")
    NSUB = int(nsub_vals[0]) if success else 0
    NSUBTOTAL = NSUB * NSUB
    
    sucsess, nact_vals = getiparm(fparm, "sys_parm.nact")
    NACT = int(nact_vals[0]) if success else 0
    
    sucsess, nwfs_vals = getiparm(fparm, "sys_parm.nwfs")
    NWFS = int(nwfs_vals[0]) if success else 0
    
    success, npixx_vals = getiparm(fparm, "wfscam_parm.npixx")
    # assumes line format like "0    120" 
    npixx = [int(line.split()[1]) for line in npixx_vals]
    
    if sum(npixx) != NWFS * npixx[0]:
        raise ValueError("Error in NPIXX")
    NPIXX = npixx[0] 
    NPIXPERSUB = [val / NSUB for val in npixx]
    
    NWFSMAX = 5
    NPIXPERFRAME = pow(NSUB * NPIXPERSUB[0], 2)
    
    # Guide star
    gs_parm = {
        "name": "",
        "mag": 0.0,
        "dRA": 0.0,
        "dDec": 0.0
    }
    
    gs_array = [gs_parm.copy() for _ in range(NWFSMAX)]
    

    return


    # Guide star parameters (gs_parm)
    gs_parm = {
        "name": "",
        "mag": 0.0,
        "dRA": 0.0,
        "dDec": 0.0
    }
    gs_array = [gs_parm.copy() for _ in range(NWFSMAX)]

    # Target parameters (target_parm)
    target_parm = {
        "name": "",
        "tRA": 0.0,
        "tDec": 0.0,
        "fRA": 0.0,
        "fDec": 0.0,
        "nGS": 0,
        "gs": gs_array
    }

    # WFS camera parameters (wfscam_parm)
    wfscam_parm = {
        "camera_sn": 0,
        "npixx": 0,
        "npixy": 0,
        "x0": 0,
        "y0": 0,
        "texp": 0.0,
        "temp": 0.0,
        "emgain": 0,
        "skyname": "",
        "flatname": ""
    }
    wfscam_array = [wfscam_parm.copy() for _ in range(NWFSMAX)]

    # WFS parameters (wfs_parm)
    wfs_parm = {
        "pixelweights": [0.0] * int(NPIXPERSUB[0]),
        "pixelthreshold": 0,
        "xsub": [0.0] * NSUBTOTAL,
        "ysub": [0.0] * NSUBTOTAL,
        "centroidoffsetname": ""
    }
    wfs_array = [wfs_parm.copy() for _ in range(NWFSMAX)]

    # DM parameters (dm_parm)
    dm_parm = {
        "xact": [0.0] * NACT,
        "yact": [0.0] * NACT,
        "voltmax": 0.0
    }

    # Loop parameters (loop_parm)
    loop_parm = {
        "niter": 0,
        "nave": 0,
        "gainp": 0.0,
        "gaini": 0.0,
        "gaind": 0.0,
        "imatname": "",
        "cmatname": [""] * 8,
        "cmat2name": ""
    }

    # System parameters (sys_parm)
    sys_parm = {
        "nwfs": NWFS,
        "nact": NACT,
        "nsub": NSUB,
        "ncb": 0,
        "sim_dm": 0,
        "sim_wfscams": 0
    }

    # Final nested parameter dictionary
    imakaparm = {
        "sys": sys_parm,
        "target": target_parm,
        "wfscam": wfscam_array,
        "wfs": wfs_array,
        "dm": dm_parm,
        "loop": loop_parm
    }

    return imakaparm


def getiparm(fparm, 
             keyword):
    """
    Searchs through the imaka parameter file for certain keywords.
    Returns a list of strings for value associated with keyword.
    
    Inputs:
    -------
    fparm : str
          Path to imaka parameter file.
          Should be the path to imakaparm.txt
    
    keyword : str
            Keyword to search for. 
            Examples: 'sys_parm.nwfs', 'target_parm.gs'

    Outputs:
    --------
    success : bool # May not need this bool feature anymore (come back)
            True if the keyword was found and returns a meaningful value.
            False otherwise
    
    values : list of str
           The values found on lines matching the keyword.
           Keyword is stripped. 
    
    Raises:
    -------
    ValueError
        If the keyword is not in the list of allowed keys. 
    
    """
    
    # 0) Setup
    # List of allowed keywords from imakaparm.txt
    allowed_keywords = {
        "sys_parm.nwfs", "sys_parm.nact", "sys_parm.nsub", "sys_parm.ncb",
        "sys_parm.ncbmax", "sys_parm.ncbskip", "sys_parm.cb_autoname",
        "sys_parm.sim_dm", "sys_parm.sim_wfscams",
        "target_parm.name", "target_parm.tRA", "target_parm.tDec",
        "target_parm.fRA", "target_parm.fDec", "target_parm.nGS", "target_parm.gs",
        "wfscam_parm.camera_sn", "wfscam_parm.npixx", "wfscam_parm.npixy",
        "wfscam_parm.x0", "wfscam_parm.y0", "wfscam_parm.texp", "wfscam_parm.temp",
        "wfscam_parm.emgain", "wfscam_parm.skyname", "wfscam_parm.flatname",
        "wfscam_parm.simcamfile",
        "wfs_parm.pixelweights", "wfs_parm.pixelthreshold",
        "wfs_parm.xsub", "wfs_parm.ysub"
    }
    
    # Raises error if given incorrect keyword. 
    if keyword not in allowed_keywords:
        raise ValueError(f"Keyword '{keyword}' is not in the allowed list for getiparm.")
    
    keyword_len = 30 # default (need to come back and make this more efficient)
    
    # 1) Tries to open the file 
    try:
        with open(fparm, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f'Cannot find {fparm}.')
        
    matches = []
    
    for line in lines:
        if keyword in line:
            parts = line.strip()
            if parts.startswith(keyword):
                matches.append(parts)
   
    values = [line[keyword_len:].strip() for line in matches]
    return True, values


def initimakadatastruct(fparam,
                        ntimes):
    
    """
    Inputs:
    -------
    fparam : str
           Path to parameter file.
    
    ntimes : int
           Number of time steps (frames).
    
    
    Optional Inputs:
    ----------------
    
    Outputs:
    --------
    data : list of dict
         List of telemetry data dictionaries, one per timestep. 
    
    """
    
    

    
    return
 