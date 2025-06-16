# FILE: make_cacophony.py

# Import packages
import numpy as np
from scipy.fft import fft
from scipy.io import wavfile
from scipy.signal.windows import hann
from astropy.io import fits
import importlib_resources as pkg_resources
from cacofoni.config import CacophonyConfig
from cacofoni.imaka_io.irdfits import irdfits

def make_cacofoni(fname, 
                  minfreq, 
                  maxfreq,
                  config=CacophonyConfig(),
                  closed=False,
                  modal=False,
                  silent=False,
                  thresh=None,
                  compute_laplacian=False):
    
    """
    Generates the interaction matrix from telemetry data.
    
    Inputs:
    -------
    fname : str
             Location of telemetry FITS file. This file contains deformable
             mirror commands and wavefront sensor centroid telemetry. 
            
    minfreq : float
            Minimum frequency to include in the analysis (Hz). 
    
    maxfreq : float
            Maximum frequency to include.
            Should be at least 3.6. Hz higher than the minfreq. 
    
    Optional Inputs:
    ----------------
    config :
    
    
    modal : bool
          If True, use modal control (i.e., via mirror modes like Zernikes).
          If False, use zonal control (i.e., actuator-by-actuator)
    
    closed : bool
           If True, use closed-loop delta voltages from DM telemetry.
           If False, use open-loop absolute voltages. 
    
    silent : bool
           If True, suppress verbose output. 
           
    
    thresh : float or None
           Optional threshold parameter for signal filtering.
    
    laplacian : bool
              If True, operate on curvature (Laplacian) of phase. 
    
    
    Outputs:
    --------
    imat:
         Interaction matrix 
    
    cmat:
    """
    
    # 0) Setup
    # Loading in the configuration file 
    # with default values and file paths
    print("Setting up make_cacofoni...\n")
    
    print("Assumptions from arguments:")
    print(f"Path to file:              {fname}")
    print(f"Minimum Frequency:         {minfreq}")
    print(f"Maximum Frequency:         {maxfreq}")
        
    # 0.1) Assigning variables from CacophonyConfig class
    num_actuators = config.num_actuators
    param_path = config.fparam_path
    mirmodes_path = config.mirror_modes_path
    
    print("\nAssumptions from configuration file:")
    print(f"Number of Actuators:       {num_actuators}")
    print(f"Path to Mirror Modes File: {mirmodes_path}")
    print(f"Path to Parameter File:    {param_path}")
   
    '''
    # 0.2) If modal mode is requested, load in the mirror mode matrix 
    if modal:
        mirmodes = fits.getdata(mirmodes_path)
        mod2act = np.linalg.pinv(mirmodes) 
        # Inverting to convert from mode amplitudes to actuator voltages
    '''  
    
    # 3) Read in the telemetry FITS file with the parameter file
    
    
    return 
    
    """
    hdul = fits.open(filename)
    dm_data = hdul[4].data # shape: (64, 2, 27000)
    wfs_x = hdul[3].data   # shape: (288, 1, 27000)
    wfs_y = hdul[8].data   # shape: (288, 1, 27000)
    
    cb = {
        'dm' : {
            'voltages': dm_data[:, 0, :], # shape: (64, 27000)
            'deltav': dm_data[:, 1, :]    # shape: (64, 27000)
        },
        'wfs': {
            # Stacks X and Y centroids
            'centroids': np.vstack([wfs_x[:, 0, :], wfs_y[:, 0, :]]) # shape: (576, 27000)
            
    """
        
        
   



"""   
function make_cacophony, filename,minfreq,maxfreq,closed=closed,modal=modal,fparm=fparm,silent=silent,thresh=thresh,laplacian=laplacian


  nmodes = 36
  
  if not(keyword_set(fparm)) then begin
     fparam='/home/imaka/python/ichigo-imaka/data/imakaparm.txt'
  endif else fparam=fparm
  if keyword_set(modal) then begin
     mirmodes=readfits('/home/imaka/python/ichigo-imaka/data/mirror_modes_20240409b.nonorm.fits')
     mod2act=invert(mirmodes)
  endif 
  cb=irdfits(filename,fparm=fparam,exten=[1,0,0,1,1,1,1,1])
  cb_dm=cb.dm
  if keyword_set(closed) then begin
     com=cb_dm.deltav[0:35,*]
  endif else com=cb_dm.voltages[0:35,*]
  
"""
