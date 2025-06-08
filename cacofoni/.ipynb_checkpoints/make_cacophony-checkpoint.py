# FILE: make_cacophony.py

# Import packages
import numpy as np
from scipy.fft import fft
from scripy.io import wavfile
from scipy.signal.windows import hann
from astropy.io import fits
import importlib_resources as pkg_resources
from cacofoni.config import CacophonyConfig
from imaka_io.irdfits import irdfits

def make_cacofoni(filename, 
                  minfreq, 
                  maxfreq,
                  config=None,
                  fparm=None,
                  closed=False,
                  modal=False,
                  silent=False,
                  thresh=None,
                  compute_laplacian=False):
    
    """
    Generates the interaction matrix from telemetry data.
    
    Inputs:
    -------
    filename : str
             Location of telemetry FITS file. 
            
    minfreq : float
            Minimum frequency. 
    
    maxfreq : float
            Maximum frequency, should be 3.6. Hz higher than minimum frequency. 
    
    Optional Inputs:
    ----------------
    config :
    
    
    modal : bool
             Modal = modulates by Zernike modes (tiptilt, focus, etc).
             Zonal = modulates actuator by actuator.
    
    closed : bool
            closed = closed loop. 
            open = open loop. 
    
    
    
    fparm : str
    
    silent : bool
    
    thresh : float 
    
    laplacian : bool
              Takes the curvature of the phase. 
    
    
    Outputs:
    --------
    imat:
    
    cmat:
    """
    
    # Loading in the configuration file for default values and file paths
    if config is None:
        config = CacophonyConfig()
    
    # Assigning variables from CacophonyConfig class
    num_actuators = config.num_actuators
    fparam = config.fparam_path
    mirmodes_path = config.mirror_modes_path
    
    # Allowing for overrides
    if fparm is not None:
        fparam = fparm
        
    if modal:
        mirmodes = fits.getdata(mirmodes_path)
        mod2act = np.linalg.pinv(mirmodes) 
        
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
            

        
        return 
   



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
