# FILE: make_cacofoni.py

# Import packages
import numpy as np
from astropy.io import fits
from cacofoni.config import CacofoniConfig
from cacofoni.imaka_io.irdfits import irdfits
from cacofoni.utils.file_utils import get_valid_path


def make_cacofoni(ftele=None, 
                  fparam=None,
                  fmirror=None,
                  config=CacofoniConfig(),
                  closed=False,
                  modal=False,
                  silent=False,
                  thresh=None,
                  compute_laplacian=False):
    
    """


    """
    
    # 0) Setup
    # Load in necessary files 
    print("Setting up make_cacofoni...")
    
    # If user gives a path, checks if path exists
    # If user does NOT give a path, checks if default path exists
    ptele   = get_valid_path(ftele, config.telemetry_filename) # (User defined path to file, default path to file)
    pparam  = get_valid_path(fparam, config.param_filename)
    
    print("\nFile Paths:")
    print(f"Telemetry File:      {ptele}")
    print(f"Parameter File:      {pparam}")
    
    print("\nAssumptions:") # From configuration file 
    print(f"Number of Actuators: {config.num_actuators}")
    print(f"Minimum Frequency:   {config.minimum_frequency}")
    print(f"Maximum Frequency:   {config.maximum_frequency}")
    
    # Checks if modal is True
    # If True, checks if mirror modes file exists (like above)
    if modal:
        print("\nMode (M/Z): Modal") # Loads in mirror mode matrix
        mirror_modes_path = get_valid_path(fmirror, config.mirror_modes_filename)
        
        print(f"Mirror Modes File:  {mirror_modes_path}")
        
        # Inverts to convert from mode amplitudes to actuator voltages
        mirmodes = fits.getdata(mirror_modes_path)
        mod2act = np.linalg.pinv(mirmodes)
    else:
        print("\nMode (M/Z): Zonal")
        
    # 1) Load the FITS telemetry and parameter structure with irdfits
    exten = config.extension
    cb = irdfits(ptele, pparam, exten)
    
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
