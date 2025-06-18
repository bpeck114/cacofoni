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
    print(f"Telemetry file         = {ptele}")
    print(f"Parameter file         = {pparam}")
    
    #print("\nAssumptions:") # From configuration file 
    #print(f"Number of actuators    = {config.num_actuators}")
    #print(f"Maximum number of wfs  = {config.nwfs_max} ")
    #print(f"Minimum frequency      = {config.minimum_frequency}")
    #print(f"Maximum frequency      = {config.maximum_frequency}")
    #print(f"Sampling Frequency.    = {config.sampling_frequency}")
        
    # 1) Load the FITS telemetry and parameter structure with irdfits
    # Makes an empty structure and fills it with telemetry data
    exten = config.extension
    cb = irdfits(ptele, pparam, exten)
    print("\nFinished loading telemetry data...")
    
    #cb_dm = cb[0]['dm']
    
    # Checks if closed is True
    # If True, uses deltav 
    # If False, uses voltages 
    #if closed:
        #print("\nMode (C/O): Closed Loop (using deltav)")
        #com = cb_dm['deltav'][0:config.num_actuators, :]
        
    #else:
        #print("\nMode (C/O): Open Loop (using voltages)")
        #com = cb_dm['voltages'][0:config.num_actuators, :]  
    
    mes1 = np.array([cb[t]['wfs'][0]['centroids'][:, 0] for t in range(len(cb))]).T
    nsub, nsamp = mes1.shape
    mes1_centered = mes1 - mes1.mean(axis=1, keepdims=True) 
    # shape: (nsub, nsamp) = (288, 27000)
    
    freq = (np.arange(nsamp // 2) + 1) / (nsamp / 2) * (config.sampling_frequency / 2)
    filter_mask = (freq >= config.minimum_frequency) & (freq <= config.maximum_frequency) # bandpass filter
    filter1 = np.tile(filter_mask, (config.num_actuators, 1))
    window = np.hanning(nsamp)
    specmes = np.empty((288, nsamp), dtype=np.complex64)
    for i in range(288):
        specmes[i, :] = np.fft.fft(window * mes1[i, :])
        
    psdmes = specmes[:, 0:(nsamp // 2)].astype(np.float32)
    
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
