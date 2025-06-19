# FILE: make_cacofoni.py

# Import packages
import numpy as np
from astropy.io import fits
from cacofoni.config import CacofoniConfig
from cacofoni.imaka_io.irdfits import irdfits
from cacofoni.utils.file_utils import get_valid_path
from cacofoni.utils.deriv2D import deriv2D 


def make_cacofoni(ftele=None, 
                  fparam=None,
                  fmirror=None,
                  silent=False):
    
    """
    """
    
    # 0) Setup
    if not silent:
        print("Setting up make_cacofoni...")
    
    config = CacofoniConfig() # Configuration file that holds assumptions
    ptele   = get_valid_path(ftele, config.telemetry_filename)
    pparam  = get_valid_path(fparam, config.param_filename)
    
    
    if config.modal:
        pmirror  = get_valid_path(fparam, config.mirror_modes_filename)
    
    if not silent:
        print("Assumptions from configuration file:") # From configuration file
        print("Override these values in the configuration file.")
        print("------------------------------------------------")
        print(f"Minimum frequency (Hz)      = {config.minimum_frequency}")
        print(f"Maximum frequency (Hz)      = {config.maximum_frequency}")
        print(f"Sampling Frequency (Hz)     = {config.sampling_frequency}")
        print(f"Maximum number of WFS       = {config.nwfs_max}")
        print(f"Number of actuators         = {config.num_actuators}")
        print(f"Closed/Open Loop            = {'Closed' if config.closed else 'Open'}")
        print(f"Modal/Zonal                 = {'Modal' if config.modal else 'Zonal'}")
        print(f"Laplacian?                  = {'Yes' if config.laplacian else 'No'}")
        print("------------------------------------------------\n")
        
        print("File Paths:")
        print("Override these paths in make_cacofoni arguments.")
        print("------------------------------------------------")
        
        print(f"Telemetry file             = {ptele}")
        print(f"Parameter file             = {pparam}")
        
        if config.modal:
            print(f"Mirror Modes file:         = {pmirror}")
        else: 
            print(f"Mirror Modes file:         = N/A")
        print("------------------------------------------------\n")
        
        
    # 1) Load the FITS telemetry and parameter structure with irdfits
    # Makes an empty structure and fills it with telemetry data
    
    if not silent:
        print("Loading telemetry data with irdfits...\n")
    
    exten = config.extension
    telemetry_frames = irdfits(ptele, pparam, exten=exten, silent=silent)
    
    if not silent:
        print("Finished loading telemetry data...")
        print("Loading WFS centroid measurements...\n")
        
    wfs_data_per_frame = [frame['wfs'] for frame in telemetry_frames]
    centroid_matrix = np.array([wfs_frame[0]['centroids'] for wfs_frame in wfs_data_per_frame])
    centroid_matrix = centroid_matrix.T # transposing to match idl code
    centered_centroids = centroid_matrix - np.mean(centroid_matrix, axis=1, keepdims=True) 
    n_centroids, n_timesteps = centered_centroids.shape
   
    if not silent:
        print(f"Found {n_centroids} centroid measurements for {n_timesteps}.\n")
        
    print("Shape:", centered_centroids.shape)                        # equivalent to IDL's SIZE
    print("Min / Max:", np.min(centered_centroids), np.max(centered_centroids))    # MIN and MAX
    print("Top-left 5x5 block:\n", centered_centroids[:5, :5])       # preview 5x5

        
    fsamp = config.sampling_frequency
    minfreq = config.minimum_frequency
    maxfreq = config.maximum_frequency
    num_actuators = config.num_actuators
    
    # Step 3: Create frequency array (match IDL’s +1 behavior and 996 Hz assumption)
    positive_freqs = (np.arange(n_timesteps // 2) + 1) / (n_timesteps / 2) * (fsamp / 2)

    # Step 4: Frequency band mask
    minfreq = config.minimum_frequency
    maxfreq = config.maximum_frequency
    freq_band_mask = (positive_freqs >= minfreq) & (positive_freqs <= maxfreq)
    freq_band_mask_2d = np.tile(freq_band_mask[:, np.newaxis], (1, config.num_actuators))  # shape: (13500, 36)

    # Step 5: Hann window and FFT
    hann_window = np.hanning(n_timesteps)
    spec_centroids = np.array([
        np.fft.fft(hann_window * centered_centroids[i, :])
        for i in range(n_centroids)
    ])

    # Step 6: PSD (abs value), include upper Nyquist bin
    psd_centroids = np.abs(spec_centroids[:, :n_timesteps // 2 + 1])  # shape (288, 13501)

    # Print debug info (match IDL's print behavior)
    print("Shape:", psd_centroids.shape)
    print("Min / Max:", np.min(psd_centroids), np.max(psd_centroids))
    print("Top-left 5x5 block:\n", psd_centroids[:5, :5])
    
    # Checks if closed is True
    # If True, uses deltav 
    # If False, uses voltages 
    # Makes 64 values because of the drivers but
    # there is actually only 36 actuators.
    
    cb_dm = [frame['dm'] for frame in cb]
    if config.closed:
        # shape will be (36, ntimes)
        com = np.array([dm['deltav'][:config.num_actuators] for dm in cb_dm]).T
        print(f"len com: {len(com)}")
        
    else:
        com = np.array([dm['voltages'][:config.num_actuators] for dm in cb_dm]).T 
        print(f"len com: {len(com)}")


    
    # Checks if modal is True
    # If True, checks if mirror modes file exists (like above)
    
    '''
    if modal:
         # Loads in mirror mode matrix
        mirror_modes_path = get_valid_path(fmirror, config.mirror_modes_filename)
        
        print(f"Mirror Modes File:  {mirror_modes_path}")
        
        # Inverts to convert from mode amplitudes to actuator voltages
        mirmodes = fits.getdata(mirror_modes_path)
        mod2act = np.linalg.pinv(mirmodes)
    else:
        print("\nMode (M/Z): Zonal")
    '''
    
    
    
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
