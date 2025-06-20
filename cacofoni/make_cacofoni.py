# FILE: make_cacofoni.py

# Import packages
import numpy as np
import matplotlib.pyplot as plt
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
        
    fsamp = config.sampling_frequency
    minfreq = config.minimum_frequency
    maxfreq = config.maximum_frequency
    num_actuators = config.num_actuators
    
    positive_freqs = (np.arange(n_timesteps // 2) + 1) / (n_timesteps / 2) * (fsamp / 2) # shape: (13500,)
    
    minfreq = config.minimum_frequency
    maxfreq = config.maximum_frequency
    
    freq_band_mask = (positive_freqs >= minfreq) & (positive_freqs <= maxfreq)  # Boolean mask, shape: (13500,)
    
    freq_band_mask_2d = np.tile(freq_band_mask[:, np.newaxis], (1, config.num_actuators))  # shape: (13500, 36)

    spec_centroids = np.zeros((n_centroids, n_timesteps), dtype=np.complex64)
    
    for i in range(n_centroids):
        spec_centroids[i, :] = np.fft.fft(np.hanning(n_timesteps) * centered_centroids[i, :]) / 27000 # To match the idl normalization
        
    psd_centroids = np.real(spec_centroids[:, :n_timesteps // 2]) # shape: (288, 13500)

    speccom = np.zeros((config.num_actuators, n_timesteps), dtype=np.complex64) # shape: (36, 27000)
    
    specmodmod = np.zeros((n_centroids, config.num_actuators, n_timesteps), dtype=np.complex64) # shape: (288, 36, 27000)
    
    cb_dm = [frame['dm'] for frame in telemetry_frames]
    
    if config.closed:
        com = np.array([dm['deltav'][:config.num_actuators] for dm in cb_dm]).T # shape: (36, 27000)
    else:
        com = np.array([dm['voltages'][:config.num_actuators] for dm in cb_dm]).T  # shape: (36, 27000)
        print("Not closed.")
      
    for i in range(config.num_actuators):
        speccom[i, :] = np.fft.fft(np.hanning(n_timesteps) * com[i, :]) / 27000
        for j in range(288):
            specmodmod[j, i, :] = spec_centroids[j, :] / speccom[i, :]
   
    # print("Shape:", specmodmod.shape)
    # print("Min / Max:", np.min(specmodmod), np.max(specmodmod))
    # print("First 5 values:\n", specmodmod[:5, :5])  

    psdmod = np.abs(speccom[:, :n_timesteps // 2])
    
    if config.thresh is None:
        thresh = np.max(psdmod[:, 5:] / 20.0)
    
    indexmax = np.zeros(config.num_actuators, dtype=np.float32)
    
    for i in range(config.num_actuators):
        if np.max(freq_band_mask * psdmod[i, :]) > thresh:
            idx = np.argmax(np.abs(freq_band_mask * psdmod[i, :]))
            print(i, idx, positive_freqs[idx], psdmod[i, idx])
            indexmax[i] = idx
        else:
            indexmax[i] = 0.
    
    imatcacophony = np.zeros((288, 36), dtype=np.float32)
    
    for i in range(config.num_actuators):
        idx = int(indexmax[i])
        imatcacophony[:, i] = -1.0 * (specmodmod[:, i, idx] + specmodmod[:, i, n_timesteps - idx]).astype(np.float32) / 2.0

    if config.closed:
        imatcacophony *= -1.0

    #if config.modal:
        #imatcacophony = np.matmul(imatcacophony, mod2act)  # or imatcacophony @ mod2act
        
    if config.laplacian:
        inffuncdx = np.zeros((12, 12, config.num_actuators), dtype=np.float32)
        inffuncdy = np.zeros((12, 12, config.num_actuators), dtype=np.float32)
        laplacian = np.zeros((12, 12, config.num_actuators), dtype=np.float32)
        for i in range(config.num_actuators):
            inffuncdx[:, :, i] = imatcacophony[0:144, i].reshape(12, 12)
            inffuncdy[:, :, i] = imatcacophony[144:n_centroids, i].reshape(12, 12)
            laplacian[:, :, i] = (deriv2D(inffuncdy[:, :, i], y=True) + deriv2D(inffuncdx[:, :, i], x=True))


    return laplacian, imatcacophony, idx
