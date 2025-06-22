# FILE: make_cacofoni_dev.py

# Import packages
import numpy as np

from cacofoni.config import CacofoniConfig
from cacofoni.utils.file_utils import get_valid_path
from cacofoni.imaka_io.load_telemetry import load_telemetry
from cacofoni.utils.calc_samp_freq import calc_samp_freq

def make_cacofoni_dev(ftele=None,
                      fparam=None,
                      fmirror=None,
                      silent=False):
    
    if not silent:
        print("Setting up make_cacofoni...\n")
    
    config = CacofoniConfig()
    nact = config.num_actuators
    nsub = config.num_subapertures
    fmin = config.minimum_frequency
    fmax = config.maximum_frequency
    fstart = config.start_frequency
        
    if not silent:
        print(f"Assumptions from configuration file:")
        print(f"------------------------------------------------")
        print(f"Number of actuators            = {nact}")
        print(f"Number of subapertures         = {nsub}")
        print(f"Minimum frequency (Hz)         = {fmin}")
        print(f"Maximum frequency (Hz)         = {fmax}")
        print(f"Starting frequency (Hz)        = {fstart}")
        print(f"Closed/Open Loop               = {'Closed' if config.closed else 'Open'}")
        print(f"Modal/Zonal                    = {'Modal' if config.modal else 'Zonal'}")
        print(f"Laplacian?                     = {'Yes' if config.laplacian else 'No'}")
        print(f"------------------------------------------------\n")
        
        print("Loading imaka telemetry file...\n")
    
    ftele = get_valid_path(ftele, config.telemetry_filename)
    data = load_telemetry(ftele=ftele, silent=True)
    
    # Data from telemetry file for interaction matrix
    clocktime = data['loop']['clocktime']    # shape: (27000,) = (ntimes, clocktime)
    xcentroids = data['wfs']['xcentroids']   # shape: (27000, 1, 144) = (ntimes, nwfs, xcentroids)
    ycentroids = data['wfs']['ycentroids']   # shape: (27000, 1, 144) = (ntimes, nwfs, ycentroids)
    deltav = data['dm']['deltav']            # shape: (27000, 36) = (ntimes, deltav, nact)
    voltage = data['dm']['voltage']          # shape: (27000, 36) = (ntimes, voltage, nact)
    
    # Derived quantities from telemetry files
    nwfs = data['wfs']['xcentroids'].shape[1]
    ntimes = data['loop']['clocktime'].shape[0]
    ncentroids = xcentroids.shape[2] + ycentroids.shape[2]
    fsamp, dt_1_sigma = calc_samp_freq(clocktime) # frequency of system
    total_runtime = int(ntimes * dt_1_sigma) # seconds
    
    if not silent:
        print(f"Assumptions from {ftele} telemetry file:")
        print(f"------------------------------------------------")
        print(f"Number of wavefront sensors    = {nwfs}")
        print(f"Number of steps                = {ntimes}")
        print(f"Number of centroids (x + y)    = {ncentroids}")
        print(f"Sampling frequency (Hz)        = {fsamp}")
        print(f"Approximate total runtime (s)  = {total_runtime}")
        print(f"------------------------------------------------\n")
        
        print("Making centroid matrix...\n")
        
    centroids = np.concatenate((xcentroids, ycentroids), axis=2)  # shape: (ntimes, nwfs, ncentroids)
    centroids_flat = centroids.reshape(ntimes, -1) # shape: (ntimes, nwfs * ncentroids) = (27000, 1 * 288)
    centered_matrix = centroids_flat - np.mean(centroids_flat, axis=1, keepdims=True) # shape: (ntimes, ncentroids) = (27000, 288)

    if not silent: 
        print("Calculating non-zero, positive FFT frequencies...\n")
    
    n_pos_freqs = ntimes // 2 # 27000 // 2 = 13500
    
    # to not have 0 Hz (DC): 0 to 13500
    fft_bin_indices = np.arange(n_pos_freqs) + 1 # shape: (nfreqs,) = (13500,)
    normalized_freqs = fft_bin_indices / n_pos_freqs # shape: (nfreqs,) = (13500,)
    
    nyquist_freq = fsamp / 2 # Nyquist frequency = 499 Hz
    positive_freqs = normalized_freqs * nyquist_freq # shape: (nfreqs,) = (13500,)
    
    mask = (positive_freqs >= fmin) & (positive_freqs <= fmax) # shape: (nfreqs,) = (13500,)
    matrix_mask = np.tile(mask[:, np.newaxis], (1, nact)) # shape: (nfreqs, nact) = (13500, 36)
    
    fft_centroids = np.zeros((ncentroids, ntimes)) # shape: (ncentroids, ntimes) = (288, 27000)
    
    for i in range(ncentroids):
        time_series = centered_matrix[:, i] # shape: (ntimes, ncentroids) = (27000, 288)
        hann_window = np.hanning(ntimes) # shape: (ntimes,) = (27000,)
        windowed_time_series = hann_window * time_series # shape: (ntimes,) = (27000,)
        fft_result = np.fft.fft(windowed_time_series) # shape: (ntimes,) = (27000,), complex
        
        # Normalizing to match how idl scales
        normalized_fft = fft_result / ntimes # shape: (ntimes,) = (27000,)
        fft_centroids[i, :] = normalized_fft 
    
    
    return