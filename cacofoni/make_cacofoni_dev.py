# FILE: make_cacofoni_dev.py

# Import packages
import numpy as np

from cacofoni.config import CacofoniConfig
from cacofoni.utils.file_utils import get_valid_path
from cacofoni.imaka_io.load_telemetry import load_telemetry
from cacofoni.utils.calc_samp_freq import calc_samp_freq
from cacofoni.utils.deriv2D import deriv2D

def make_cacofoni_dev(ftele=None,
                      fparam=None,
                      fmirror=None,
                      hanning=False,
                      silent=False,
                      debug=True):
    """
    Brianna Notes:
        Centroid = position of a spot on the wavefront sensor 
            - We take the data over time (27000 steps)
            - We have the x- and y-coordinate on the wavefront sensor
            
        Power spectrum = how much power at that frequency
        Hann window = reduce edge effects 
        
        To make python match idl better:
            fsamp = 996 
            hann_window = hann_window * (ntimes / np.sum(hann_window)) # matches idl better
            
        Figure out the DC stuff?
        Ask about thresh, why 5? Why 20?
        
        Why taking the float in one place (psdmes) but abs in the other (psdmod)?
        
        Translations for idl variable names:
            freq = freq_pos
            filter = freq_mask
            filter1 = freq_mask_2d
            
            specmes = fft_centroids_all
            psdmes = fft_centroids_real
            
            specmodmod = fft_actuators_all
            
            
    """
    
    if not silent:
        print("Setting up make_cacofoni...\n")
    
    # Assumed configuration from configuration file
    config = CacofoniConfig()
    nact = config.num_actuators # 36
    num_dimx_subap = config.num_dimx_subapertures # 12
    num_dimy_subap = config.num_dimx_subapertures # 12
    fmin = config.minimum_frequency # 4
    fmax = config.maximum_frequency # 10
    fstart = config.start_frequency # 5
    nsub = num_dimx_subap *  num_dimy_subap
    
    # Getting file paths
    ftele = get_valid_path(ftele, config.telemetry_filename)
    
    if not silent:
        print(f"Assumptions from configuration file:")
        print(f"------------------------------------------------")
        print(f"Number of actuators            = {nact}")
        print(f"Number of subapertures in x    = {num_dimx_subap}")
        print(f"Number of subapertures in y    = {num_dimy_subap}")
        print(f"Number of total subapertures   = {nsub}")
        print(f"Minimum frequency (Hz)         = {fmin}")
        print(f"Maximum frequency (Hz)         = {fmax}")
        print(f"Starting frequency (Hz)        = {fstart}")
        print(f"Closed/Open Loop               = {'Closed' if config.closed else 'Open'}")
        print(f"Modal/Zonal                    = {'Modal' if config.modal else 'Zonal'}")
        print(f"Laplacian?                     = {'Yes' if config.laplacian else 'No'}")
        print(f"Telemetry File Path            = {ftele}")
        print(f"------------------------------------------------\n")
        
        print("Loading imaka telemetry file...\n")
    
    # Getting data from telemetry file for interaction matrix
    data = load_telemetry(ftele=ftele, silent=True)
    
    clocktime = data['loop']['clocktime']    # shape: (27000,) = (ntimes, clocktime)
    xcentroids = data['wfs']['xcentroids']   # shape: (27000, 1, 144) = (ntimes, nwfs, xcentroids)
    ycentroids = data['wfs']['ycentroids']   # shape: (27000, 1, 144) = (ntimes, nwfs, ycentroids)
    deltav = data['dm']['deltav']            # shape: (27000, 36) = (ntimes, deltav, nact) # Note: only if closed
    voltage = data['dm']['voltage']          # shape: (27000, 36) = (ntimes, voltage, nact) # Note: only if open
    
    # Derived quantities from telemetry files
    nwfs = data['wfs']['xcentroids'].shape[1]
    ntimes = data['loop']['clocktime'].shape[0]
    ncentroids = xcentroids.shape[2] + ycentroids.shape[2]
    
    fsamp, dt_1_sigma = calc_samp_freq(clocktime) # frequency of the AO system
    hann_window = np.hanning(ntimes) # calcuating hann window
    total_runtime = int(ntimes * dt_1_sigma) # seconds
    
    if debug is True:
        # For debugging, matches idl better 
        fsamp = 996  
        tmp = np.arange(ntimes, dtype=np.float32)
        hann_idl = 0.5 * (1.0 - np.cos(2.0 * np.pi * tmp / (ntimes - 1)))
        hann_window = hann_idl.astype(np.float32)
    
    nyquist_freq = fsamp / 2 # nyquist frequency
    
    if not silent:
        print(f"Assumptions from {ftele} telemetry file:")
        print(f"------------------------------------------------")
        print(f"Number of wavefront sensors    = {nwfs}")
        print(f"Number of steps                = {ntimes}")
        print(f"Number of centroids (x + y)    = {ncentroids}")
        print(f"Sampling frequency (Hz)        = {fsamp}")
        print(f"Nyquist frequncy (Hz)          = {nyquist_freq}")
        print(f"Approximate total runtime (s)  = {total_runtime}")
        print(f"------------------------------------------------\n")
        
        print(f"Combining x/y centroid arrays...") 
        
    centroids = np.concatenate((xcentroids, ycentroids), axis=2) 
    # shape: (ntimes, nwfs, ncentroids = (27000, 1, 288)
    
    if debug is True:
        centroids = centroids.astype(np.float32) # to match idl better
    
    if not silent:
        print(f"    x   : {xcentroids.shape}") # shape: (27000, 1, 144)
        print(f"    y   : {ycentroids.shape}") # shape:  (27000, 1, 144)
        print(f"    x/y : {centroids.shape}\n")  # shape: (27000, 1, 288)
    
    centroids_flat = centroids.reshape(ntimes, -1) 
    # shape: (ntimes, nwfs * ncentroids) = (27000, 1 * 288) = (27000, 288)
        
    if not silent:
        print(f"Flattening centroid array from {centroids.shape} to {centroids_flat.shape}...\n")
        print("Centering centroid array (removing DC)....\n")
    
    centroid_means = np.mean(centroids_flat, axis=0, keepdims=True)
    centroids_centered = centroids_flat - centroid_means
        
    if not silent:
        print(f"Computing FFT for {ncentroids} centroids...")
     
    # Empty array to hold FFT centroid data
    fft_centroids_all = np.zeros((ntimes, ncentroids), dtype=np.complex64) # shape: (ntimes, ncentroids) = (27000, 288) 
    
    for i in range(ncentroids):
        # Take the position of the centroid at all timesteps
        centroid_signal = centroids_centered[:, i] # shape at i=0: (ntimes,) = (27000,)
        
        if hanning:
            if i == 0:
                print(f"Applying a Hann window filter to centroids to reduce edge effects...\n")
            centroid_signal = centroid_signal * hann_window # shape at i=0: (ntimes,) = (27000,)
        
        # Compute FFT
        fft_centroid_unnorm = np.fft.fft(centroid_signal) # shape at i=0: (ntimes,) = (27000,)
        
        # Normalize FFT like idl
        fft_centroids_all[:, i] = fft_centroid_unnorm / ntimes # shape at i=0: (ntimes,) = (27000,)
    
    if not silent:
        print(f"    FFT centroids shape   : {fft_centroids_all.shape}\n") # shape: (ntimes, ncentroids) = (27000, 288) 
        print(f"Computing FFT for {nact} actuators...")
        
    fft_actuators_all = np.zeros((ntimes, nact), dtype=np.complex64) # shape: (ntimes, nact) = (27000, 36)
    
    for i in range(nact):
        actuator_signal = voltage[:, i] # shape at i=0: (ntimes,) = (27000,)
        
        if hanning:
            if i == 0:
                print(f"Applying a Hann window filter to actuators reduce edge effects...\n")
            actuator_signal = actuator_signal * hann_window # shape at i=0: (ntimes,) = (27000,)
        
        fft_actuator_unnorm = np.fft.fft(actuator_signal) # shape at i=0: (ntimes,) = (27000,)
        fft_actuators_all[:, i] = fft_actuator_unnorm / ntimes # shape at i=0: (ntimes,) = (27000,)
    
    if not silent:
        print(f"    FFT actuators shape   : {fft_actuators_all.shape}\n") 
        # shape: (ntimes, nact) = (27000, 36)
        
        print(f"Computing response matrix...")
        
    fft_response_all = np.zeros((ntimes, ncentroids, nact), dtype=np.complex64) 
    # shape: (ntimes, ncentroids, nact) = (27000, 288, 36)
        
    for i in range(ncentroids):
        for j in range(nact):
            fft_response_all[:, i, j] = fft_centroids_all[:, i] / fft_actuators_all[:, j] 
    
    if not silent: 
        print(f"    FFT response shape   : {fft_response_all.shape}\n") 
        # shape: (ntimes, nact) = (27000, 288, 36) 
        
        print(f"Extracting positive or real frequencies...\n")
        
    n_pos_freq_bins = ntimes // 2 # 27000 / 2 = 13500 positive frequencies
    freq_pos = (np.arange(n_pos_freq_bins) + 1) / n_pos_freq_bins * nyquist_freq 
    # shape: (n_pos_freq_bins) = (13500,)
    
    if debug is True:
        freq_pos = freq_pos.astype(np.float32)    
        
    freq_mask = (freq_pos >= fmin) & (freq_pos <= fmax) # shape: (13500,)
    # freq_band_mask_2d = np.tile(freq_mask[:, np.newaxis], (1, nact)) # shape: (13500, 36)
    fft_actuators_pos = np.abs(fft_actuators_all[0:n_pos_freq_bins, :]) # shape: (13499, 36) 
    fft_centroids_real = np.real(fft_centroids_all[0:n_pos_freq_bins, :])
    thresh = np.max(fft_actuators_pos[5:, :]) / 20.0 # arbitrary?
    
    if not silent: 
        print(f"Threshold                      = {thresh}\n")
        print("Preparing for interaction matrix...\n")
    
    peak_freq_indices = np.zeros(nact, dtype=np.int32) # shape: (nact,) = (36,)
    for i in range(nact):
        # Apply frequency mask to PSD of the actuator
        psd_filtered = freq_mask * fft_actuators_pos[:, i] # shape: (ntimes,) = (13500,)
        
        if np.max(psd_filtered) > thresh:
            idx = np.argmax(psd_filtered) 
            peak_freq_indices[i] = idx
            
            if not silent:
                print(f"Mode {i}: index {idx}, freq {freq_pos[idx]}, psd {fft_actuators_pos[idx, i]}")
            
            peak_freq_indices[i] = np.argmax(psd_filtered) # shape: (nact,) = (36,)
        else:
            peak_freq_indices[i] = 0 # shape: (nact,) = (36,)
    
    if not silent:
        print("Computing interaction matrix...\n")
        
    interaction_matrix = np.zeros((ncentroids, nact), dtype=np.float32) # shape: (ncentroids, nact) = (288, 36)
    for i in range(nact):
        positive_idx = peak_freq_indices[i]
        negative_idx = ntimes - positive_idx
        
        if negative_idx >= ntimes:
            negative_idx = ntimes -1
        
        interaction_matrix[:, i] = -1.0 * (
            fft_response_all[positive_idx, :, i].real +
            fft_response_all[negative_idx, :, i].real
        ) / 2.0
            
    if not silent:
        print("Computing laplcian...")
            
    if config.laplacian:
        inffuncdx = np.zeros((num_dimx_subap, num_dimy_subap, nact), dtype=float) # shape: (12, 12, 36)
        inffuncdy = np.zeros((num_dimx_subap, num_dimy_subap, nact), dtype=float) # shape: (12, 12, 36)
        laplacian = np.zeros((num_dimx_subap, num_dimy_subap, nact), dtype=float) # shape: (12, 12, 36)


        for i in range(nact):
            inffuncdx[:, :, i] = interaction_matrix[0:ncentroids//2, i].reshape(12, 12, order='F')
            inffuncdy[:, :, i] = interaction_matrix[ncentroids//2:ncentroids, i].reshape(12, 12, order='F')
   
            laplacian[:, :, i] = (
                deriv2D(inffuncdy[:, :, i], y=True) +
                deriv2D(inffuncdx[:, :, i], x=True)
            )
            
    return freq_pos, fft_centroids_real, fft_actuators_pos, freq_mask, interaction_matrix, laplacian