# FILE: make_cacofoni.py

# Import packages
import numpy as np
from astropy.io import fits

from cacofoni.config import CacofoniConfig
from cacofoni.calc_samp_freq import calc_samp_freq
from cacofoni.deriv2D import deriv2D
from cacofoni._load_and_log_telemetry import _load_and_log_telemetry
from cacofoni.file_utils import check_default_path, check_user_path, resolve_filepath
from cacofoni._compute_fft import _compute_fft
from cacofoni.find_peak_freq import find_peak_actuator_frequencies

def make_cacofoni(telemetry_filepath=None,
                  fparam=None,
                  modal_filepath=None,
                  modal=False,
                  apply_hanning=None,
                  closed=None,
                  silent=False,
                  debug=None,
                  config=None):
    """   
    """
    
    if not silent:
        print("Setting up make_cacofoni...\n")
        
    config = config or CacofoniConfig()
    if not silent:
        print(f"[Config] Assuming {config.n_actuators} actuators from config for loading telemetry data.")
        print(f"[Config] Assuming {config.n_xsubapertures} 'x' subapertures from config for loading telemetry data.")
        print(f"[Config] Assuming {config.n_ysubapertures} 'y' subapertures from config for loading telemetry data.\n")
        
        print(f"[Config] Assuming {config.minimum_freq_hz} Hz for minimum frequency.")
        print(f"[Config] Assuming {config.maximum_freq_hz} Hz for maximum frequency.\n")
    
    apply_hanning = config.apply_hanning if apply_hanning is None else apply_hanning
    silent = config.silent if silent is None else silent
    modal = config.modal if modal is None else modal
    
    telemetry_filepath = resolve_filepath(telemetry_filepath, config.telemetry_filename, silent)   
    telemetry_data = _load_and_log_telemetry(telemetry_filepath, silent)
    xcentroids, ycentroids, deltav, voltage, n_steps, n_act, n_centroids, sampling_freq_hz, nyquist_freq_hz = telemetry_data
    commands = deltav if config.closed else voltage
    
    if modal:
        modal_filepath = resolve_filepath(modal_filepath, config.modal_filename, silent) 
    
    if debug:
        # For debugging, matches idl better 
        sampling_freq_hz = 996  
        tmp = np.arange(n_steps, dtype=np.float32)
        hann_idl = 0.5 * (1.0 - np.cos(2.0 * np.pi * tmp / (n_steps - 1)))
        hann_window = hann_idl.astype(np.float32)
        nyquist_freq = sampling_freq_hz / 2
        
    centroids = np.concatenate((xcentroids, ycentroids), axis=2).astype(np.float32)
    centroids_flat = centroids.reshape(n_steps, -1) 
    centroid_means = np.mean(centroids_flat, axis=0, keepdims=True)
    centroids_centered = centroids_flat - centroid_means
    
    print("centroids_centered.shape", centroids_centered.shape)
    print("np.min(centroids_centered)", np.min(centroids_centered))
    print("np.max(centroids_centered)", np.max(centroids_centered))
    print("centroids_centered[0:5, 0]", centroids_centered[0:5, 0])
    print("centroids_centered[0, 0:5]", centroids_centered[0, 0:5])
        
    fft_centroids_all = _compute_fft(centroids_centered, apply_hanning, hann_window)
    
    print("fft_centroids_all.shape", fft_centroids_all.shape)
    print("np.min(fft_centroids_all)", np.min(fft_centroids_all))
    print("np.max(fft_centroids_all)", np.max(fft_centroids_all))
    print("fft_centroids_all[0:5, 0]", fft_centroids_all[0:5, 0])
    print("fft_centroids_all[0, 0:5]", fft_centroids_all[0, 0:5])
    
    if modal:
        if not silent:
            print("Loading modal file...\n")
    
        with fits.open(modal_filepath) as hdul:
            mirmodes = hdul[0].data # shape: (n_modes, n_act)? = (36, 36)
            
        mod2act = np.linalg.inv(mirmodes) # shape (36, 36)
        modcom = np.dot(commands, mirmodes.T)
        
        fft_actuators_all = np.zeros((n_steps, n_act), dtype=np.complex64)
        fft_response_all = np.zeros((n_steps, n_centroids, n_act), dtype=np.complex64)
            
        for i in range(n_act):
            fft_actuators_all = np.fft.fft(modcom, axis=0).astype(np.complex64)
            
        print("fft_actuators_all.shape", fft_actuators_all.shape)
        print("np.min(fft_actuators_all)", np.min(fft_actuators_all))
        print("np.max(fft_actuators_all)", np.max(fft_actuators_all))
        print("fft_actuators_all[0:5, 0]", fft_actuators_all[0:5, 0])
        print("fft_actuators_all[0, 0:5]", fft_actuators_all[0, 0:5])
        
        fft_response_all = fft_centroids_all[:, :, np.newaxis] / fft_actuators_all[:, np.newaxis, :]
        
    else:
        fft_actuators_all = _compute_fft(commands, apply_hanning, hann_window)
        fft_response_all = fft_centroids_all[:, :, np.newaxis] / fft_actuators_all[:, np.newaxis, :]
        
    print("fft_response_all.shape", fft_response_all.shape)
    print("np.min(fft_response_all)", np.min(fft_response_all))
    print("np.max(fft_response_all)", np.max(fft_response_all))
    print("fft_response_all[0:5, 4, 4]", fft_response_all[0:5, 4, 4])
    print("fft_response_all[4, 0:5, 4]", fft_response_all[4, 0:5, 4])
    print("fft_response_all[4, 4, 0:5]", fft_response_all[4, 4, 0:5])
        
    n_pos_freq_bins = n_steps // 2
    freq_pos = ((np.arange(n_pos_freq_bins) + 1) / n_pos_freq_bins * nyquist_freq).astype(np.float32) 
        
    freq_mask = (freq_pos >= config.minimum_freq_hz) & (freq_pos <= config.maximum_freq_hz)
    actuator_fft_magnitude = np.abs(fft_actuators_all[0:n_pos_freq_bins, :]) 
    fft_centroids_real = np.real(fft_centroids_all[0:n_pos_freq_bins, :])
    
    if not silent: 
        print("Preparing for interaction matrix...\n")
        
    peak_freq_indices = find_peak_actuator_frequencies(actuator_fft_magnitude, freq_mask, freq_pos, silent)
    
    print("peak_freq_indices.shape", peak_freq_indices.shape)
    print("np.min(peak_freq_indices)", np.min(peak_freq_indices))
    print("np.max(peak_freq_indices)", np.max(peak_freq_indices))
    print("peak_freq_indices[0:5]", peak_freq_indices[0:5])
    
    if not silent:
        print("Computing interaction matrix...\n")
        
    interaction_matrix = np.zeros((n_centroids, n_act), dtype=np.float32)
    for i in range(n_act):
        positive_idx = peak_freq_indices[i]
        negative_idx = n_steps - positive_idx
        
        if negative_idx >= n_steps:
            negative_idx = n_steps -1
        
        interaction_matrix[:, i] = -1.0 * (
            fft_response_all[positive_idx, :, i].real +
            fft_response_all[negative_idx, :, i].real
        ) / 2.0
        
    print("interaction_matrix.shape", interaction_matrix.shape)
    print("np.min(interaction_matrix)", np.min(interaction_matrix))
    print("np.max(interaction_matrix)", np.max(interaction_matrix))
    print("interaction_matrix[0:5]", interaction_matrix[0, 0:5])
        
    if modal:
        interaction_matrix = interaction_matrix @ mod2act.T
        
    if config.closed:
        interaction_matrix *= -1.0
            
    if not silent:
        print("Computing laplcian...")
            
    if config.laplacian:
        inffuncdx = np.zeros((config.n_xsubapertures, config.n_ysubapertures, n_act), dtype=float)
        inffuncdy = np.zeros((config.n_xsubapertures, config.n_ysubapertures, n_act), dtype=float)
        laplacian = np.zeros((config.n_xsubapertures, config.n_ysubapertures, n_act), dtype=float)


        for i in range(n_act):
            inffuncdx[:, :, i] = interaction_matrix[0:n_centroids//2, i].reshape(config.n_xsubapertures, config.n_xsubapertures, order='F')
            inffuncdy[:, :, i] = interaction_matrix[n_centroids//2:n_centroids, i].reshape(config.n_xsubapertures, config.n_xsubapertures, order='F')
   
            laplacian[:, :, i] = (
                deriv2D(inffuncdy[:, :, i], y=True) +
                deriv2D(inffuncdx[:, :, i], x=True)
            )
            
    return freq_pos, fft_centroids_real, actuator_fft_magnitude, freq_mask, interaction_matrix, laplacian