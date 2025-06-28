# FILE: _load_and_log_telemetry.py

# Import packages 
import numpy as np
from astropy.io import fits

from cacofoni.calc_samp_freq import calc_samp_freq
from cacofoni.config import CacofoniConfig
from cacofoni.file_utils import check_default_path, check_user_path

def _load_and_log_telemetry(telemetry_filepath, silent):
    if not silent:
        print("Loading imaka telemetry file...\n")
        
    config = CacofoniConfig()
    telemetry_data = load_telemetry(telemetry_filepath, silent)
    
    clocktimes = telemetry_data['loop']['clocktime']   
    xcentroids = telemetry_data['wfs']['xcentroids']  
    ycentroids = telemetry_data['wfs']['ycentroids']
    deltav = telemetry_data['dm']['deltav']
    voltage = telemetry_data['dm']['voltage']
    
    n_steps = clocktimes.shape[0]
    n_act = deltav.shape[1]
    n_wfs = xcentroids.shape[1]
    n_xcentroids = xcentroids.shape[2]
    n_ycentroids = ycentroids.shape[2]
    n_centroids = n_xcentroids + n_ycentroids
    
    n_xsub = int(np.sqrt(n_xcentroids))
    n_ysub = int(np.sqrt(n_ycentroids))
    n_sub = n_xsub * n_ysub
    
    sampling_freq_hz, avg_time_step = calc_samp_freq(clocktimes)
    nyquist_freq_hz = sampling_freq_hz / 2
    estimated_duration_s = int(n_steps * avg_time_step)
        
    if not silent:
        print("Assumptions from telemetry file:")
        print("------------------------------------------------")
        print(f"Number of time steps           = {n_steps}")
        print(f"Average time interval (s)      = {avg_time_step:.4f}")
        print(f"Approximate runtime (s)        = {estimated_duration_s} ")
        print()
        print(f"Number of wavefront sensors    = {n_wfs}")
        print(f"Number of actuators            = {n_act}")
        print()    
        print(f"Number of x centroids          = {n_xcentroids}")
        print(f"Number of y centroids          = {n_ycentroids}")
        print(f"Total number of centroids      = {n_centroids}")
        print()
        print(f"Sampling frequency (Hz)        = {sampling_freq_hz:.2f}")
        print(f"Nyquist frequency (Hz)         = {nyquist_freq_hz:.2f}")
        print("------------------------------------------------\n")
        
    return xcentroids, ycentroids, deltav, voltage, n_steps, n_act, n_centroids, sampling_freq_hz, nyquist_freq_hz

def load_telemetry(telemetry_filepath, 
                   silent):
    """
    Simple loader for imaka telemetry FITS file.
    
    Parameters
    ----------
    ftele : str
        Path to the telemetry FITS file.

    silent : bool
        Suppress print statements.

    Returns
    -------
    data : dict
        Dictionary of numpy arrays with keys:
        'loop', 'centroids', 'dm'
    """

    telemetry_data = {}
    telemetry_data['loop'] = {}
    telemetry_data['wfs'] = {}
    telemetry_data['dm'] = {}
    
    # Note: Actual number of actuators from 
    # imaka telemetry file does not match 
    # size of telemetry array
    
    config = CacofoniConfig
    n_act = config.n_actuators
    n_xsub = config.n_xsubapertures 
    n_ysub = config.n_ysubapertures
    n_sub = n_xsub *  n_ysub

    with fits.open(telemetry_filepath) as hdul:

        if not silent: 
            print(f"Loading Extension 0: Loop state")
            print(f"------------------------------------------------")
        
        d0 = hdul[0].data # shape: (27000, 5) = (ntimes, parameter)
        # parameter = (0) counter, (1) state, (2) clocktime, (3) dTime, (4) WFScam time
        
        if not silent:
            print(f"Shape of Extension             = {d0.shape}")
        
        telemetry_data['loop']['clocktime'] = d0[:, 2] # shape: (27000,) = (ntimes, clocktime)

        if not silent: 
            for main_key in telemetry_data:
                for sub_key in telemetry_data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = '{sub_key}': {telemetry_data['loop']['clocktime'].shape}")
                    warn_if_zero("loop.clocktime", telemetry_data['loop']['clocktime'], silent) # checks if array is all zeroes

            print(f"------------------------------------------------\n")

            print(f"Loading Extension 3: Centroids")
            print(f"------------------------------------------------")
        
        d3 = hdul[3].data  # shape: (27000, 1, 288) = (ntimes, nwfs, ncentroids)
        # ncentroids = 144 x-slopes + 144 y-slopes
        # 144 = 12x12 n_sub
        
        if not silent: 
            print(f"Shape of Extension             = {d3.shape}")

        telemetry_data['wfs']['xcentroids'] = d3[:, :, :n_sub] # shape: (27000, 1, 144) = (ntimes, nwfs, xcentroids)
        telemetry_data['wfs']['ycentroids'] = d3[:, :, n_sub:] # shape: (27000, 1, 144) = (ntimes, nwfs, ycentroids)
        
        if not silent: 
            for main_key in list(telemetry_data.keys())[1:]:
                for sub_key in telemetry_data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = {telemetry_data[main_key][sub_key].shape}")
                    warn_if_zero("loop.clocktime", telemetry_data['loop']['clocktime'], silent) # checks if array is all zeroes
            print(f"------------------------------------------------\n")
            
            print(f"Loading Extension 4: DM Voltages")
            print(f"------------------------------------------------")
        
        d4 = hdul[4].data # shape: (27000, 2, 64) = (ntimes, voltage_mes, nact_driver)
        # voltage_mes = deltav, voltage
        # nact = 36 actuators with 36 padded zeroes 
        
        if not silent: 
            print(f"Shape of Extension             = {d4.shape}")
        
        telemetry_data['dm']['deltav'] = d4[:, 0, :n_act] # shape: (27000, 36) = (ntimes, deltav, nact)
        telemetry_data['dm']['voltage'] = d4[:, 1, :n_act] # shape: (27000, 36) = (ntimes, voltage, nact)
        
        if not silent: 
            for main_key in list(telemetry_data.keys())[2:]:
                for sub_key in telemetry_data[main_key]:
                    print(f"Loading                        = '{sub_key}' in '{main_key}'")
                    print(f"Key Shape                      = {telemetry_data[main_key][sub_key].shape}")
                    warn_if_zero("loop.clocktime", telemetry_data['loop']['clocktime'], silent) # checks if array is all zeroes
            print(f"------------------------------------------------\n")
        
    return telemetry_data

def warn_if_zero(name, 
                 arr, 
                 silent):
    
    if not silent: 
        if np.all(arr == 0):
            print(f"WARNING: The array '{name}' contains all zeroes.")
            
    return