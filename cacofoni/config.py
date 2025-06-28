# FILE: config.py
# Contains default configuration values and path
# Note: The files are assumed to be in data/.

# Import packages
from dataclasses import dataclass, field

@dataclass
class CacofoniConfig:
    telemetry_filename: str = "aocb0090.fits"
    param_filename: str = "imakaparm.txt"
    modal_filename: str = "mm2a_norm.fits"
    
    minimum_freq_hz: float = 4.0
    maximum_freq_hz: float = 10.0
    
    n_actuators: int = 36
    n_xsubapertures: int = 12
    n_ysubapertures: int = 12
    
    closed: bool = False
    modal: bool = False
    thresh: bool = None
    apply_hanning: bool = False
    laplacian: bool = True
    silent: bool = False