# FILE: config.py
# Contains default configuration values and path
# Note: The files are assumed to be in data/.

from dataclasses import dataclass, field

@dataclass
class CacofoniConfig:
    telemetry_filename: str = "aocb0090.fits"
    param_filename: str = "imakaparm.txt"
    mirror_modes_filename: str = "mm2a_norm.fits"
    minimum_frequency: float = 4.0
    maximum_frequency: float = 10.0
    
    num_actuators: int = 36
    nwfs_max: int = 5
    
    extension: list = None  

    def __post_init__(self):
        if self.extension is None:
            self.extension = [1, 0, 0, 1, 1, 1, 1, 1]
    