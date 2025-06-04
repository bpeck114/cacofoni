# FILE: config.py

from dataclasses import dataclass
from pathlib import Path
import importlib_resources
from cacofoni import data

@dataclass
class CacophonyConfig:
    param_filename: str = "imakaparm.txt"
    mirror_modes_filename: str = "mm2a_norm.fits"
    num_actuators: int = 36
    num_centroids = 288
    grid_shape: tuple = (12, 12)
    sampling_rate_hz: float = 996.0
    output_wavfile: str = "cacophony.wav"
    
    # Points path for deafult files to data/
    def resolve_data_path(self, filename: str):
        return str(importlib_resources.files(data) / filename)
    
    @property
    def fparam_path(self):
        return self.resolve_data_path(self.fparam_filename)
    
    @property
    def mirror_modes_path(self):
        return self.resolve_data_path(self.mirror_modes_filename)