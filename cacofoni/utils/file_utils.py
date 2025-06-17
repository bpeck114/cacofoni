# FILE: file_utils.py

from pathlib import Path
import importlib_resources
from cacofoni import data

def get_valid_path(user_path,
                   default_filename):
    """
    """
    
    if user_path is not None:
        path = Path(user_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"File not found at user-provided path: {path}")
        return path
    
    else:
        path = importlib_resources.files(data) / default_filename
        if not path.is_file():
            raise FileNotFoundError(f"File not found at default path inside package: {path}")
        return path
