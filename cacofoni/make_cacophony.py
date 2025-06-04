# FILE: make_cacophony.py

# Import packages
import numpy as np
from astropy.io import fits
from scipy.fft import fft
from scipy.signal import windows
from scipy.io.wavfile import write as write_wav
import importlib_resources as resources
import cacofoni.data

n_modes = 36 # Number of actuators on ASM

def deriv2D(ima, x, y):
    return

def make_cacofoni(filename, minfreq, maxfreq, fparm, silent, thresh, closed=closed, modal=modal, laplacian=1, n_modes=36):
    """
    Generates the interaction and control matrix from wavefront
    sensor telemetry data.
    
    Inputs:
    -------
    filename : str
             Location of telemetry file. 
            
    minfreq : float
            Minimum frequency. 
    
    maxfreq : float
            Maximum frequency, should be 3.6. Hz higher than minimum frequency. 
    
    closed : bool
             closed = closed loop. 
             open = open loop. 
    
    modal : bool
             Modal = modulates actuator by actuator.
             Zonal = modulates by Zernike modes (tiptilt, focus, etc). 
    
    fparm : str
    
    silent : bool
    
    thresh : float 
    
    laplacian : bool
              Takes the curvature of the phase. 
              
    Optional Inputs:
    ----------------
    n_modes : 36
             Number of actuators on ASM. 
    
    
    Outputs:
    --------
    imat:
    
    cmat:
    """
    
    if fparam = fparm or resources.files("cacofoni.data").joinpath("imakaparm.txt")
    else: 
        fparm = fparm 
    
    with fits.open(filename) as hdul:
        
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
