# FILE: initimakadatastruct.py
# Under construction

# Import packages
import numpy as np
from cacofoni.imaka_io.get_param_values import get_param_values
from cacofoni.config import CacofoniConfig
import copy

def initimakadatastruct(fparam,
                        ntimes):
    '''
    Setting up an empty data structure
    to eventually hold telemetry data.
    
    Inputs:
    -------
    fparam : str
           Path to parameter file.
           
    ntimes : int
           Number of time steps. 
           
    Outputs:
    --------
    data : list of dict
         A list of 'ntimes' dictionaries,
         each representing a telemetry frame.
    '''
    
    
    # Contains hard-coded information on max number of wfs
    config = CacofoniConfig()
    nwfs_max = config.nwfs_max
    
    # Get the basic system parameters: subaperatures, actuator, number of WFS
    # get_param_values produces a list of lists, we want the first value inside
    nsub = int(get_param_values(fparam, "sys_parm.nsub")[0][0]) 
    nact = int(get_param_values(fparam, "sys_parm.nact")[0][0])
    
    # Assuming that the WFS is square (equal pixels in x- and y-direction)
    npixx = int(get_param_values(fparam, "wfscam_parm.npixx", which_column=1)[0])
    
    # Compute derived quantities
    nsub_total = nsub * nsub # Total subapertures (assumed square)
    npix_total = int(npixx ** 2)
    npix_per_sub = npixx / nsub # How many pixels wide is each subaperture
    
    print(f"WFS: {npixx}x{npixx} ({npix_total}) total pixels.")
    print(f"WFS: {nsub}x{nsub} ({nsub_total}) total subapertures.")
    print(f"WFS: {int(npix_per_sub)} pixels per subaperture.")
    
    
    # Defining an empty frame (data)
    # A = loop data (once per time step)
    # B = dm data (per actuator)
    # C = wfs data (per subaperture)
    # D = wfs cam data
    
    # A) 
    loop_data = {
        "state": 0, # 0 = idle, 1 = running
        "cntr": 0, # counter
    }
    
    # B)
    dm_data = {
        "deltav": [0.0] * nact,
        "voltages": [0.0] * nact,
        "avevoltages": [0.0] * nact,
    }
    
    # C)
    # Holds x and y data
    wfs_data = []
    for _ in range(nwfs_max):
        wfs_data.append({
            "raw_centroids": [0.0] * (2 * nsub_total),
            "centroids": [0.0] * (2 * nsub_total),
            "avecentroids": [0.0] * (2 * nsub_total),
        })
        
    # D)
    wfscam_data = []
    for _ in range(nwfs_max):
        wfscam_data.append({
            "timestamp": 0, # time when image was taken
            "fieldcount": 0, # frame number 
            "tsample": 0.0, # exposure time
            "rawpixels": np.zeros((npixx, npixx), dtype=float), # raw pixel image 
            "pixels": np.zeros((npixx, npixx), dtype=float), # processed image
            "avepixels": np.zeros((npixx, npixx), dtype=float), # average image
        })
        
    # One full frame template
    frame = {
        "loop": loop_data, # Holds idle/running + counter
        "dm": dm_data, # Holds deltav + voltages + avevoltages
        "wfs": wfs_data, # Holds raw centroids + centroids + avecentroids
        "wfscam": wfscam_data,
    }
    
    # Replicate across all of the time steps
    data = [copy.copy(frame) for _ in range(ntimes)]

    
    return data 


'''
; 
; NAME:  initimakadatastruct
; DESCRIPTION:  sets up an IDL structure for imaka data
; HISTORY:  
;  2015-06-25 - v1.0 
;+-----------------------------------------------------------------------------
FUNCTION initimakadatastruct, ntimes=NTIMES, fparm=FPARM
    
    ;; MAYBE GET THESE FROM lp_parm READ..
    IF ( getiparm(fparm, 'nsub', str) ) THEN NSUB = (FIX(str))[0]
    NSUBTOTAL = NSUB*NSUB
    IF ( getiparm(fparm, 'nact', str) ) THEN NACT = (FIX(str))[0]
    IF ( getiparm(fparm, 'nwfs', str) ) THEN NWFS = (FIX(str))[0]
    npixx = INTARR(NWFS)
    rc = getiparm(fparm, 'npixx', str) 
    FOR i=0,NWFS-1 DO npixx[i] =  (FIX(STRSPLIT(str[i],/EXTRACT)))[1]
    IF ( TOTAL(npixx) NE NWFS*npixx[0] ) THEN PRINT, 'Error in NPIXX'
    NPIXX = npixx[0]

    NPIXPERSUB = FLOAT(npixx)/nsub
    NPIXPERFRAME = (NSUB*NPIXPERSUB)^2.
    ;IF ( getiparm(fparm, 'ncb', str) ) THEN NTIMES = (FIX(str))[0]
    NWFSMAX = 5

    ;IF ( N_ELEMENTS(NTIMES) EQ 0 ) THEN NTIMES = 10 ; NCBMAXLENGTH.

    ;; FIRST WE SETUP STRUCTURES THAT MATCH THE imaka C CODE
    ;typedef struct wfscam_data {		     // Data for individual WFS Camera
    ;	long timestamp;
    ;	long	fieldcount;		          // index of buffer data
    ;	uint16_t	pixels[NPIXPERFRAME];		// processed pixels
    ;} wfscam_data
;    wfscam_data = {wfscams, timestamp:LONG64(0), fieldcount: LONG64(0), $
    wfscam_data = {wfscams, timestamp:LONG64(0), fieldcount: LONG64(0), tsample:0.0, $
                            rawpixels: UINTARR(NSUB*NPIXPERSUB,NSUB*NPIXPERSUB), $
                            pixels: FLTARR(NSUB*NPIXPERSUB,NSUB*NPIXPERSUB), $
					  avepixels: FLTARR(NSUB*NPIXPERSUB,NSUB*NPIXPERSUB) }

    ;typedef struct wfs_data {		// Data for individual WFS
    ;	float	 raw_centroids[2*NSUBTOTAL];	// raw centroids
    ;	float centroids[2*NSUBTOTAL];		// offset subtracted centroids
    ;} wfs_data
    wfs_data = {wfss, raw_centroids: FLTARR(2*NSUBTOTAL),  $
                    centroids: FLTARR(2*NSUBTOTAL), $
				avecentroids: FLTARR(2*NSUBTOTAL) }
    
    ;typedef struct dm_data {			// Data for individual DM
    ;	float deltav[NACT];	// delta voltages
    ;	float voltages[NACT];		// integrated voltages
    ;} dm_data;
    dm_data = {dms, deltav: FLTARR(NACT), voltages: FLTARR(NACT), avevoltages: FLTARR(NACT)}

    ;typedef struct loop_state_data {   // Loop state data
    ; int state;     // idle, running, etc (should replace gLoopThreadRunning)
    ; unsigned long cntr;     // increments for each loop
    ;} loop_state_data;
    lp_data = {lps, state: FIX(0), cntr: ULONG(0)}

    ;typedef struct lp_data {
    ;	long loop_cntr;
    ;	wfscam_data wfscam[NWFSMAX];
    ;	wfs_data wfs[NWFSMAX];	
    ;	dm_data dm;
    ;} lp_data;
    wfscam = REPLICATE({wfscams}, NWFSMAX) 
    wfs    = REPLICATE({wfss}, NWFSMAX)

    datas = {imakadata, loop: lp_data, wfscam: wfscam, wfs: wfs, dm: dm_data}
    data = REPLICATE(datas, NTIMES)

RETURN, dataf
END
'''