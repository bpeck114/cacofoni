# FILE: initimakadatastruct.py
# Under construction

# Import packages
import numpy as np
from cacofoni.imaka_io.get_param_values import get_param_values
from cacofoni.config import CacofoniConfig

def initimakadatastruct(fparam,
                        ntimes):
    '''
    '''
    
    config=CacofoniConfig()
    
    # Store core parameters from the file
    nsub = int(get_param_values(fparam, "sys_parm.nsub")[0][0])
    nact = int(get_param_values(fparam, "sys_parm.nact")[0][0])
    nwfs = int(get_param_values(fparam, "sys_parm.nwfs")[0][0])
    
    # Read the camera pixel size for each WFS
    npixx_list = [int(x) for x in get_param_values(fparam, "wfscam_parm.npixx", which_column=1)]
    
    if sum(npixx_list) != nwfs * npixx_list[0]:
        print("Warning: inconsistent npixx across WFSs")
        
    npixx = npixx_list[0]
    
    # Compute derived quantities
    nsub_total = nsub * nsub
    npix_per_sub = npixx / nsub
    npix_per_frame = int((nsub * npix_per_sub) ** 2)
    
    # Defining an empty frame
    
    # Loop state info
    loop_data = {
        "state": 0,
        "cntr": 0,
    }
    
    # DM info (voltages and deltas)
    dm_data = {
        "deltav": [0.0] * nact,
        "voltages": [0.0] * nact,
        "avevoltages": [0.0] * nact,
    }
    
    # WFS sensor data
    wfs_data = []
    for _ in range(config.nwfs_max):
        wfs_data.append({
            "raw_centroids": [0.0] * (2 * nsub_total),
            "centroids": [0.0] * (2 * nsub_total),
            "avecentroids": [0.0] * (2 * nsub_total),
        })
        
    # WFS camera data
    wfscam_data = []
    for _ in range(config.nwfs_max):
        wfscam_data.append({
            "timestamp": 0,
            "fieldcount": 0,
            "tsample": 0.0,
            "rawpixels": [[0] * int(npix_per_sub * nsub) for _ in range(int(npix_per_sub * nsub))],
            "pixels":     [[0.0] * int(npix_per_sub * nsub) for _ in range(int(npix_per_sub * nsub))],
            "avepixels":  [[0.0] * int(npix_per_sub * nsub) for _ in range(int(npix_per_sub * nsub))],
        })
        
    # One full frame template
    frame = {
        "loop": loop_data,
        "wfscam": wfscam_data,
        "wfs": wfs_data,
        "dm": dm_data,
    }
    
    # Replicate across all of the time steps
    data = [frame.copy() for _ in range(ntimes)]
    
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