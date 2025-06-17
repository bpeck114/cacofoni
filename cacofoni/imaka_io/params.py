# FILE: params.py 
# Still needs to be tested and made more efficient

from dataclasses import dataclass, field
import numpy as np
from typing import List
from math import pow

# Need to go back through and understand what is going on
def initimakaparmstruct(fparm):
    """
    Initializes a nested parameter structure from
    the imaka parameter file. 
    
    Inputs:
    -------
    fparm : str
          Path to parameter file. 
    
    Outputs:
    --------
    structure : dict 
              A dictionary representing the full 
              imaka parameter structure. 
    
    """
    
    # may not need to be a bool anymore
    success, nsub_vals = getiparm(fparm, "sys_parm.nsub")
    NSUB = int(nsub_vals[0]) if success else 0
    NSUBTOTAL = NSUB * NSUB
    
    sucsess, nact_vals = getiparm(fparm, "sys_parm.nact")
    NACT = int(nact_vals[0]) if success else 0
    
    sucsess, nwfs_vals = getiparm(fparm, "sys_parm.nwfs")
    NWFS = int(nwfs_vals[0]) if success else 0
    
    success, npixx_vals = getiparm(fparm, "wfscam_parm.npixx")
    # assumes line format like "0    120" 
    npixx = [int(line.split()[1]) for line in npixx_vals]
    
    if sum(npixx) != NWFS * npixx[0]:
        raise ValueError("Error in NPIXX")
    NPIXX = npixx[0] 
    NPIXPERSUB = [val / NSUB for val in npixx]
    
    NWFSMAX = 5
    NPIXPERFRAME = pow(NSUB * NPIXPERSUB[0], 2)
    
    # Guide star
    gs_parm = {
        "name": "",
        "mag": 0.0,
        "dRA": 0.0,
        "dDec": 0.0
    }
    
    gs_array = [gs_parm.copy() for _ in range(NWFSMAX)]
    

    return


    # Guide star parameters (gs_parm)
    gs_parm = {
        "name": "",
        "mag": 0.0,
        "dRA": 0.0,
        "dDec": 0.0
    }
    gs_array = [gs_parm.copy() for _ in range(NWFSMAX)]

    # Target parameters (target_parm)
    target_parm = {
        "name": "",
        "tRA": 0.0,
        "tDec": 0.0,
        "fRA": 0.0,
        "fDec": 0.0,
        "nGS": 0,
        "gs": gs_array
    }

    # WFS camera parameters (wfscam_parm)
    wfscam_parm = {
        "camera_sn": 0,
        "npixx": 0,
        "npixy": 0,
        "x0": 0,
        "y0": 0,
        "texp": 0.0,
        "temp": 0.0,
        "emgain": 0,
        "skyname": "",
        "flatname": ""
    }
    wfscam_array = [wfscam_parm.copy() for _ in range(NWFSMAX)]

    # WFS parameters (wfs_parm)
    wfs_parm = {
        "pixelweights": [0.0] * int(NPIXPERSUB[0]),
        "pixelthreshold": 0,
        "xsub": [0.0] * NSUBTOTAL,
        "ysub": [0.0] * NSUBTOTAL,
        "centroidoffsetname": ""
    }
    wfs_array = [wfs_parm.copy() for _ in range(NWFSMAX)]

    # DM parameters (dm_parm)
    dm_parm = {
        "xact": [0.0] * NACT,
        "yact": [0.0] * NACT,
        "voltmax": 0.0
    }

    # Loop parameters (loop_parm)
    loop_parm = {
        "niter": 0,
        "nave": 0,
        "gainp": 0.0,
        "gaini": 0.0,
        "gaind": 0.0,
        "imatname": "",
        "cmatname": [""] * 8,
        "cmat2name": ""
    }

    # System parameters (sys_parm)
    sys_parm = {
        "nwfs": NWFS,
        "nact": NACT,
        "nsub": NSUB,
        "ncb": 0,
        "sim_dm": 0,
        "sim_wfscams": 0
    }

    # Final nested parameter dictionary
    imakaparm = {
        "sys": sys_parm,
        "target": target_parm,
        "wfscam": wfscam_array,
        "wfs": wfs_array,
        "dm": dm_parm,
        "loop": loop_parm
    }

    return imakaparm


def getiparm(fparm, 
             keyword):
    """
    Searchs through the imaka parameter file for certain keywords.
    Returns a list of strings for value associated with keyword.
    
    Inputs:
    -------
    fparm : str
          Path to imaka parameter file.
          Should be the path to imakaparm.txt
    
    keyword : str
            Keyword to search for. 
            Examples: 'sys_parm.nwfs', 'target_parm.gs'

    Outputs:
    --------
    success : bool # May not need this bool feature anymore (come back)
            True if the keyword was found and returns a meaningful value.
            False otherwise
    
    values : list of str
           The values found on lines matching the keyword.
           Keyword is stripped. 
    
    Raises:
    -------
    ValueError
        If the keyword is not in the list of allowed keys. 
    
    """
    
    # 0) Setup
    # List of allowed keywords from imakaparm.txt
    allowed_keywords = {
        "sys_parm.nwfs", "sys_parm.nact", "sys_parm.nsub", "sys_parm.ncb",
        "sys_parm.ncbmax", "sys_parm.ncbskip", "sys_parm.cb_autoname",
        "sys_parm.sim_dm", "sys_parm.sim_wfscams",
        "target_parm.name", "target_parm.tRA", "target_parm.tDec",
        "target_parm.fRA", "target_parm.fDec", "target_parm.nGS", "target_parm.gs",
        "wfscam_parm.camera_sn", "wfscam_parm.npixx", "wfscam_parm.npixy",
        "wfscam_parm.x0", "wfscam_parm.y0", "wfscam_parm.texp", "wfscam_parm.temp",
        "wfscam_parm.emgain", "wfscam_parm.skyname", "wfscam_parm.flatname",
        "wfscam_parm.simcamfile",
        "wfs_parm.pixelweights", "wfs_parm.pixelthreshold",
        "wfs_parm.xsub", "wfs_parm.ysub"
    }
    
    # Raises error if given incorrect keyword. 
    if keyword not in allowed_keywords:
        raise ValueError(f"Keyword '{keyword}' is not in the allowed list for getiparm.")
    
    keyword_len = 30 # default (need to come back and make this more efficient)
    
    # 1) Tries to open the file 
    try:
        with open(fparm, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f'Cannot find {fparm}.')
        
    matches = []
    
    for line in lines:
        if keyword in line:
            parts = line.strip()
            if parts.startswith(keyword):
                matches.append(parts)
   
    values = [line[keyword_len:].strip() for line in matches]
    return True, values



"""
; 
; NAME:  initimakaparmstruct
; DESCRIPTION:  sets up an IDL structure for imaka parameters
; HISTORY:  
;  2015-06-25 - v0.25  - in sync with v0.25 of imaka.h
;+-----------------------------------------------------------------------------
FUNCTION initimakaparmstruct, fname

    	IF ( getiparm(fname, 'nsub', str) ) THEN NSUB = (FIX(str))[0]
	NSUBTOTAL = NSUB*NSUB
    	IF ( getiparm(fname, 'nact', str) ) THEN NACT = (FIX(str))[0]
    	IF ( getiparm(fname, 'nwfs', str) ) THEN NWFS = (FIX(str))[0]
	npixx = INTARR(NWFS)
	rc = getiparm(fname, 'npixx', str) 
    	FOR i=0,NWFS-1 DO npixx[i] =  (FIX(STRSPLIT(str[i],/EXTRACT)))[1]
	IF ( TOTAL(npixx) NE NWFS*npixx[0] ) THEN PRINT, 'Error in NPIXX'
	NPIXX = npixx[0]
    	NPIXPERSUB = FLOAT(npixx)/nsub
    
    	;; MAYBE GET THESE FROM lp_parm READ...
    	NWFSMAX = 5
    	NPIXPERFRAME = (NSUB*NPIXPERSUB)^2.

    ;; FIRST WE SETUP STRUCTURES THAT MATCH THE imaka C CODE

	;typedef struct gs_parm {
	;	char  name[NAMEMAXCHAR];
	;	float mag;
	;	float dRA;			// Offset in arcseconds to GS
	;	float dDec;			// (dRA, dDec)
	;} gs_parm;
	gs_parm = {gsparms, name:"", mag:0.0, dRA:0.0, dDec:0.0}
    	gs = REPLICATE({gsparms}, NWFSMAX) 

	;typedef struct target_parm {
	;	char name[NAMEMAXCHAR];  // Name of target
	;	float tRA;			// RA, Dec of center target
	;	float tDec;				
	;	float fRA;			// RA,Dec of center of field
	;	float fDec;
	;	int nGS;			     // number of GSs
	;	gs_parm	gs[NWFSMAX];	// structure for each GS
	;} target_parm;
	target_parm = {targetparms, name:"", tRA:0., tDec:0., fRA:0., fDec:0., nGS:0, gs: gs}

	;typedef struct wfscam_parm {
	;	int camera_sn;		     // Serial number of WFS camera
	;	int npixx;			// camera width of ROI
	;	int npixy;  			// camera height of ROI
	;	int x0;				// camera left of ROI
	;	int y0;				// camera bottom of ROI
	;	float texp;			// camera exposure time of WFS camera (msec)
	;	float temp;			// camera set temperature of WFS camera
	;	int emgain;			// camera numeric EM gain (not true EM gain)
	;	char skyname[NAMEMAXCHAR];    // name of sky frame
	;	char flatname[NAMEMAXCHAR];   // name of flat calibration file
	;} wfscam_parm;
	wfscam_parm = {wfscamparms, camera_sn: 0, npixx: 0, npixy: 0, x0: 0, y0: 0, $
					texp:0.0, temp:0.0, emgain: 0, skyname:"", flatname:"" }
    	wfscam = REPLICATE({wfscamparms}, NWFSMAX) 

	;typedef struct wfs_parm {
	;     float pixelweights[NPIXPERSUB];   // pixel centroid weights
	;     int pixelthreshold;;            // pixel intensity threshold
	;	float xsub[NSUB*NSUB];		// pixel for the lower-left corner of subap
	;	float ysub[NSUB*NSUB];		// subap LL pixel locations 
	;	char centroidoffsetname[NAMEMAXCHAR];	// name of centroid offset
	;} wfs_parm;
	wfs_parm = {wfsparms, pixelweights: FLTARR(NPIXPERSUB), pixelthreshold: 0., $
				xsub: FLTARR(NSUB*NSUB), ysub: FLTARR(NSUB*NSUB), centroidoffsetname: ""}
    	wfs = REPLICATE({wfsparms}, NWFSMAX) 

	;typedef struct dm_parm {
	;	float xact[NACT];    // x,y positions of each actuator in pupil space
	;	float yact[NACT];
	;} dm_parm;
	dm_parm = {dmparms, xact: FLTARR(NACT), yact: FLTARR(NACT), voltmax:0.0}

	;typedef struct loop_parm {
	;	unsigned long niter;          // maximum number of loop cycles to run
     ;	float gaini;     	// servo integrator gain (normally =1 for no leaky)
	;	float gainp;		// servo proportional gain (gainp*error)
	;	float gaind;		// servo derivative gain (gaind*(derror/dt) Not used.
	;	char imatname[NAMEMAXCHAR];
	;	char cmatname[NAMEMAXCHAR];
	;	char cmat2name[NAMEMAXCHAR];
	;} loop_parm;
	loop_parm = {lpparms, niter: 0UL, nave: 0UL, gainp: 0.0, gaini: 0.0, gaind:0.0, imatname: "", cmatname:STRARR(8), cmat2name:""}

	;typedef struct sys_parm {
	;	int nwfs;			// This is the actual number of WFSs being used (not 
	;					// to be confused with NWFSMAX)...maybe this is just nGS...
	;					// This should be set to the nGS when a target is loaded...
	;} sys_parm;
	sys_parm = {sysparms, nwfs:0, nact:0, nsub:0, ncb:0, sim_dm:0, sim_wfscams:0 }

	parm = {imakaparm, sys:sys_parm, target:target_parm, wfscam:wfscam, wfs:wfs, dm:dm_parm, loop: loop_parm}
RETURN, parm
END


; NAME: getiparm, filename, keyword, svalue
; DESCRIPTION: returns string corresponding to a keywork in an imaka parameter TEXT file.  
;   This is a helper routine for irptxt.
; OUTPUTS:  sets svalue to the parameter values (as a string with the keyword cropped off)
;           RETURNS an error code = TRUE (no errors) = FALSE (error empty string, didn't find keyword)
; HISTORY:
;+-----------------------------------------------------------------------------
FUNCTION getiparm, fname, keyword, svalue
	KEYWORD_LEN = 30
	
	SPAWN, 'grep  "'+keyword+'" '+fname, str
	IF ( STRLEN(str[0]) LT KEYWORD_LEN ) THEN svalue='' ELSE BEGIN
		;; Do case of more than one occurence of this keyword (e.g. target_parm.gs)
		ni = N_ELEMENTS(str)
		svalue = STRARR(ni)
		FOR i=0,ni-1 DO BEGIN
			svalue[i] = STRMID(str[i],KEYWORD_LEN-1, STRLEN(str[i])-KEYWORD_LEN+1)
		ENDFOR
	ENDELSE

RETURN, STRLEN(svalue[0]) NE 0
END


"""
