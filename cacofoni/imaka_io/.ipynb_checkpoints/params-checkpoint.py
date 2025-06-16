# FILE: params.py 
# Still needs to be tested 

from dataclasses import dataclass, field
import numpy as np
from typing import List

def initimakaparmstruct():
    """
    
    
    Inputs:
    -------
    
    Optional Inputs:
    ----------------
    
    Outputs:
    --------
    
    """

    return


def getiparm(fname, keyword):
    """
    Search through the parameter file for certain keywords.
    Returns the remaining content of those lines. 
    
    Inputs:
    -------
    
    Optional Inputs:
    ----------------
    
    Outputs:
    --------
    
    """
    
    




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

"""
