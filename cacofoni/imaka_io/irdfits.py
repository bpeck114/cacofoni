# FILE: irdfits.py
# Still needs to be tested 

# Import packages
from astropy.io import fits
import numpy as np
import os 
from cacofoni.config import CacophonyConfig

def irdfits(fname, 
            fparm=None, 
            silent=False, 
            exten=None):
    """
    Reads a imaka telemetry FITS file and returns
    a structured object with WFS and DM telemetry. 
    
    Inputs:
    -------
    fname : str
          Path to FITS file.
          Example: "/data/asm/20250326/ao/aocb0090.fits"
    
    Optional Inputs:
    ----------------
    fparm : str or None
          Path to parameter file.
          Default is located in: "cacofoni/data/imakaparm.txt"
          
    silent : bool
           If True, suppress warning messages. 
    
    exten : 
          Specifies which FITS extensions to load.
          - If None: loads all 8 extensions
          - If int: load only that exetension
          - If list: interpreted as 8-element on/off mask
          - Note: Extension '0' is read every time, primary
          header/data unit
    
    Outputs:
    ---------
    """
    
    # 1) If no parameter file path given, goes with the default
    if fparm is None:
        fparm = CacophonyConfig.fparam_path
        print(f"Parameter File: {CacophonyConfig.fparam_path}")
        
    # 2) Interpret the `exten' parameter
    if exten is None:
        # If no extension list is given, defaults to all 8 on
        exten = [1] * 8
    elif isinstance(exten, int):
        # If a single extension number is given, make a full list
        # with only that one on. 
        tmp_exten = [0] * 8
        tmp_exten[0] = 1  # always read extension 0
        tmp_exten[exten] = 1
        exten = tmp_exten
    elif isinstance(exten, list) and len(exten) !=8:
        raise ValueError("exten must be an 8-element list, an integer, or None")
        
    # 3) Read in extension 0
    with fits.open(fname) as hdul:
        ext0 = hdul[0]
        d0 = ext0.data # shape: (2, NTIMES)
        h0 = ext0.header 
        
        ntimes = h0.get('NAXIS2')
        if ntimes is None:
            raise ValueError("Missing NAXIS2 in primary header — cannot determine number of time steps.")
            
        # 4) Initalize the full telemetry data structure 
        
             
        return
    

def initimakadatastruct():
            

"""
FUNCTION irdfits, fname, fparm=fparm, silent=SILENT, exten=exten
	IF ( N_ELEMENTS(fparm) EQ 0 ) THEN fparm='/tmp/imakaparm.txt'
     IF ( N_ELEMENTS(exten) EQ 0 ) THEN exten=[1,1,1,1,1,1,1,1]
     IF ( N_ELEMENTS(exten) EQ 1 ) THEN BEGIN
        ext = INTARR(8)
        ext[0] = 1
        ext[exten] = 1
        exten = ext
     ENDIF
     
	;; exten=0 (loop data)
	;; This extension is guaranteed to exist
	d = mrdfits(fname, 0, h, silent=SILENT)
     	NTIMES = SXPAR(h,'NAXIS2')
	data = initimakadatastruct(ntimes=ntimes,fparm=fparm)
	
	;; Read all the extensions with headers

	data[0:NTIMES-1].loop.state = REFORM(FIX(d[0,*]))
	data[0:NTIMES-1].loop.cntr = REFORM(ULONG(d[1,*]))
	nwfs = SXPAR(h,'NWFS')

	;; exten=1 (wfscam data - raw pixels)
	IF ( exten[1] ) THEN BEGIN
       d = mrdfits(fname, 1, h, /UNSIGNED, silent=SILENT)
	  IF ( WHERE(STRPOS(h,"No wfscam (rawpixels) data saved") NE -1) EQ -1 ) THEN BEGIN
            	FOR i=0,nwfs-1 DO data.wfscam[i].rawpixels = REFORM(d[*,*,i,*])
		FOR i=0,nwfs-1 DO BEGIN
    		;;timestamp:LONG64(0), fieldcount: LONG64(0), pixels: UINTARR(NPIXPERFRAME)
			data.wfscam[i].timestamp = SXPAR(h,'TSTAMPA'+STRING(i,FORMAT='(I1)'))
			data.wfscam[i].fieldcount = 0
               		data.wfscam[i].tsample = SXPAR(h,'TSAMPLE'+STRING(i,FORMAT='(I1)'))
		ENDFOR
	  ENDIF ELSE PRINT, 'No wfs camera data saved'
     ENDIF

	;; exten=2 (wfscam data - processed pixels)
	IF ( exten[2] ) THEN BEGIN
       d = mrdfits(fname, 2, h, silent=SILENT)
	  IF ( WHERE(STRPOS(h,"No wfscam (processed pixels) data saved") NE -1) EQ -1 ) THEN BEGIN
            	FOR i=0,nwfs-1 DO data.wfscam[i].pixels = REFORM(d[*,*,i,*])
		FOR i=0,nwfs-1 DO BEGIN
    		;;timestamp:LONG64(0), fieldcount: LONG64(0), pixels: UINTARR(NPIXPERFRAME)
			data.wfscam[i].timestamp = SXPAR(h,'TSTAMPA'+STRING(i,FORMAT='(I1)'))
			data.wfscam[i].fieldcount = 0
               		data.wfscam[i].tsample = SXPAR(h,'TSAMPLE'+STRING(i,FORMAT='(I1)'))
		ENDFOR
	  ENDIF ELSE PRINT, 'No wfs camera (processed pixels) data saved'
     ENDIF
     
	;; exten=3 (wfs data)
	IF ( exten[3] ) THEN BEGIN
       d = mrdfits(fname,3,h, silent=SILENT)
	  IF ( WHERE(STRPOS(h,"No wfs data saved") NE -1) EQ -1 ) THEN BEGIN
            	FOR i=0,nwfs-1 DO data.wfs[i].centroids = REFORM(d[*,i,*])
		;; raw_centroids: FLTARR(2*NSUBTOTAL),  centroids: FLTARR(2*NSUBTOTAL)
	  ENDIF ELSE PRINT, 'No wfs data saved'
     ENDIF

	;; exten=4 (dm data)
	IF ( exten[4] ) THEN BEGIN
       d = mrdfits(fname,4,h, silent=SILENT)
	  IF ( WHERE(STRPOS(h,"No dm data saved") NE -1) EQ -1 ) THEN BEGIN
     		;; dm_data = {dms, deltav: FLTARR(NACT), voltages: FLTARR(NACT)}
		data.dm.deltav = REFORM(d[*,0,*])
		data.dm.voltages = REFORM(d[*,1,*])	
	  ENDIF ELSE PRINT, 'No dm data saved'
	ENDIF

	;; exten=5 (ave wfscam data)
	IF ( exten[5] ) THEN BEGIN 
       d = mrdfits(fname,5,h, silent=SILENT, status=STATUS,/UNSIGNED)
	  IF ( status EQ 0 ) THEN BEGIN
		IF ( WHERE(STRPOS(h,"No average wfscam data saved") NE -1) EQ -1 ) THEN $
            		FOR i=0,nwfs-1 DO data.wfscam[i].avepixels = REFORM(d[*,*,i])
            		;data.wfscam.avepixels = d
	  ENDIF ELSE PRINT, 'No average wfs camera data saved'
	ENDIF

	;; exten=6 (wfs data)
	IF ( exten[6] ) THEN BEGIN
       d = mrdfits(fname,6,h, silent=SILENT, status=STATUS)
	  IF ( status EQ 0 ) THEN BEGIN
		IF ( WHERE(STRPOS(h,"No average wfs data saved") NE -1) EQ -1 ) THEN $
            		FOR i=0,nwfs-1 DO data.wfs[i].avecentroids = REFORM(d[*,i])
            	;data.wfs.avecentroids = d
	  ENDIF ELSE PRINT, 'No average wfs data saved'
     ENDIF
     
	;; exten=7 (dm data)
	IF ( exten[7] ) THEN BEGIN
       d = mrdfits(fname,7,h, silent=SILENT, status=STATUS)
	  IF ( status EQ 0 ) THEN BEGIN
		IF ( WHERE(STRPOS(h,"No average dm data saved") NE -1) EQ -1 ) THEN $
			data.dm.avevoltages = d
	  ENDIF ELSE PRINT, 'No average dm data saved'
     ENDIF
       
RETURN, data[0:NTIMES-1]

END
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

; NAME: irdtxt, filename
; DESCRIPTION: reads a "standard" imaka parameter TEXT file into an imaka IDL data structure
; HISTORY:
;	2015-06-25 - v1.0 - reads FITS file from dataclient v.?
;+-----------------------------------------------------------------------------
"""