# FILE: irdfits.py
# Still needs to be tested 

# Import packages
import numpy as np
from astropy.io import fits
from cacofoni.imaka_io.initimakadatastruct import initimakadatastruct
from cacofoni.imaka_io.get_param_values import get_param_values
from cacofoni.config import CacofoniConfig


def irdfits(ftele, 
            fparam,
            exten=None):
    """
    Populate the structured array 
    from initimakadatastruct.
    
    Extensions:
      0    = Loop control (always used)
      1    = Raw wfscam pixels
      2    = Processed wfscam pixels 
      3    = WFS centroids 
      4    = DM voltages
      5    = Average wfscam pixels
      6    = Average wfs centroids
      7    = Average dm voltages
    """
    
    # Default extension loading
    if exten is None:
        exten = [1] * 8
    elif isinstance(exten, int):
        tmp = [0] * 8
        tmp[exten] = 1
        exten = tmp
        
    # Load basic structure from extension 0
    with fits.open(ftele) as hdul:
        print("\nLoading telemetry FITS file...")
        h0 = hdul[0].header
        d0 = hdul[0].data  # raw loop control data
        ntimes = h0.get('NAXIS2', d0.shape[0]) 
        nwfs = int(h0['NWFS'])
        
        print(f"Number of time steps: {ntimes}")
        print(f"Number of WFS: {nwfs}")
        
        print("\nSetting up empty data structure...")
        cb = initimakadatastruct(ntimes=ntimes, fparam=fparam)
        
        print("\nChecking extensions...")
        print("\nRunning E0: ")
        d0 = fits.getdata(ftele, ext=0).T  # shape (5, ntimes)

        loop_state_array = np.int16(d0[0, :])  # equivalent to FIX(d[0,*])
        loop_cntr_array  = np.uint32(d0[1, :]) # equivalent to ULONG(d[1,*])

        # Fill into cb structure
        for i in range(ntimes):
            cb[i]['loop']['state'] = loop_state_array[i]
            cb[i]['loop']['cntr']  = loop_cntr_array[i]
        
        
        if exten[1]:
            print("Running E1: Raw wfscam data")
        else:
            print("Skipping E1: Raw wfscam data")
        
        
        if exten[2]:
            print("Running E2: Processed wfscam data")
        else:
            print("Skipping E2: Processed wfscam data")
        
        if exten[3]:
            print("Running E3: WFS centroids")
            d3 = hdul[3].data 
            # d3 shape: (nsamp, nwfs, 2*nsub) = (27000, 1, 288)
            
            # Transpose to match idl logic (for now)
            d3T = np.transpose(d3, (2, 1, 0)) 
            # d3T shape : (2*nsub, nwfs, nsamp) = (288, 1, 27000)
            
            for i in range(nwfs):
                cb[0]['wfs'][i]['centroids'] = d3T[:, i, :] 
                # shape: (2*nsub, nsamp) = (288, 27000)
       
        '''
        if exten[4]:
            print("Running E4: DM data")
            d4 = hdul[4].data  # shape: (NACT, 2, ntimes)
            for t in range(ntimes):
                cb[t]['dm']['deltav']   = d4[:, 0, t]
                cb[t]['dm']['voltages'] = d4[:, 1, t]
            
        if exten[5]:
            print("Running E5: Average wfscam")
            d5 = hdul[5].data  # (nwfs, y, x)
            for i in range(nwfs):
                cb[0]['wfscam'][i]['avepixels'] = d5[i, :, :].T  # (x, y)
                
        if exten[6]:
            print("Running E6: Average WFS centroids")
            d6 = hdul[6].data  # shape: (nwfs, 2*nsub)
            for i in range(nwfs):
                cb[0]['wfs'][i]['avecentroids'] = d6[i, :]

        
        if exten[7]:
            print("Running E7: Average DM voltages")
            d7 = hdul[7].data
            cb[0]['dm']['avevoltages'] = d7
        '''
        
    return cb


    
            

"""
; NAME: irdfits, filename
; DESCRIPTION: reads a "standard" imaka data FITS file into an imaka IDL data structure
; HISTORY:
;	2015-06-25 - v1.0 - reads FITS file from dataclient v.?
;+-----------------------------------------------------------------------------
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


