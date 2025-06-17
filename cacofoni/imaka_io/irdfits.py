# FILE: irdfits.py
# Still needs to be tested 

# Import packages
from astropy.io import fits
import numpy as np
import os 
from dataclasses import dataclass, field
from cacofoni.config import CacophonyConfig

def irdfits(fname, 
            fparm=None, 
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
    
        return
  

'''

    """

    Parameters
    ----------
    fname : str
        Path to the FITS file.
    fparm : str
        Path to the parameter file used for structure setup.
    silent : bool
        Suppress warnings if True.
    exten : list of int (length 8)
        Flags indicating which FITS extensions to read.

    Returns
    -------
    list of dict
        List of telemetry data structures, one per timestep.
    """
    from .initimakadatastruct import initimakadatastruct  # assumed implemented
    from .sxpar import sxpar  # helper function to extract header keywords
    
    if exten is None:
        exten = [1] * 8
    elif isinstance(exten, int):
        tmp = [0] * 8
        tmp[0] = 1
        tmp[exten] = 1
        exten = tmp

    hdul = fits.open(fname)
    
    # Extension 0: loop data
    d0 = hdul[0].data
    h0 = hdul[0].header
    ntimes = h0['NAXIS2']
    data = initimakadatastruct(ntimes=ntimes, fparm=fparm)
    
    for i in range(ntimes):
        data[i]['loop']['state'] = d0[0][i]
        data[i]['loop']['cntr'] = d0[1][i]
    
    nwfs = h0.get('NWFS', 0)

    # Extension 1: wfscam rawpixels
    if exten[1]:
        try:
            d1 = hdul[1].data
            h1 = hdul[1].header
            if not any("No wfscam" in card for card in h1.values()):
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['rawpixels'] = d1[:, :, i, t]
                        data[t]['wfscam'][i]['timestamp'] = h1.get(f'TSTAMPA{i}', 0)
                        data[t]['wfscam'][i]['fieldcount'] = 0
                        data[t]['wfscam'][i]['tsample'] = h1.get(f'TSAMPLE{i}', 0)
            else:
                if not silent: print("No wfs camera data saved")
        except Exception as e:
            if not silent: print("Extension 1 missing:", e)

    # Extension 2: wfscam processed pixels
    if exten[2]:
        try:
            d2 = hdul[2].data
            h2 = hdul[2].header
            if not any("No wfscam" in card for card in h2.values()):
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['pixels'] = d2[:, :, i, t]
                        data[t]['wfscam'][i]['timestamp'] = h2.get(f'TSTAMPA{i}', 0)
                        data[t]['wfscam'][i]['fieldcount'] = 0
                        data[t]['wfscam'][i]['tsample'] = h2.get(f'TSAMPLE{i}', 0)
            else:
                if not silent: print("No wfs camera (processed pixels) data saved")
        except Exception as e:
            if not silent: print("Extension 2 missing:", e)

    # Extension 3: wfs centroids
    if exten[3]:
        try:
            d3 = hdul[3].data
            h3 = hdul[3].header
            if not any("No wfs data" in card for card in h3.values()):
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfs'][i]['centroids'] = d3[:, i, t]
            else:
                if not silent: print("No wfs data saved")
        except Exception as e:
            if not silent: print("Extension 3 missing:", e)

    # Extension 4: dm voltages
    if exten[4]:
        try:
            d4 = hdul[4].data
            h4 = hdul[4].header
            if not any("No dm data" in card for card in h4.values()):
                for t in range(ntimes):
                    data[t]['dm']['deltav'] = d4[:, 0, t]
                    data[t]['dm']['voltages'] = d4[:, 1, t]
            else:
                if not silent: print("No dm data saved")
        except Exception as e:
            if not silent: print("Extension 4 missing:", e)

    # Extension 5: average wfscam pixels
    if exten[5]:
        try:
            d5 = hdul[5].data
            h5 = hdul[5].header
            if not any("No average wfscam" in card for card in h5.values()):
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['avepixels'] = d5[:, :, i]
            else:
                if not silent: print("No average wfs camera data saved")
        except Exception as e:
            if not silent: print("Extension 5 missing:", e)

    # Extension 6: average wfs centroids
    if exten[6]:
        try:
            d6 = hdul[6].data
            h6 = hdul[6].header
            if not any("No average wfs" in card for card in h6.values()):
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfs'][i]['avecentroids'] = d6[:, i]
            else:
                if not silent: print("No average wfs data saved")
        except Exception as e:
            if not silent: print("Extension 6 missing:", e)

    # Extension 7: average dm voltages
    if exten[7]:
        try:
            d7 = hdul[7].data
            h7 = hdul[7].header
            if not any("No average dm" in card for card in h7.values()):
                for t in range(ntimes):
                    data[t]['dm']['avevoltages'] = d7
            else:
                if not silent: print("No average dm data saved")
        except Exception as e:
            if not silent: print("Extension 7 missing:", e)

    hdul.close()
    return data
'''
    
            

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


