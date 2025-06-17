# FILE: irdfits.py
# Still needs to be tested 

# Import packages
import numpy as np
from astropy.io import fits
from cacofoni.imaka_io.initimakadatastruct import initimakadatastruct
from cacofoni.config import CacofoniConfig


def irdfits(ftele, 
            fparam,
            exten=None):
    """
    """
    
    # Always use the first (0 index) extension
    if exten is None:
        exten = [1] * 8
    elif isinstance(exten, int):
        temp = [0] * 8
        temp[exten] = 1
        exten = temp
    
    # Extension 0: Loop state data (always exists)
    with fits.open(ftele, no_scale=True, ignore_end=True) as hdul:
        h0 = hdul[0].header
        d0 = hdul[0].data
        
        ntimes = h0['NAXIS2']
        nwfs = h0['NWFS']
        
        if d0.shape[0] == ntimes and d0.shape[1] == 5:
            d0 = d0.T  # Now d0.shape = (5, ntimes)
        
        data = initimakadatastruct(ntimes=ntimes, fparam=fparam)
        
        for t in range(ntimes):
            data[t]['loop']['state'] = int(d0[0][t])
            data[t]['loop']['cntr'] = int(d0[1][t])
            
        # Extension 1: Raw pixels from WFS cameras
        if exten[1] and len(hdul) > 1:
            print("E1: Loading raw pixels from WFS cameras.")
            if not any("No wfscam (rawpixels) data saved" in c for c in hdul[1].header.cards):
                d = hdul[1].data
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['rawpixels'] = d[:, :, i, t]
                        data[t]['wfscam'][i]['timestamp'] = h0.get(f'TSTAMPA{i}', 0)
                        data[t]['wfscam'][i]['tsample'] = h0.get(f'TSAMPLE{i}', 0.0)

        # Extension 2: Processed pixels 
        if exten[2] and len(hdul) > 2:
            print("E2: Loading processed pixels.")
            if not any("No wfscam (processed pixels) data saved" in c for c in hdul[2].header.cards):
                d = hdul[2].data
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['pixels'] = d[:, :, i, t]
                        data[t]['wfscam'][i]['timestamp'] = h0.get(f'TSTAMPA{i}', 0)
                        data[t]['wfscam'][i]['tsample'] = h0.get(f'TSAMPLE{i}', 0.0)
        else:
            print("SKIPPING. E2: Loading processed pixels.")
                        
        # Extension 3: Centroids
        if exten[3] and len(hdul) > 3:
            print("E3: Loading centroids.")
            if not any("No wfs data saved" in c for c in hdul[3].header.cards):
                d = hdul[3].data
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfs'][i]['centroids'] = d[:, i, t]
        else:
            print("SKIPPING. E3: Loading centroids.")
                        
        
        # Extension 4: DM voltages
        if exten[4] and len(hdul) > 4:
            print("E4: Loading DM voltages.")
            if not any("No dm data saved" in c for c in hdul[4].header.cards):
                d = hdul[4].data
                for t in range(ntimes):
                    data[t]['dm']['deltav'] = d[:, 0, t]
                    data[t]['dm']['voltages'] = d[:, 1, t]
        else:
            print("SKIPPING. E4: Loading DM voltages.")
                    
        # Extension 5: Average wfscam data
        if exten[5] and len(hdul) > 5:
            print("E5: Loading average WFS camera data.")
            if not any("No average wfscam data saved" in c for c in hdul[5].header.cards):
                d = hdul[5].data
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfscam'][i]['avepixels'] = d[:, :, i]
        else:
            print("SKIPPING. E5: Loading average WFS camera data.")
                        
        
        # Extension 6: Average WFS data
        if exten[6] and len(hdul) > 6:
            print("E6: Loading average WFS data.")
            if not any("No average wfs data saved" in c for c in hdul[6].header.cards):
                d = hdul[6].data
                for i in range(nwfs):
                    for t in range(ntimes):
                        data[t]['wfs'][i]['avecentroids'] = d[:, i]
        else:
            print("SKIPPING. E6: Loading average WFS data.")
                        
        # Extension 7: Average dm data
        if exten[7] and len(hdul) > 7:
            print("E7: Loading average DM data.")
            if not any("No average dm data saved" in c for c in hdul[7].header.cards):
                d = hdul[7].data
                for t in range(ntimes):
                    data[t]['dm']['avevoltages'] = d
                    
        else:
            print("SKIPPING. E7: Loading average DM data.")
    
    return data


    
            

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


