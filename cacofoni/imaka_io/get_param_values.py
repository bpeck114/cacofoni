# FILE: read_param_values.py 
# Helper functions to read imaka-style parameter file and get values

# Set of allowed keywords in the parameter file
# Helps protect against typos and unsupported entries when called 
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

def get_param_values(fparam,
                     keyword,
                     which_column=None):
    """
    Returns processed values for a given parameter keyword.
    
    Inputs:
    -------
    fparam : str
           Path to imaka-style parameter file
    
    keyword : str
            The keyword to search for at the start of each line.
            Example: 'sys_parm.nact'
            
    Optional Inputs:
    ----------------
    which_column : int or None
                 If None (default), splits each matching line into parts
                    (as a list of strings).
                 
                 If an integer is provided, returns only that specific column
                    from each line.
                 
                 Note: Columns are 0-indexed. 
            
    Outputs:
    --------
    result : list
        If which_column is None: return list of lists, each split by whitespace.
        If which_column is given: returns a list of individual string values 
           from that column.
           
    Raises:
    -------
    ValueError
         If the requested column index is out of range for any line.
    """
    
    # Get the trimmed line(s) for the keyword
    raw_lines = get_raw_param_lines(fparam, keyword)
    result = []
    
    for line in raw_lines:
        # Split the line into tokens with whitespace
        parts = line.split()
        
        if which_column is None:
            # Return all tokens on the line
            result.append(parts)
            
        elif 0 <= which_column < len(parts):
            # Return only the token at the specified column index 
            result.append(parts[which_column])
            
        else:
            # Raise an error if wrong format
            raise ValueError(f"Column ({which_column}) out of range for {parts}.")
    
    return result

def get_raw_param_lines(fparam, keyword):
    """
    Finds all the lines in the parameter file that start with the given keyword.
    Returns the portion of the line after the keyword. 
    
    Inputs:
    -------
    fparam : str
           Path to imaka-style parameter file.
           Example: 'imakaparm.txt'
           
    keyword : str
            The keyword to search for at the start of each line.
            Example: 'sys_parm.nact'
            
    Outputs:
    --------
    list of str
         A list of strings containing the part of the line
         after the keyword.
         
         For example, if the line is: "wfscam_parm.npixx     0    120"
         It will return: "0    120"
         
    Raises:
    -------
    ValueError
         If the keyword is not in predefined `allowed_keywords' set.
    """
    
    # Checks for invalid keywords (typos)
    if keyword not in allowed_keywords:
        raise ValueError(f"Keyword '{keyword}' is not allowed.")

    values = []
    
    # Open the parameter file and read line-by-line
    with open(fparam, 'r') as file:
        for line in file:
            # Strip the leading or trailing whitespace
            # Check if line starts with the keyword
            if line.strip().startswith(keyword):
                # Split the line into 2 parts: keyword and everything else
                parts = line.strip().split(maxsplit=1)
                
                # Only keep the 'everything after' part if it exists
                if len(parts) == 2:
                    values.append(parts[1].strip())

    return values



                    
                    
'''

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
'''