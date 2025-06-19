# FILE: getiparm.py
# Used to be idl getiparm function

def get_param_values(filepath, 
                     keyword, 
                     which_column=None, 
                     cast_type=str):
    """
    Reads an imaka-style parameter file and gets the value(s) for a keyword.

    Parameters
    ----------
    filepath : str
        Path to the parameter file (e.g. 'imakaparm.txt')
    
    keyword : str
        Parameter name to search for (e.g. 'sys_parm.nact')
    
    which_column : int or None
        If None, returns all parts of the line after keyword.
        If int, only return that one column (0-indexed within the value part).
    
    cast_type : type
        Type to convert the value(s) to (e.g., int, float, str)

    Returns
    -------
    list
        List of values (as strings or casted to type).
    """
    
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

    if keyword not in allowed_keywords:
        raise ValueError(f"Keyword '{keyword}' is not allowed. Check for typos.")

    results = [] # Collects all matching results from keyword 

    # Open and read file line-by-line
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()  # Remove surronding whitespace 
            if line.startswith(keyword): # Check if line contains desired keyword
                parts = line.split() # Split line into keyword and value
                values = parts[1:]  # Skip the keyword

                if which_column is None:
                    # Return all values assigned to the keyword 
                    # Cast to the desired type 
                    casted = [cast_type(v) for v in values]
                    results.append(casted)
                else:
                    # Return only the specified column
                    if which_column < len(values):
                        results.append(cast_type(values[which_column]))
                    else:
                        # If the column doesn't exist on that line, raise an error
                        raise ValueError(f"Column {which_column} is out of range for line: {line}")
    
    if not results:
        raise ValueError(f"No values found for keyword: {keyword}")
    
    return results


def get_single_value(filepath, 
                     keyword, 
                     column=0, 
                     cast_type=int):
    """
    Shortcut to get a single value from a parameter file.

    Example:
    >>> get_single_value("imakaparm.txt", "sys_parm.nact")
    >>> 64
    """
    # Use the general function to get the desired values 
    values = get_param_values(filepath, keyword, which_column=column, cast_type=cast_type)
    
    return values[0] # Return the first match only
