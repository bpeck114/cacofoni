# FILE: getiparm.py 
# Has been tested in python 

def getiparm(fparam, keyword):
    """
    Reads a text file and searches for lines
    starting with a given keyword.
    
    Returns whether any matching lines were
    found, and a list of the values found. 
    
    Inputs:
    -------
    fparam : str
          The path to the parameter file.
    keyword : str
          The keyword to search for.
          
    Outputs:
    --------
    success : bool
            True if at least one matching line was found.
            False otherwise.
    svalue : lst
           A list of strings, one for each matching line
           with the keyword cropped off. 
    """
    
    svalue = [] # Empty list will hold the results
    
    try:
        # Will try to open the file and read all the lines
        with open(fparam, 'r') as f:
            lines = f.readlines()
            
    except Exception as e:
        # If the file doesn't exist, returns an error
        print("Error reading file:", e)
        return False, []
    
    # Loops through each line in the file
    for line in lines:
        # Removes leading/trailing spaces and newlines
        line = line.strip()
        
        # Skips blank lines and comments
        if not line or line.startswith("#"):
            continue
            
        # Checks if the line starts with the given keyword
        if line.startswith(keyword):
            # First, removes keyword part from beginning of the line
            # Second, keeps only the values that follow it
            cropped = line[len(keyword):].strip() 
            
            if cropped: # Only keeps line is something is present after keyword
                svalue.append(cropped)
     
    # Will return True if the function found something with the list of values
    success = len(svalue) > 0 and len(svalue[0]) > 0    
    return (success, svalue)