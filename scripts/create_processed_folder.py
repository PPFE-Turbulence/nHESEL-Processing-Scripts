# Imports
import os
import sys
from datetime import date
import shutil
import warnings

def create_processed_folder(inputName, overwrite = False):
    """
    Creates a processing folder with BOUT-inp within the processing directory defined in the paths.py file, using a .inp filename as input.
    The filename will be used as the name for the folder, and the file itself will be copied into the folder as BOUT.inp, ready to run a nHESEL simulation.


    Parameters
    ----------
    inputName : str
        Name of the BOUT input file found in the input folder
    overwrite : Bool. Optional.
        Bool indicating whether the data folder should be overwritten, if it exists. The default is False.

    Raises
    ------
    Exception
        In case the folder already exists, and overwrite is set to false, an exception is raised.

    Returns
    -------
    None.

    """
    
    # Returns the current local date
    today = date.today()
    dd = today.strftime("%d")
    mm = today.strftime("%m")
    yy = today.strftime("%y")
    
    # Import paths file
    pathPath = f"{os.environ['HOME']}/{'nHESEL-Processing-Scripts'}"
    sys.path.append(pathPath)
    import paths
    
    # Ensure folder does not exist already
    exists = os.path.exists(path=f"{paths.processedPath}/{inputName}")
    if exists == True and overwrite == False:
        raise Exception(f"{paths.processedPath}/{inputName} already exists, and overwrite is set to {overwrite}. Choose another input file, or set overwrite to True. Note setting overwrite = True will delete all data from the folder with the same name as the input file.")
    elif exists == True and overwrite == True:
        warnings.warn(f"The folder {paths.processedPath}/{inputName} exists and overwrite is set to {overwrite}. The folder and its files will be deleted.")
        shutil.rmtree(f"{paths.processedPath}/{inputName}")
            
    # Generate new processed folder with input filename
    os.mkdir(path=f"{paths.processedPath}/{inputName}")
    
    # Generate file saving the time the folder was generated
    metaFile = open(f"{paths.processedPath}/{inputName}/meta.txt","w+")
    metaFile.write(f"This folder was generated on {dd}-{mm}-{yy} (dd-mm-yy)")
    
    # For good measure
    metaFile.close()

    return 0
