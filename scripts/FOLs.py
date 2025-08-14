### Libraries imports
import sys
import numpy as np
import os
import warnings
from scipy.integrate import simpson

### Locals imports
# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Import scripts folder
sys.path.append(paths.scriptPath)

# Scripts imports
from dumpProcessing import dumpProcessing
from heselProcessing import heselProcessing

def FOLs(dataName, tindexRange = [75, -1], relativeToLCFS = True, info = False, eOri = None):
    """
    Function to calculate Power Fall-off Length and Particle/Density Fall-off Length, of the *AVERAGE* profiles, using processed or unprocessed data.

    In case the scrapeoff-layer processed data is available, it will be utilized. Otherwise, the required data will be loaded from the raw data, and saved for later use.

    Input:
     - dataName: str with name of data folder or equivalent processed data folder
     - relativeToLCFS: Whether the FOLs returned are given relative to the LCFS position (which means the actual PFOL and DFOL quantity), or not (which will then simply be the resulting PFOL and DFOL positions in the system). Turning this off effectively changes the outputs from Fall-off *Lengths* to Fall-off *Positions*
         - Default is true
     - info: Whether to print simple information about the data loading process, to see progress.

     Output:
     - PFOL: float, Power Fall-off Length
     - DFOL: float, Density Fall-off Length

    """
    ### Load the required data
    if info == True:
        print(f"Getting FOLs for {dataName}")
    ## Initialization
    # Basic definitions
    minIndex = tindexRange[0]
    
    # Start by defining paths
    prePath = f"{paths.processedPath}/{dataName}"
    paraePath = f"{prePath}/q_parae_averaged_z_scrapeoff.npy"
    paraiPath = f"{prePath}/q_parai_averaged_z_scrapeoff.npy"
    gammaparaPath = f"{prePath}/Gamma_para_averaged_z_scrapeoff.npy"

    # Then also define a heselProcessing object
    processingObj = heselProcessing(dataName = dataName)

    # And dumpProcessing for x grid, and LCFS coordinate
    dumpObj = dumpProcessing(dataName = dataName)

    ## Actual loading and/or processing of required data
    # Start by getting LCFS index
    LCFS = dumpObj.get_lcfs_index()
    
    # Check if already processed, or process if applicable
    
    # First for the energy flux
    # q parallel electrons first
    # if info == True:
    #     print("Getting q_parae")
    if os.path.isfile(paraePath):
        q_parae = np.load(paraePath)[tindexRange[0]:tindexRange[1],:]
    else:
        q_parae = processingObj.zaveragedVariables(["q_parae"], xindexRange = [LCFS,-1], tindexRange = tindexRange, returnData = True, saveData = True, customName = "scrapeoff")[0] # No reason not to save it for next time

    # Then q parallel ions
    # if info == True:
    #     print("Getting q_parai")
    if os.path.isfile(paraiPath):
        q_parai = np.load(paraiPath)[tindexRange[0]:tindexRange[1],:]
    else:
        q_parai = processingObj.zaveragedVariables(["q_parai"], xindexRange = [LCFS,-1], tindexRange = tindexRange, returnData = True, saveData = True, customName = "scrapeoff")[0]

    # Then combine the data, taking into account if we only need one of the two
    if eOri != None:
        if eOri == "e":
            q_para = q_parae
        elif eOri == "i":
            q_para = q_parai
    else:
        q_para = q_parae + q_parai

    # Now get the mean along the time-axis, to get the tz-average
    q_para = np.mean(q_para, axis=0)

    # Then for the particle flux
    # And the parallel particle flux
    # if info == True:
    #     print("Getting Gamma_para")
    if os.path.isfile(gammaparaPath):
        Gamma_para = np.load(gammaparaPath)[tindexRange[0]:tindexRange[1],:]
    else:
        Gamma_para = processingObj.zaveragedVariables(["Gamma_para"], xindexRange = [LCFS,-1], tindexRange = tindexRange, returnData = True, saveData = True, customName = "scrapeoff")[0]

    # And now get the time-mean
    Gamma_para = np.mean(Gamma_para, axis=0)

    ## Then calculate the fall-off lengths
    # Start by getting the x-scrapeoff grid
    xgrid = (dumpObj.get_grids(scrapeoff = True)[0])[1:] # [1:] to get rid of '0' coordinate
    
    # Now for the energy flux (Power fall-off length)
    PFOL = simpson(xgrid*q_para,xgrid)/simpson(q_para, xgrid)

    # Then for the particle flux (Density fall-off length)
    DFOL = simpson(xgrid*Gamma_para,xgrid)/simpson(Gamma_para, xgrid)

    # In case we do not wish to get the coordinates relative to the LCFS position
    if relativeToLCFS == False:
        LCFSpos = dumpObj.get_lcfs()[1]
        PFOL += LCFSpos
        DFOL += LCFSpos

    # Fix for special case where data is too short to gain data from
    if (np.isnan(PFOL)) or (np.isnan(PFOL)):
        PFOL = 0
        DFOL = 0
        warnings.warn("PFOL or DFOL had value nan, indicating the lower bound on the tindexRange was above the amount of data present in the simulation run. Returning PFOL and DFOL values as 0.")

    return PFOL, DFOL

def FOLs_index(dataName, relativeToLCFS = True, info = False):
    """
    Function to calculate Power Fall-off Length and Particle/Density Fall-off Length indix locations (so not actual PFOL and DFOL quantities, but resulting positions given in indices instead).

    In case the scrapeoff-layer processed data is available, it will be utilized. Otherwise, the required data will be loaded from the raw data, and saved for later use.

    Input:
     - dataName: str with name of data folder or equivalent processed data folder
     - relativeToLCFS: Whether the FOLs returned are given relative to the LCFS index, or not
         - Default is true
     - info: Whether to print simple information about the data loading process, to see progress.

     Output:
     - PFOL: int, Power Fall-off index position
     - DFOL: int, Density Fall-off index position

    """
    # First get the current positions (relative to LCFS)
    PFOL, DFOL = FOLs(dataName = dataName, relativeToLCFS = relativeToLCFS, info = info)

    # Then get a grid (with the 0 removed), and try to match the closest index to it
    # First get the grid
    dumpObj = dumpProcessing(dataName = dataName)
    xgrid = (dumpObj.get_grids(scrapeoff = relativeToLCFS)[0])[1:]

    # Then match the two, by looking for the minimum of the squared difference
    PFOL_index = np.argmin((xgrid - PFOL)**2)
    DFOL_index = np.argmin((xgrid - DFOL)**2)

    return PFOL_index, DFOL_index
