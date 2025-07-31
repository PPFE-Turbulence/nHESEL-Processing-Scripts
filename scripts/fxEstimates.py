# Imports
import numpy as np
import sys
import warnings

### Locals imports
# Import paths file
pathPath = "/home/s194112/"
sys.path.append(pathPath)
import paths

# Import job submission script
sys.path.append(paths.scriptPath)

# Scripts import
from FOLs import FOLs
from downstreamMap import downstreamMap
from dumpProcessing import dumpProcessing

def fxEstimates(dataName = None, eOri = None, constMode = False):
    """
    Function which calculates two flux line spreading estimates (fx), by mapping the PFOL and DFOL values downstream to the diverter, and taking the downstream/upstream fraction.

    The returned values are fx_PFOL, fx_DFOL.
    """

    if (constMode == False) and (dataName != None):
        # First get the PFOL and DFOL for the data
        PFOL, DFOL = FOLs(dataName, eOri = eOri)
        
        # Then map them both
        PFOLmapped, DFOLmapped = downstreamMap([PFOL, DFOL])
        
        # And now divide the downstream with the upstream PFOL and DFOL to get fx estimates for this data
        fx_PFOL = PFOLmapped/PFOL
        fx_DFOL = DFOLmapped/DFOL
    
        # Check to ensure NaN issues are taken care of
        if np.isnan(fx_PFOL) or np.isnan(fx_DFOL) or fx_PFOL == -np.inf or fx_DFOL == -np.inf:
            warnings.warn("Either fx_PFOL or fx_DFOL are NaN or -inf, indicating something is wrong with the loaded data. Returning 0 for fx_PFOL and fx_DFOL.")
            fx_PFOL = 0
            fx_DFOL = 0
    
        return fx_PFOL, fx_DFOL

    elif (constMode == True):
        # First define the value right between LCFS and wall
        # dumpObj = dumpProcessing(dataName = dataName)
        # wallPos = dumpObj.get_wall(scrapeoff = True)[1] # wall position relative to LCFS
        # midPos = wallPos/2 # Midway point

        # Check which constant (mean) value to use
        if eOri == None:
            pos = 0.019822226299382444 # Combined
        elif eOri == "e":
            pos = 0.013165301625718714 # Electrons
        elif eOri == "i":
            pos = 0.02567993412089612 # Ions

        # Then get how long this point is projected
        projectedPos = downstreamMap([pos])[0]

        # And now get the spreading
        fx_Pos = projectedPos/pos

        return fx_Pos
        
    elif (constMode == "DFOL"):
        # Define the position to use
        pos = 0.044838766120814645

        # Then get how long this point is projected
        projectedPos = downstreamMap([pos])[0]

        # And now get the spreading
        fx_Pos = projectedPos/pos

        return fx_Pos
        

    else:
        raise Exception("Something went wrong with your inputs. If using constmode, you can ignore the dataName input. Otherwise, it will be required.")
        