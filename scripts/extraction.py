# Module import
import sys

# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Load required scripts/libraries
from heselProcessing import heselProcessing

def extraction(dataName)
    # Create required object
    processingObj = heselProcessing(dataName = dataName)
    
    # Load and save data
    extractVars = ["n","Te","Ti","phi","u","Gamma","Qe","Qi", "pe", "pi"]
    processingObj.tzaveragedVariables(varNames = extractVars, returnData = False, saveData = True, info = False)
    processingObj.tzaveragedVariables(varNames = extractVars, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")
    print("Finished")
