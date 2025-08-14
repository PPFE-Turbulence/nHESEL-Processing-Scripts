# Module import
import sys
import datetime
import os

# Import paths file
pathPath = os.environ["HOME"]
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Function to process data
from heselProcessing import heselProcessing
from dumpProcessing import dumpProcessing
from create_processed_folder import create_processed_folder

def processFolder(dataName, plain = False, zMean = False, tMean = False, tzMean = False, zSum = False, tSum = False, tzSum = False, Init = False, scrapeoff = False):
    # Call processedfolder creation in case it does not exist
    if os.path.isdir(f"{paths.processedPath}/{dataName}") == False:
        create_processed_folder(inputName = dataName, overwrite = False)
    
    # Create the processingObj required to load the data
    processingObj = heselProcessing(dataName = dataName)

    # Define variables to load
    variables = ["n","Te","Ti","pe","pi","phi","u","Gamma_perp","Gamma_para","q_perpe","q_perpi","q_parae","q_parai"]

    # Do dumpProcessing pre-processing (part of the folder processing - speeds up EVERYTHING after this where dumpProcessing's built-in methods are required)
    dumpObj = dumpProcessing(dataName = dataName)
    dumpObj.preprocess_data()
    
    # Get LCFS index
    LCFS_index = dumpObj.get_lcfs_index()

    # Create the meta-file for which data extraction parameters were included
    f = open(f"{paths.processedPath}/{dataName}/meta.txt","a+")
    metaString = f"""
    =========================================================================================
    The following parameters were extracted on {datetime.datetime.utcnow()}
    The dataName was: {dataName}
    Calculate plain data: {plain}
    Calculate z-meaned data: {zMean}
    Calculate tz-meaned data: {tzMean}
    Calculate z-summed data: {zSum}
    Calculate tz-summed data: {tzSum}
    Evaluate initial conditions of above data: {Init}
    Evaluate for scrape-off layer of above data (including Init): {scrapeoff}

    The amount of successful nouts are: {len(dumpObj.get_times())}
    =========================================================================================\n\n
    """
    f.write(metaString)
    f.close()
    
    ### Load plain data ###
    if plain == True:
        print("Loading plain data")
        processingObj.plainVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.plainVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.plainVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.plainVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")

    ### Load z averaged data ###
    if zMean == True:
        print("Loading z averaged data")
        processingObj.zaveragedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.zaveragedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.zaveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.zaveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")

    ### Load t averaged data ###
    if tMean == True:
        print("Loading t averaged data")
        processingObj.taveragedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.taveragedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.taveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.taveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")
    
    ### Load tz averaged data ###
    if tzMean == True:
        print("Loading tz averaged data")
        processingObj.tzaveragedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.tzaveragedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.tzaveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.tzaveragedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")

    ### Load z summed data ###
    if zSum == True:
        print("Loading z summed data")
        processingObj.zsummedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.zsummedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.zsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.zsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")
    
    ### Load t summed data ###
    if tSum == True:
        print("Loading t summed data")
        processingObj.tsummedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.tsummedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.tsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.tsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")
    
    ### Load tz summed data ###
    if tzSum == True:
        print("Loading tz summed data")
        processingObj.tzsummedVariables(varNames = variables, returnData = False, saveData = True, info = False)

        if Init == True:
            processingObj.tzsummedVariables(varNames = variables, tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "Init")

            if scrapeoff == True:
                processingObj.tzsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], tindexRange = [0,0], returnData = False, saveData = True, info = False, customName = "scrapeoff_Init")
        
        if scrapeoff == True:
            processingObj.tzsummedVariables(varNames = variables, xindexRange = [LCFS_index, -1], returnData = False, saveData = True, info = False, customName = "scrapeoff")

    return 0