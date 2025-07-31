# Module import
import sys
import numpy as np
import os
import warnings

# Import paths file
pathPath = "/home/s194112"
sys.path.append(pathPath)
import paths

# Append scripts folder
sys.path.append(paths.scriptPath)

# Scripts imports
from dumpProcessing import dumpProcessing
from heselProcessing import heselProcessing
from maxMeans import maxMeansEnergyflux

def lcfsEnergyFlux(dataName, tindexRange = [100,-1], info=False, perpOrpara = "perp", eOri = None, seriesMode = False):
    ## Interpret perpOrpara
    # Map it
    if type(perpOrpara) == str:
        perpOrpara = perpOrpara.lower()
        perpOrparaDict = {"perp" : 0, "para" : 1}
        mode = perpOrparaDict[perpOrpara]

    # Check it
    if mode not in [0, 1]:
        raise Exception(f"perpOrpara is currently set to {perpOrpara}, which is an invalid value. Choose either perp or para, or their corresponding value 0 and 1.")
    
    ## Get the (mean) q_perp for the LCFS index
    # First get the LCFS index
    dumpObj = dumpProcessing(dataName)
    LCFS_index = dumpObj.get_lcfs_index(includeLCFS = True)

    # Also get the grid positions for the two indices closest to the actual LCFS position, as well as the actual position
    grid = dumpObj.get_grids(scrapeoff = False)[0][LCFS_index:LCFS_index+2]
    lcfsPos = dumpObj.get_lcfs()

    # Calculate the lengths from the lcfs position
    grid -= lcfsPos

    # Then make the distance to lcfs into weights
    weights = 1/abs(grid)

    # Then get the average energy flux for this index in one of two ways
    if info == True:
        if mode == 0:
            print(f"LCFS Mode: Getting q_perpe and q_perpi for {dataName}")
        if mode == 1:
            print(f"LCFS Mode: Getting q_parae and q_parai for {dataName}")
    
    # In case the data has already been processed
    prePath = f"{paths.processedPath}/{dataName}"
    if mode == 0:
        perpePath = f"{prePath}/q_perpe_averaged_z.npy"
        perpiPath = f"{prePath}/q_perpi_averaged_z.npy"
    if mode == 1:
        perpePath = f"{prePath}/q_parae_averaged_z.npy"
        perpiPath = f"{prePath}/q_parai_averaged_z.npy"
    if os.path.isfile(perpePath) and os.path.isfile(perpiPath):
        # if info == True:
        #     print("Loading from processed data")
        
        # First load the data
        q_perpe = np.load(perpePath)
        q_perpi = np.load(perpiPath)

        # Take the index range
        if tindexRange[1] == -1:
            q_perpe = q_perpe[tindexRange[0]:, LCFS_index:LCFS_index+2]
            q_perpi = q_perpi[tindexRange[0]:, LCFS_index:LCFS_index+2]
        else:
            q_perpe = q_perpe[tindexRange[0]:tindexRange[1], LCFS_index:LCFS_index+2]
            q_perpi = q_perpi[tindexRange[0]:tindexRange[1], LCFS_index:LCFS_index+2]

        # Then take the mean over the t-axis, if we are not running a series
        if seriesMode == False:
            q_perpe = np.mean(q_perpe, axis=0)
            q_perpi = np.mean(q_perpi, axis=0)
        else:
            pass
    
    # In case it hasn't
    else:
        if info == True:
            print("Loading from RAW data")
        
        # First process the data (includes tindexRange)
        data = processingObj.tzaveragedVariables(["q_perpe","q_perpi"], tsliceSize = 100, tindexRange = tindexRange, xindexRange = [LCFS_index, LCFS_index+1], returnData = True, info = False)
        q_perpe, q_perpi = data

        # Then take the mean over the t-axis, if we are not running a series
        if seriesMode == False:
            q_perpe = np.mean(q_perpe, axis=0)
            q_perpi = np.mean(q_perpi, axis=0)
        else:
            pass

    # Two indices, so take the average of them, using the weights based on the inverse of how close to LCFS the two indices are. Remember to include axis of seriesMode
    if seriesMode == False: # [t,x,z] -> [x] -> []
        q_perpe = np.average(q_perpe, weights = weights)
        q_perpi = np.average(q_perpi, weights = weights)
    else: # [t,x,z] -> [t,x] -> [t]
        q_perpe = np.average(q_perpe, axis=1, weights = weights)
        q_perpi = np.average(q_perpi, axis=1, weights = weights)

    # And then add the values together
    if eOri != None:
        if eOri == "e":
            q_perp = q_perpe
        elif eOri == "i":
            q_perp = q_perpi
    else:
        q_perp = q_perpe + q_perpi

    return q_perp

def P_loss(dataName, tindexRange = [100,-1], info = False, mode = 1, extendedOutput = False, diffusionCorrection = 1.1, eOri = None, seriesMode = False):
    """
    Calculates the z- and t-averaged P_loss over the LCFS, in units of watt.

    The tindexRange indices which t-indices to utilize during calculations. Default if from 600 to -1.

    The possible modes are: LCFS, maxMeans, max (no case sensitivity)
     - LCFS: Take the approximate value of the energy fluxes at the LCFS (slow as hell if dumpProcessing has yet to extract values from this)
     - maxMeans [DEFAULT]: Take the mean value of q_perp in a range around the max value. See maxMeans function for further details.
     - max: Take the max value of q_perp within the interval

    Modes can also be given as '0', '1', and '2' corresponding to LCFS, maxMeans, and max, in order.

    extendedOutput set to True will result in outputting not only the P_loss value, but also the P_loss values for the ions and electrons, as well as the ballooning area. The output will basically change from simply P_loss to [P_loss, P_loss_e, P_loss_i, A_b].

    As the actual LCFS position does not fit exactly with the grid, a weighted average of the two closest grid points are used,
    where the weights are based upon the absolute inverse of the distance between the gridpoints and the actual LCFS position in the simulation.

    The formula used is P_loss = (q_perpe + q_perpi)*A_b, where q_perpe and q_perpi is the average outward energy-flux
    across the LCFS in [J/m^2*s], and A_b is the total ballooning area in [m^2], given as A_b = 2*pi*(Rmajor+Rminor)*2*pi*Rminor/6.
    Here Rmajor is the major radius of the tokamak, and Rminor the minor radius.

    seriesMode allows one to ge tthe calculated P_loss for each timestep, meaning only using a z-average instead of tz-average

    """
    # First convert modes to numbers, if need be
    if type(mode) == str:
        modeMap = {"lcfs" : 0, "maxmeans" : 1, "max" : 2}
        mode = modeMap[mode.lower()]

    # Now check that the mode is valid
    if mode not in [0,1,2]:
        raise Exception(f"{mode} is not a valid mode. Choose either 0, 1, 2, or their corresponding strings.")
    
    # Simple function to get the total effect over the LCFS
    ## Start by defining the required objects
    dumpObj = dumpProcessing(dataName = dataName)
    processingObj = heselProcessing(dataName = dataName)
    
    ## First calculate the ballooning area
    # Get the minor and major radii
    Rmajor = dumpObj.get_option_value("Rmajor")
    Rminor = dumpObj.get_option_value("Rminor")

    # Then calculate the ballooning area
    A_b = 2*np.pi*(Rmajor+Rminor)*2*np.pi*Rminor/6

    ### Now get the values for the LCFS (0) mode
    if mode == 0:
        # Call lcfsEnergyFlux
        q_perp = lcfsEnergyFlux(dataName, tindexRange = tindexRange, info = info, perpOrpara = "perp", eOri = eOri, seriesMode = seriesMode)

    ### Now for mode 1 & 2 (maxMeans & max)
    if (mode == 1) or (mode == 2):
        # First get the max and maxMean values
        maxVal, maxMean, toleranceBroken = maxMeansEnergyflux(dataName = dataName, eOri = eOri, seriesMode = seriesMode, tindexRange = tindexRange)

        # Check the tolerance
        if toleranceBroken == True:
            warnings.warn(f"The tolerance has been broken for {dataName}. Continuing script.")

        # Then get the q_perp value for the P_loss calculation
        if mode == 1:
            if info == True:
                print(f"maxMeans Mode: Getting q_perp for {dataName}")
            q_perp = maxMean
        if mode == 2:
            if info == True:
                print(f"max Mode: Getting q_perp for {dataName}")
            q_perp = maxVal
        

    ## Finally, add and multiply the quantities together to get the mean loss across the ~60 degree area
    P_loss = q_perp*A_b

    # Fix for special case where data is too short to gain data from
    if seriesMode == False:
        if (np.isnan(P_loss)):
            P_loss = 0
            q_perp = 0
            q_perpi = 0
            q_perpe = 0
            warnings.warn("P_loss had value nan, indicating the lower bound on the tindexRange was above the amount of data present in the simulation run. Returning P_loss and q_perp values as 0.")

    # extendedOutput
    if extendedOutput == True:
        ## In case we have mode 1 or 2, we need to get the q_perpe and q_perpi values first
        if (mode == 1) or (mode == 2):
            # Get the max and maxMean values for eelectrons and ions
            maxVale, maxMeane, toleranceBrokene = maxMeansEnergyflux(dataName = dataName, electronsOrIons = "e", seriesMode = seriesMode, tindexRange = tindexRange)
            maxVali, maxMeani, toleranceBrokeni = maxMeansEnergyflux(dataName = dataName, electronsOrIons = "i", seriesMode = seriesMode, tindexRange = tindexRange)

            # Take our the values depending on mode
            if mode == 1:
                q_perpe = maxMeane
                q_perpi = maxMeani
            if mode == 2:
                q_perpe = maxVale
                q_perpi = maxVali
            
        ## Now add them to the P_loss output
        P_loss = [P_loss, q_perpe*A_b, q_perpi*A_b, A_b]

        # And diffusioncorrection
        P_loss = [p*diffusionCorrection for p in P_loss]

    else:
        P_loss = P_loss*diffusionCorrection

    return P_loss
